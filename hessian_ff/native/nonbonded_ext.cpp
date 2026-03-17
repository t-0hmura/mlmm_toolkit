#include <torch/extension.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

// References:
// - hessian_ff/terms/nonbonded.py formulas in this repository.
// - OpenMM NonbondedForce conventions (Coulomb + LJ/HB style terms).
// - OpenMM CPU threading model/property (reference for thread-count behavior):
//   https://github.com/openmm/openmm/blob/master/platforms/cpu/src/CpuPlatform.cpp
//   https://github.com/openmm/openmm/blob/master/platforms/cpu/src/CpuPlatform.h
// Note:
// - The adaptive pair-count-based thread selection implemented below is NOT a
//   byte-for-byte OpenMM policy. It is a hessian_ff-specific optimization
//   inspired by OpenMM's CPU-threaded execution model.

namespace {

constexpr double kCoulomb = 332.0637132991921;

int64_t min_pairs_per_thread_threshold(int64_t requested) {
  if (requested > 0) {
    return requested;
  }
  // Tuned via benchmark sweep on bundled small+large data.
  // Current default prioritizes aggregate throughput across representative sizes.
  constexpr int64_t kDefault = 100000;
  const char* env = std::getenv("HESSIAN_FF_MIN_PAIRS_PER_THREAD");
  if (env == nullptr || env[0] == '\0') {
    return kDefault;
  }
  char* end = nullptr;
  const long long v = std::strtoll(env, &end, 10);
  if (end == env || v <= 0) {
    return kDefault;
  }
  return static_cast<int64_t>(v);
}

int effective_max_threads() {
#ifdef _OPENMP
  int max_threads = std::max(1, omp_get_max_threads());
  const int nprocs = std::max(1, omp_get_num_procs());

  const char* env = std::getenv("HESSIAN_FF_NONBONDED_MAX_THREADS");
  if (env != nullptr && env[0] != '\0') {
    char* end = nullptr;
    const long long v = std::strtoll(env, &end, 10);
    if (end != env && v > 0) {
      max_threads = static_cast<int>(v);
    }
  } else {
    // Prefer all OpenMP-visible processors by default. This avoids being
    // unintentionally capped by unrelated intra-op thread settings.
    max_threads = std::max(max_threads, nprocs);
  }
  return std::max(1, max_threads);
#else
  return 1;
#endif
}

std::tuple<torch::Tensor, torch::Tensor, torch::Tensor> pair_energy_force_aten(
    const torch::Tensor& coords,
    const torch::Tensor& charge,
    const torch::Tensor& atom_type,
    const torch::Tensor& lj_acoef,
    const torch::Tensor& lj_bcoef,
    const torch::Tensor& hb_acoef,
    const torch::Tensor& hb_bcoef,
    const torch::Tensor& nb_index,
    const torch::Tensor& ii,
    const torch::Tensor& jj,
    const c10::optional<torch::Tensor>& inv_scee_opt,
    const c10::optional<torch::Tensor>& inv_scnb_opt) {
  auto rij = coords.index_select(0, jj) - coords.index_select(0, ii);
  auto r2 = (rij * rij).sum(-1).clamp_min(1.0e-24);
  auto inv_r = torch::rsqrt(r2);
  auto inv_r2 = inv_r * inv_r;

  auto qq = charge.index_select(0, ii) * charge.index_select(0, jj);
  auto e_coul = qq * (double)kCoulomb * inv_r;
  auto fscale_coul = -(qq * (double)kCoulomb) * inv_r2 * inv_r;

  auto ti = atom_type.index_select(0, ii);
  auto tj = atom_type.index_select(0, jj);
  auto raw_idx = nb_index.index({ti, tj});

  auto inv_r4 = inv_r2 * inv_r2;
  auto inv_r6 = inv_r4 * inv_r2;
  auto inv_r8 = inv_r6 * inv_r2;
  auto inv_r10 = inv_r8 * inv_r2;
  auto inv_r12 = inv_r6 * inv_r6;
  auto inv_r14 = inv_r12 * inv_r2;

  auto e_lj = torch::zeros_like(inv_r);
  auto fscale_lj = torch::zeros_like(inv_r);

  auto lj_mask = raw_idx > 0;
  if (lj_mask.any().item<bool>()) {
    auto lj_idx = raw_idx.index({lj_mask}) - 1;
    auto a = lj_acoef.index_select(0, lj_idx);
    auto b = lj_bcoef.index_select(0, lj_idx);
    e_lj.index_put_({lj_mask}, a * inv_r12.index({lj_mask}) - b * inv_r6.index({lj_mask}));
    fscale_lj.index_put_(
        {lj_mask},
        -12.0 * a * inv_r14.index({lj_mask}) + 6.0 * b * inv_r8.index({lj_mask}));
  }

  auto hb_mask = raw_idx < 0;
  if (hb_mask.any().item<bool>() && hb_acoef.numel() > 0 && hb_bcoef.numel() > 0) {
    auto hb_idx_full = (-raw_idx.index({hb_mask})) - 1;
    auto valid = (hb_idx_full >= 0) & (hb_idx_full < hb_acoef.numel()) &
                 (hb_idx_full < hb_bcoef.numel());
    if (valid.any().item<bool>()) {
      auto hb_idx = hb_idx_full.index({valid});
      auto hb_pos = torch::nonzero(hb_mask).reshape(-1).index({valid});
      auto a = hb_acoef.index_select(0, hb_idx);
      auto b = hb_bcoef.index_select(0, hb_idx);
      e_lj.index_put_({hb_pos}, a * inv_r12.index({hb_pos}) - b * inv_r10.index({hb_pos}));
      fscale_lj.index_put_(
          {hb_pos},
          -12.0 * a * inv_r14.index({hb_pos}) + 10.0 * b * inv_r12.index({hb_pos}));
    }
  }

  if (inv_scee_opt.has_value()) {
    auto inv_scee = inv_scee_opt.value();
    e_coul = e_coul * inv_scee;
    fscale_coul = fscale_coul * inv_scee;
  }
  if (inv_scnb_opt.has_value()) {
    auto inv_scnb = inv_scnb_opt.value();
    e_lj = e_lj * inv_scnb;
    fscale_lj = fscale_lj * inv_scnb;
  }

  auto fij = (fscale_coul + fscale_lj).unsqueeze(-1) * rij;
  return {e_coul.sum(), e_lj.sum(), fij};
}

std::tuple<torch::Tensor, torch::Tensor, torch::Tensor> pair_energy_force_preparam_aten(
    const torch::Tensor& coords,
    const torch::Tensor& ii,
    const torch::Tensor& jj,
    const torch::Tensor& coul_coeff,
    const torch::Tensor& a12_coeff,
    const torch::Tensor& b6_coeff,
    const torch::Tensor& b10_coeff) {
  auto rij = coords.index_select(0, jj) - coords.index_select(0, ii);
  auto r2 = (rij * rij).sum(-1).clamp_min(1.0e-24);
  auto inv_r = torch::rsqrt(r2);
  auto inv_r2 = inv_r * inv_r;
  auto inv_r4 = inv_r2 * inv_r2;
  auto inv_r6 = inv_r4 * inv_r2;
  auto inv_r8 = inv_r6 * inv_r2;
  auto inv_r10 = inv_r8 * inv_r2;
  auto inv_r12 = inv_r6 * inv_r6;
  auto inv_r14 = inv_r12 * inv_r2;

  auto e_coul = coul_coeff * inv_r;
  auto e_lj = a12_coeff * inv_r12 - b6_coeff * inv_r6 - b10_coeff * inv_r10;
  auto fscale = -coul_coeff * inv_r2 * inv_r - 12.0 * a12_coeff * inv_r14 +
                6.0 * b6_coeff * inv_r8 + 10.0 * b10_coeff * inv_r12;
  auto fij = fscale.unsqueeze(-1) * rij;
  return {e_coul.sum(), e_lj.sum(), fij};
}

std::vector<torch::Tensor> nonbonded_energy_force_aten(
    const torch::Tensor& coords,
    const torch::Tensor& charge,
    const torch::Tensor& atom_type,
    const torch::Tensor& lj_acoef,
    const torch::Tensor& lj_bcoef,
    const torch::Tensor& hb_acoef,
    const torch::Tensor& hb_bcoef,
    const torch::Tensor& nb_index,
    const torch::Tensor& pair_i,
    const torch::Tensor& pair_j,
    const torch::Tensor& pair14_i,
    const torch::Tensor& pair14_j,
    const torch::Tensor& pair14_inv_scee,
    const torch::Tensor& pair14_inv_scnb,
    int64_t chunk_size) {
  TORCH_CHECK(chunk_size > 0, "chunk_size must be > 0");
  TORCH_CHECK(coords.dim() == 2 && coords.size(1) == 3, "coords must be [N,3]");
  TORCH_CHECK(coords.is_floating_point(), "coords must be floating tensor");
  TORCH_CHECK(charge.is_floating_point(), "charge must be floating tensor");
  TORCH_CHECK(
      coords.scalar_type() == charge.scalar_type(),
      "coords/charge dtype mismatch in nonbonded_energy_force");
  TORCH_CHECK(
      coords.scalar_type() == lj_acoef.scalar_type() && coords.scalar_type() == lj_bcoef.scalar_type(),
      "coords/LJ coeff dtype mismatch in nonbonded_energy_force");
  TORCH_CHECK(
      coords.scalar_type() == hb_acoef.scalar_type() && coords.scalar_type() == hb_bcoef.scalar_type(),
      "coords/HB coeff dtype mismatch in nonbonded_energy_force");
  TORCH_CHECK(atom_type.scalar_type() == torch::kInt64, "atom_type must be int64");
  TORCH_CHECK(nb_index.scalar_type() == torch::kInt64, "nb_index must be int64");
  TORCH_CHECK(pair_i.scalar_type() == torch::kInt64, "pair_i must be int64");
  TORCH_CHECK(pair_j.scalar_type() == torch::kInt64, "pair_j must be int64");
  TORCH_CHECK(pair14_i.scalar_type() == torch::kInt64, "pair14_i must be int64");
  TORCH_CHECK(pair14_j.scalar_type() == torch::kInt64, "pair14_j must be int64");
  TORCH_CHECK(
      pair_i.device() == coords.device() && pair_j.device() == coords.device() &&
          pair14_i.device() == coords.device() && pair14_j.device() == coords.device(),
      "pair index tensors must be on same device as coords");
  TORCH_CHECK(
      charge.device() == coords.device() && atom_type.device() == coords.device() &&
          lj_acoef.device() == coords.device() && lj_bcoef.device() == coords.device() &&
          hb_acoef.device() == coords.device() && hb_bcoef.device() == coords.device() &&
          nb_index.device() == coords.device() && pair14_inv_scee.device() == coords.device() &&
          pair14_inv_scnb.device() == coords.device(),
      "nonbonded parameter tensors must be on same device as coords");
  TORCH_CHECK(coords.size(0) == charge.numel(), "coords/charge atom count mismatch");
  TORCH_CHECK(atom_type.numel() == charge.numel(), "atom_type/charge atom count mismatch");
  TORCH_CHECK(pair_i.numel() == pair_j.numel(), "pair_i/pair_j size mismatch");
  TORCH_CHECK(pair14_i.numel() == pair14_j.numel(), "pair14_i/pair14_j size mismatch");
  TORCH_CHECK(
      pair14_i.numel() == pair14_inv_scee.numel() && pair14_i.numel() == pair14_inv_scnb.numel(),
      "1-4 pair/scale size mismatch");

  auto e_coul = torch::zeros({}, coords.options());
  auto e_lj = torch::zeros({}, coords.options());
  auto e_coul14 = torch::zeros({}, coords.options());
  auto e_lj14 = torch::zeros({}, coords.options());
  auto force = torch::zeros_like(coords);

  const int64_t n = pair_i.numel();
  for (int64_t start = 0; start < n; start += chunk_size) {
    int64_t end = std::min(start + chunk_size, n);
    auto ii = pair_i.slice(0, start, end);
    auto jj = pair_j.slice(0, start, end);
    auto out = pair_energy_force_aten(
        coords,
        charge,
        atom_type,
        lj_acoef,
        lj_bcoef,
        hb_acoef,
        hb_bcoef,
        nb_index,
        ii,
        jj,
        c10::nullopt,
        c10::nullopt);
    auto ce = std::get<0>(out);
    auto le = std::get<1>(out);
    auto fij = std::get<2>(out);
    e_coul = e_coul + ce;
    e_lj = e_lj + le;
    force.index_add_(0, ii, fij);
    force.index_add_(0, jj, -fij);
  }

  const int64_t n14 = pair14_i.numel();
  for (int64_t start = 0; start < n14; start += chunk_size) {
    int64_t end = std::min(start + chunk_size, n14);
    auto ii = pair14_i.slice(0, start, end);
    auto jj = pair14_j.slice(0, start, end);
    auto inv_scee = pair14_inv_scee.slice(0, start, end);
    auto inv_scnb = pair14_inv_scnb.slice(0, start, end);
    auto out = pair_energy_force_aten(
        coords,
        charge,
        atom_type,
        lj_acoef,
        lj_bcoef,
        hb_acoef,
        hb_bcoef,
        nb_index,
        ii,
        jj,
        inv_scee,
        inv_scnb);
    auto ce = std::get<0>(out);
    auto le = std::get<1>(out);
    auto fij = std::get<2>(out);
    e_coul14 = e_coul14 + ce;
    e_lj14 = e_lj14 + le;
    force.index_add_(0, ii, fij);
    force.index_add_(0, jj, -fij);
  }

  return {e_coul, e_lj, e_coul14, e_lj14, force};
}

std::vector<torch::Tensor> nonbonded_energy_force_preparam_aten(
    const torch::Tensor& coords,
    const torch::Tensor& pair_i,
    const torch::Tensor& pair_j,
    const torch::Tensor& pair_coul_coeff,
    const torch::Tensor& pair_a12_coeff,
    const torch::Tensor& pair_b6_coeff,
    const torch::Tensor& pair_b10_coeff,
    const torch::Tensor& pair14_i,
    const torch::Tensor& pair14_j,
    const torch::Tensor& pair14_coul_coeff,
    const torch::Tensor& pair14_a12_coeff,
    const torch::Tensor& pair14_b6_coeff,
    const torch::Tensor& pair14_b10_coeff,
    int64_t chunk_size) {
  TORCH_CHECK(chunk_size > 0, "chunk_size must be > 0");
  TORCH_CHECK(coords.dim() == 2 && coords.size(1) == 3, "coords must be [N,3]");
  TORCH_CHECK(coords.is_floating_point(), "coords must be floating tensor");
  TORCH_CHECK(pair_i.scalar_type() == torch::kInt64, "pair_i must be int64");
  TORCH_CHECK(pair_j.scalar_type() == torch::kInt64, "pair_j must be int64");
  TORCH_CHECK(pair14_i.scalar_type() == torch::kInt64, "pair14_i must be int64");
  TORCH_CHECK(pair14_j.scalar_type() == torch::kInt64, "pair14_j must be int64");
  TORCH_CHECK(
      pair_i.device() == coords.device() && pair_j.device() == coords.device() &&
          pair14_i.device() == coords.device() && pair14_j.device() == coords.device(),
      "pair index tensors must be on same device as coords");
  TORCH_CHECK(
      pair_coul_coeff.device() == coords.device() && pair_a12_coeff.device() == coords.device() &&
          pair_b6_coeff.device() == coords.device() && pair_b10_coeff.device() == coords.device() &&
          pair14_coul_coeff.device() == coords.device() && pair14_a12_coeff.device() == coords.device() &&
          pair14_b6_coeff.device() == coords.device() && pair14_b10_coeff.device() == coords.device(),
      "precomputed pair coeff tensors must be on same device as coords");
  TORCH_CHECK(
      pair_coul_coeff.scalar_type() == coords.scalar_type() &&
          pair_a12_coeff.scalar_type() == coords.scalar_type() &&
          pair_b6_coeff.scalar_type() == coords.scalar_type() &&
          pair_b10_coeff.scalar_type() == coords.scalar_type() &&
          pair14_coul_coeff.scalar_type() == coords.scalar_type() &&
          pair14_a12_coeff.scalar_type() == coords.scalar_type() &&
          pair14_b6_coeff.scalar_type() == coords.scalar_type() &&
          pair14_b10_coeff.scalar_type() == coords.scalar_type(),
      "precomputed pair coeff dtype mismatch with coords");
  TORCH_CHECK(pair_i.numel() == pair_j.numel(), "pair_i/pair_j size mismatch");
  TORCH_CHECK(pair14_i.numel() == pair14_j.numel(), "pair14_i/pair14_j size mismatch");
  TORCH_CHECK(
      pair_i.numel() == pair_coul_coeff.numel() && pair_i.numel() == pair_a12_coeff.numel() &&
          pair_i.numel() == pair_b6_coeff.numel() && pair_i.numel() == pair_b10_coeff.numel(),
      "general pair/coeff size mismatch");
  TORCH_CHECK(
      pair14_i.numel() == pair14_coul_coeff.numel() && pair14_i.numel() == pair14_a12_coeff.numel() &&
          pair14_i.numel() == pair14_b6_coeff.numel() && pair14_i.numel() == pair14_b10_coeff.numel(),
      "1-4 pair/coeff size mismatch");

  auto e_coul = torch::zeros({}, coords.options());
  auto e_lj = torch::zeros({}, coords.options());
  auto e_coul14 = torch::zeros({}, coords.options());
  auto e_lj14 = torch::zeros({}, coords.options());
  auto force = torch::zeros_like(coords);

  const int64_t n = pair_i.numel();
  for (int64_t start = 0; start < n; start += chunk_size) {
    int64_t end = std::min(start + chunk_size, n);
    auto ii = pair_i.slice(0, start, end);
    auto jj = pair_j.slice(0, start, end);
    auto cc = pair_coul_coeff.slice(0, start, end);
    auto a12 = pair_a12_coeff.slice(0, start, end);
    auto b6 = pair_b6_coeff.slice(0, start, end);
    auto b10 = pair_b10_coeff.slice(0, start, end);
    auto out = pair_energy_force_preparam_aten(coords, ii, jj, cc, a12, b6, b10);
    auto ce = std::get<0>(out);
    auto le = std::get<1>(out);
    auto fij = std::get<2>(out);
    e_coul = e_coul + ce;
    e_lj = e_lj + le;
    force.index_add_(0, ii, fij);
    force.index_add_(0, jj, -fij);
  }

  const int64_t n14 = pair14_i.numel();
  for (int64_t start = 0; start < n14; start += chunk_size) {
    int64_t end = std::min(start + chunk_size, n14);
    auto ii = pair14_i.slice(0, start, end);
    auto jj = pair14_j.slice(0, start, end);
    auto cc = pair14_coul_coeff.slice(0, start, end);
    auto a12 = pair14_a12_coeff.slice(0, start, end);
    auto b6 = pair14_b6_coeff.slice(0, start, end);
    auto b10 = pair14_b10_coeff.slice(0, start, end);
    auto out = pair_energy_force_preparam_aten(coords, ii, jj, cc, a12, b6, b10);
    auto ce = std::get<0>(out);
    auto le = std::get<1>(out);
    auto fij = std::get<2>(out);
    e_coul14 = e_coul14 + ce;
    e_lj14 = e_lj14 + le;
    force.index_add_(0, ii, fij);
    force.index_add_(0, jj, -fij);
  }

  return {e_coul, e_lj, e_coul14, e_lj14, force};
}

std::vector<torch::Tensor> nonbonded_energy_force_preparam_cpu(
    const torch::Tensor& coords,
    const torch::Tensor& pair_i,
    const torch::Tensor& pair_j,
    const torch::Tensor& pair_coul_coeff,
    const torch::Tensor& pair_a12_coeff,
    const torch::Tensor& pair_b6_coeff,
    const torch::Tensor& pair_b10_coeff,
    const torch::Tensor& pair14_i,
    const torch::Tensor& pair14_j,
    const torch::Tensor& pair14_coul_coeff,
    const torch::Tensor& pair14_a12_coeff,
    const torch::Tensor& pair14_b6_coeff,
    const torch::Tensor& pair14_b10_coeff,
    int64_t min_pairs_per_thread);

std::vector<torch::Tensor> nonbonded_energy_force_preparam(
    const torch::Tensor& coords,
    const torch::Tensor& pair_i,
    const torch::Tensor& pair_j,
    const torch::Tensor& pair_coul_coeff,
    const torch::Tensor& pair_a12_coeff,
    const torch::Tensor& pair_b6_coeff,
    const torch::Tensor& pair_b10_coeff,
    const torch::Tensor& pair14_i,
    const torch::Tensor& pair14_j,
    const torch::Tensor& pair14_coul_coeff,
    const torch::Tensor& pair14_a12_coeff,
    const torch::Tensor& pair14_b6_coeff,
    const torch::Tensor& pair14_b10_coeff,
    int64_t chunk_size,
    bool cpu_fast,
    int64_t min_pairs_per_thread) {
  if (coords.device().is_cpu() && cpu_fast) {
    return nonbonded_energy_force_preparam_cpu(
        coords,
        pair_i,
        pair_j,
        pair_coul_coeff,
        pair_a12_coeff,
        pair_b6_coeff,
        pair_b10_coeff,
        pair14_i,
        pair14_j,
        pair14_coul_coeff,
        pair14_a12_coeff,
        pair14_b6_coeff,
        pair14_b10_coeff,
        min_pairs_per_thread);
  }
  return nonbonded_energy_force_preparam_aten(
      coords,
      pair_i,
      pair_j,
      pair_coul_coeff,
      pair_a12_coeff,
      pair_b6_coeff,
      pair_b10_coeff,
      pair14_i,
      pair14_j,
      pair14_coul_coeff,
      pair14_a12_coeff,
      pair14_b6_coeff,
      pair14_b10_coeff,
      chunk_size);
}

template <typename scalar_t>
inline void accumulate_pair_one(
    const scalar_t* coords,
    const scalar_t* charge,
    const int64_t* atom_type,
    const scalar_t* lj_acoef,
    const scalar_t* lj_bcoef,
    int64_t n_lj,
    const scalar_t* hb_acoef,
    const scalar_t* hb_bcoef,
    int64_t n_hb,
    const int64_t* nb_index,
    int64_t ntypes,
    int64_t i,
    int64_t j,
    const scalar_t* inv_scee,
    const scalar_t* inv_scnb,
    int64_t p,
    scalar_t* force,
    double& e_coul,
    double& e_lj) {
  const double dx = static_cast<double>(coords[3 * j + 0] - coords[3 * i + 0]);
  const double dy = static_cast<double>(coords[3 * j + 1] - coords[3 * i + 1]);
  const double dz = static_cast<double>(coords[3 * j + 2] - coords[3 * i + 2]);
  const double r2 = std::max(dx * dx + dy * dy + dz * dz, 1.0e-24);
  const double inv_r = 1.0 / std::sqrt(r2);
  const double inv_r2 = inv_r * inv_r;

  const double qq = static_cast<double>(charge[i]) * static_cast<double>(charge[j]);
  double ec = qq * kCoulomb * inv_r;
  double f_coul = -(qq * kCoulomb) * inv_r2 * inv_r;

  const int64_t ti = atom_type[i];
  const int64_t tj = atom_type[j];
  const int64_t raw = nb_index[ti * ntypes + tj];

  const double inv_r4 = inv_r2 * inv_r2;
  const double inv_r6 = inv_r4 * inv_r2;
  const double inv_r8 = inv_r6 * inv_r2;
  const double inv_r10 = inv_r8 * inv_r2;
  const double inv_r12 = inv_r6 * inv_r6;
  const double inv_r14 = inv_r12 * inv_r2;

  double el = 0.0;
  double f_lj = 0.0;
  if (raw > 0) {
    const int64_t idx = raw - 1;
    if (idx >= 0 && idx < n_lj) {
      const double a = static_cast<double>(lj_acoef[idx]);
      const double b = static_cast<double>(lj_bcoef[idx]);
      el = a * inv_r12 - b * inv_r6;
      f_lj = -12.0 * a * inv_r14 + 6.0 * b * inv_r8;
    }
  } else if (raw < 0 && n_hb > 0) {
    const int64_t idx = -raw - 1;
    if (idx >= 0 && idx < n_hb) {
      const double a = static_cast<double>(hb_acoef[idx]);
      const double b = static_cast<double>(hb_bcoef[idx]);
      el = a * inv_r12 - b * inv_r10;
      f_lj = -12.0 * a * inv_r14 + 10.0 * b * inv_r12;
    }
  }

  if (inv_scee != nullptr) {
    const double s = static_cast<double>(inv_scee[p]);
    ec *= s;
    f_coul *= s;
  }
  if (inv_scnb != nullptr) {
    const double s = static_cast<double>(inv_scnb[p]);
    el *= s;
    f_lj *= s;
  }

  const double fs = f_coul + f_lj;
  const double fx = fs * dx;
  const double fy = fs * dy;
  const double fz = fs * dz;

  force[3 * i + 0] += static_cast<scalar_t>(fx);
  force[3 * i + 1] += static_cast<scalar_t>(fy);
  force[3 * i + 2] += static_cast<scalar_t>(fz);
  force[3 * j + 0] -= static_cast<scalar_t>(fx);
  force[3 * j + 1] -= static_cast<scalar_t>(fy);
  force[3 * j + 2] -= static_cast<scalar_t>(fz);

  e_coul += ec;
  e_lj += el;
}

template <typename scalar_t>
void accumulate_pairs_cpu(
    const scalar_t* coords,
    const scalar_t* charge,
    const int64_t* atom_type,
    const scalar_t* lj_acoef,
    const scalar_t* lj_bcoef,
    int64_t n_lj,
    const scalar_t* hb_acoef,
    const scalar_t* hb_bcoef,
    int64_t n_hb,
    const int64_t* nb_index,
    int64_t ntypes,
    const int64_t* ii,
    const int64_t* jj,
    const scalar_t* inv_scee,
    const scalar_t* inv_scnb,
    int64_t npairs,
    int64_t natom,
    int64_t min_pairs_per_thread,
    std::vector<scalar_t>& force_out,
    double& e_coul_out,
    double& e_lj_out) {
  const int64_t stride = natom * 3;
  if (npairs <= 0) {
    return;
  }

  int nthreads = 1;
#ifdef _OPENMP
  const int max_threads = effective_max_threads();
  // Adaptive threading:
  // For small pair counts, OpenMP launch/reduction overhead dominates.
  // Use fewer threads (or 1) to keep latency low on toy/small systems.
  // This threshold rule is an hessian_ff heuristic (not OpenMM-identical).
  // Priority: API argument -> env var -> internal default.
  const int64_t kMinPairsPerThread = min_pairs_per_thread_threshold(min_pairs_per_thread);
  nthreads = std::min<int>(
      max_threads,
      std::max<int64_t>(1, npairs / kMinPairsPerThread));
#endif

  if (nthreads <= 1) {
    scalar_t* f = force_out.data();
    for (int64_t p = 0; p < npairs; ++p) {
      accumulate_pair_one<scalar_t>(
          coords,
          charge,
          atom_type,
          lj_acoef,
          lj_bcoef,
          n_lj,
          hb_acoef,
          hb_bcoef,
          n_hb,
          nb_index,
          ntypes,
          ii[p],
          jj[p],
          inv_scee,
          inv_scnb,
          p,
          f,
          e_coul_out,
          e_lj_out);
    }
    return;
  }

  std::vector<scalar_t> force_tls(static_cast<size_t>(nthreads) * static_cast<size_t>(stride), static_cast<scalar_t>(0));
  std::vector<double> ec_tls(static_cast<size_t>(nthreads), 0.0);
  std::vector<double> el_tls(static_cast<size_t>(nthreads), 0.0);

#ifdef _OPENMP
#pragma omp parallel num_threads(nthreads)
#endif
  {
    int tid = 0;
#ifdef _OPENMP
    tid = omp_get_thread_num();
#endif
    auto* f = force_tls.data() + static_cast<size_t>(tid) * static_cast<size_t>(stride);
    double ec_loc = 0.0;
    double el_loc = 0.0;
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (int64_t p = 0; p < npairs; ++p) {
      accumulate_pair_one<scalar_t>(
          coords,
          charge,
          atom_type,
          lj_acoef,
          lj_bcoef,
          n_lj,
          hb_acoef,
          hb_bcoef,
          n_hb,
          nb_index,
          ntypes,
          ii[p],
          jj[p],
          inv_scee,
          inv_scnb,
          p,
          f,
          ec_loc,
          el_loc);
    }
    ec_tls[static_cast<size_t>(tid)] = ec_loc;
    el_tls[static_cast<size_t>(tid)] = el_loc;
  }

  for (int tid = 0; tid < nthreads; ++tid) {
    const auto* f = force_tls.data() + static_cast<size_t>(tid) * static_cast<size_t>(stride);
    for (int64_t k = 0; k < stride; ++k) {
      force_out[static_cast<size_t>(k)] += f[static_cast<size_t>(k)];
    }
    e_coul_out += ec_tls[static_cast<size_t>(tid)];
    e_lj_out += el_tls[static_cast<size_t>(tid)];
  }
}

template <typename scalar_t>
inline void accumulate_pair_preparam_one(
    const scalar_t* coords,
    int64_t i,
    int64_t j,
    const scalar_t* coul_coeff,
    const scalar_t* a12_coeff,
    const scalar_t* b6_coeff,
    const scalar_t* b10_coeff,
    int64_t p,
    scalar_t* force,
    double& e_coul,
    double& e_lj) {
  const double dx = static_cast<double>(coords[3 * j + 0] - coords[3 * i + 0]);
  const double dy = static_cast<double>(coords[3 * j + 1] - coords[3 * i + 1]);
  const double dz = static_cast<double>(coords[3 * j + 2] - coords[3 * i + 2]);
  const double r2 = std::max(dx * dx + dy * dy + dz * dz, 1.0e-24);
  const double inv_r = 1.0 / std::sqrt(r2);
  const double inv_r2 = inv_r * inv_r;
  const double inv_r4 = inv_r2 * inv_r2;
  const double inv_r6 = inv_r4 * inv_r2;
  const double inv_r8 = inv_r6 * inv_r2;
  const double inv_r10 = inv_r8 * inv_r2;
  const double inv_r12 = inv_r6 * inv_r6;
  const double inv_r14 = inv_r12 * inv_r2;

  const double c = static_cast<double>(coul_coeff[p]);
  const double a12 = static_cast<double>(a12_coeff[p]);
  const double b6 = static_cast<double>(b6_coeff[p]);
  const double b10 = static_cast<double>(b10_coeff[p]);

  const double ec = c * inv_r;
  const double el = a12 * inv_r12 - b6 * inv_r6 - b10 * inv_r10;
  const double fs = -c * inv_r2 * inv_r - 12.0 * a12 * inv_r14 + 6.0 * b6 * inv_r8 + 10.0 * b10 * inv_r12;

  const double fx = fs * dx;
  const double fy = fs * dy;
  const double fz = fs * dz;

  force[3 * i + 0] += static_cast<scalar_t>(fx);
  force[3 * i + 1] += static_cast<scalar_t>(fy);
  force[3 * i + 2] += static_cast<scalar_t>(fz);
  force[3 * j + 0] -= static_cast<scalar_t>(fx);
  force[3 * j + 1] -= static_cast<scalar_t>(fy);
  force[3 * j + 2] -= static_cast<scalar_t>(fz);

  e_coul += ec;
  e_lj += el;
}

template <typename scalar_t>
void accumulate_pairs_preparam_cpu(
    const scalar_t* coords,
    const int64_t* ii,
    const int64_t* jj,
    const scalar_t* coul_coeff,
    const scalar_t* a12_coeff,
    const scalar_t* b6_coeff,
    const scalar_t* b10_coeff,
    int64_t npairs,
    int64_t natom,
    int64_t min_pairs_per_thread,
    std::vector<scalar_t>& force_out,
    double& e_coul_out,
    double& e_lj_out) {
  const int64_t stride = natom * 3;
  if (npairs <= 0) {
    return;
  }

  int nthreads = 1;
#ifdef _OPENMP
  const int max_threads = effective_max_threads();
  const int64_t kMinPairsPerThread = min_pairs_per_thread_threshold(min_pairs_per_thread);
  nthreads = std::min<int>(
      max_threads,
      std::max<int64_t>(1, npairs / kMinPairsPerThread));
#endif

  if (nthreads <= 1) {
    scalar_t* f = force_out.data();
    for (int64_t p = 0; p < npairs; ++p) {
      accumulate_pair_preparam_one<scalar_t>(
          coords,
          ii[p],
          jj[p],
          coul_coeff,
          a12_coeff,
          b6_coeff,
          b10_coeff,
          p,
          f,
          e_coul_out,
          e_lj_out);
    }
    return;
  }

  std::vector<scalar_t> force_tls(
      static_cast<size_t>(nthreads) * static_cast<size_t>(stride),
      static_cast<scalar_t>(0));
  std::vector<double> ec_tls(static_cast<size_t>(nthreads), 0.0);
  std::vector<double> el_tls(static_cast<size_t>(nthreads), 0.0);

#ifdef _OPENMP
#pragma omp parallel num_threads(nthreads)
#endif
  {
    int tid = 0;
#ifdef _OPENMP
    tid = omp_get_thread_num();
#endif
    auto* f = force_tls.data() + static_cast<size_t>(tid) * static_cast<size_t>(stride);
    double ec_loc = 0.0;
    double el_loc = 0.0;
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (int64_t p = 0; p < npairs; ++p) {
      accumulate_pair_preparam_one<scalar_t>(
          coords,
          ii[p],
          jj[p],
          coul_coeff,
          a12_coeff,
          b6_coeff,
          b10_coeff,
          p,
          f,
          ec_loc,
          el_loc);
    }
    ec_tls[static_cast<size_t>(tid)] = ec_loc;
    el_tls[static_cast<size_t>(tid)] = el_loc;
  }

  for (int tid = 0; tid < nthreads; ++tid) {
    const auto* f = force_tls.data() + static_cast<size_t>(tid) * static_cast<size_t>(stride);
    for (int64_t k = 0; k < stride; ++k) {
      force_out[static_cast<size_t>(k)] += f[static_cast<size_t>(k)];
    }
    e_coul_out += ec_tls[static_cast<size_t>(tid)];
    e_lj_out += el_tls[static_cast<size_t>(tid)];
  }
}

std::vector<torch::Tensor> nonbonded_energy_force_cpu(
    const torch::Tensor& coords,
    const torch::Tensor& charge,
    const torch::Tensor& atom_type,
    const torch::Tensor& lj_acoef,
    const torch::Tensor& lj_bcoef,
    const torch::Tensor& hb_acoef,
    const torch::Tensor& hb_bcoef,
    const torch::Tensor& nb_index,
    const torch::Tensor& pair_i,
    const torch::Tensor& pair_j,
    const torch::Tensor& pair14_i,
    const torch::Tensor& pair14_j,
    const torch::Tensor& pair14_inv_scee,
    const torch::Tensor& pair14_inv_scnb,
    int64_t min_pairs_per_thread) {
  auto c = coords.contiguous();
  auto q = charge.contiguous();
  auto at = atom_type.contiguous();
  auto acoef = lj_acoef.contiguous();
  auto bcoef = lj_bcoef.contiguous();
  auto ha = hb_acoef.contiguous();
  auto hb = hb_bcoef.contiguous();
  auto nbi = nb_index.contiguous();
  auto pi = pair_i.contiguous();
  auto pj = pair_j.contiguous();
  auto p14i = pair14_i.contiguous();
  auto p14j = pair14_j.contiguous();
  auto scee = pair14_inv_scee.contiguous();
  auto scnb = pair14_inv_scnb.contiguous();

  TORCH_CHECK(c.device().is_cpu(), "CPU kernel expects CPU coords tensor");
  TORCH_CHECK(c.dim() == 2 && c.size(1) == 3, "coords must be [N,3]");
  TORCH_CHECK(c.scalar_type() == q.scalar_type(), "coords/charge dtype mismatch");
  TORCH_CHECK(c.scalar_type() == acoef.scalar_type(), "coords/lj dtype mismatch");
  TORCH_CHECK(c.scalar_type() == bcoef.scalar_type(), "coords/lj dtype mismatch");
  TORCH_CHECK(nbi.scalar_type() == torch::kInt64, "nb_index must be int64");
  TORCH_CHECK(pi.scalar_type() == torch::kInt64, "pair_i must be int64");
  TORCH_CHECK(pj.scalar_type() == torch::kInt64, "pair_j must be int64");
  TORCH_CHECK(p14i.scalar_type() == torch::kInt64, "pair14_i must be int64");
  TORCH_CHECK(p14j.scalar_type() == torch::kInt64, "pair14_j must be int64");

  const int64_t natom = c.size(0);
  const int64_t ntypes = nbi.size(0);
  auto e_coul = torch::zeros({}, c.options());
  auto e_lj = torch::zeros({}, c.options());
  auto e_coul14 = torch::zeros({}, c.options());
  auto e_lj14 = torch::zeros({}, c.options());
  auto force = torch::zeros_like(c);

  AT_DISPATCH_FLOATING_TYPES(c.scalar_type(), "nonbonded_energy_force_cpu", [&] {
    std::vector<scalar_t> force_acc(static_cast<size_t>(natom) * 3, static_cast<scalar_t>(0));
    double ec = 0.0;
    double el = 0.0;
    double ec14 = 0.0;
    double el14 = 0.0;

    const scalar_t* hb_a_ptr = ha.numel() > 0 ? ha.data_ptr<scalar_t>() : nullptr;
    const scalar_t* hb_b_ptr = hb.numel() > 0 ? hb.data_ptr<scalar_t>() : nullptr;

    accumulate_pairs_cpu<scalar_t>(
        c.data_ptr<scalar_t>(),
        q.data_ptr<scalar_t>(),
        at.data_ptr<int64_t>(),
        acoef.data_ptr<scalar_t>(),
        bcoef.data_ptr<scalar_t>(),
        acoef.numel(),
        hb_a_ptr,
        hb_b_ptr,
        ha.numel(),
        nbi.data_ptr<int64_t>(),
        ntypes,
        pi.data_ptr<int64_t>(),
        pj.data_ptr<int64_t>(),
        nullptr,
        nullptr,
        pi.numel(),
        natom,
        min_pairs_per_thread,
        force_acc,
        ec,
        el);

    const scalar_t* scee_ptr = scee.numel() > 0 ? scee.data_ptr<scalar_t>() : nullptr;
    const scalar_t* scnb_ptr = scnb.numel() > 0 ? scnb.data_ptr<scalar_t>() : nullptr;
    accumulate_pairs_cpu<scalar_t>(
        c.data_ptr<scalar_t>(),
        q.data_ptr<scalar_t>(),
        at.data_ptr<int64_t>(),
        acoef.data_ptr<scalar_t>(),
        bcoef.data_ptr<scalar_t>(),
        acoef.numel(),
        hb_a_ptr,
        hb_b_ptr,
        ha.numel(),
        nbi.data_ptr<int64_t>(),
        ntypes,
        p14i.data_ptr<int64_t>(),
        p14j.data_ptr<int64_t>(),
        scee_ptr,
        scnb_ptr,
        p14i.numel(),
        natom,
        min_pairs_per_thread,
        force_acc,
        ec14,
        el14);

    auto* force_ptr = force.data_ptr<scalar_t>();
    std::copy(force_acc.begin(), force_acc.end(), force_ptr);
    e_coul.fill_(static_cast<scalar_t>(ec));
    e_lj.fill_(static_cast<scalar_t>(el));
    e_coul14.fill_(static_cast<scalar_t>(ec14));
    e_lj14.fill_(static_cast<scalar_t>(el14));
  });

  return {e_coul, e_lj, e_coul14, e_lj14, force};
}

std::vector<torch::Tensor> nonbonded_energy_force_preparam_cpu(
    const torch::Tensor& coords,
    const torch::Tensor& pair_i,
    const torch::Tensor& pair_j,
    const torch::Tensor& pair_coul_coeff,
    const torch::Tensor& pair_a12_coeff,
    const torch::Tensor& pair_b6_coeff,
    const torch::Tensor& pair_b10_coeff,
    const torch::Tensor& pair14_i,
    const torch::Tensor& pair14_j,
    const torch::Tensor& pair14_coul_coeff,
    const torch::Tensor& pair14_a12_coeff,
    const torch::Tensor& pair14_b6_coeff,
    const torch::Tensor& pair14_b10_coeff,
    int64_t min_pairs_per_thread) {
  auto c = coords.contiguous();
  auto pi = pair_i.contiguous();
  auto pj = pair_j.contiguous();
  auto cc = pair_coul_coeff.contiguous();
  auto a12 = pair_a12_coeff.contiguous();
  auto b6 = pair_b6_coeff.contiguous();
  auto b10 = pair_b10_coeff.contiguous();
  auto p14i = pair14_i.contiguous();
  auto p14j = pair14_j.contiguous();
  auto cc14 = pair14_coul_coeff.contiguous();
  auto a1214 = pair14_a12_coeff.contiguous();
  auto b614 = pair14_b6_coeff.contiguous();
  auto b1014 = pair14_b10_coeff.contiguous();

  TORCH_CHECK(c.device().is_cpu(), "CPU kernel expects CPU coords tensor");
  TORCH_CHECK(c.dim() == 2 && c.size(1) == 3, "coords must be [N,3]");
  TORCH_CHECK(pi.scalar_type() == torch::kInt64, "pair_i must be int64");
  TORCH_CHECK(pj.scalar_type() == torch::kInt64, "pair_j must be int64");
  TORCH_CHECK(p14i.scalar_type() == torch::kInt64, "pair14_i must be int64");
  TORCH_CHECK(p14j.scalar_type() == torch::kInt64, "pair14_j must be int64");
  TORCH_CHECK(cc.scalar_type() == c.scalar_type(), "pair coeff dtype mismatch");
  TORCH_CHECK(a12.scalar_type() == c.scalar_type(), "pair coeff dtype mismatch");
  TORCH_CHECK(b6.scalar_type() == c.scalar_type(), "pair coeff dtype mismatch");
  TORCH_CHECK(b10.scalar_type() == c.scalar_type(), "pair coeff dtype mismatch");
  TORCH_CHECK(cc14.scalar_type() == c.scalar_type(), "pair14 coeff dtype mismatch");
  TORCH_CHECK(a1214.scalar_type() == c.scalar_type(), "pair14 coeff dtype mismatch");
  TORCH_CHECK(b614.scalar_type() == c.scalar_type(), "pair14 coeff dtype mismatch");
  TORCH_CHECK(b1014.scalar_type() == c.scalar_type(), "pair14 coeff dtype mismatch");
  TORCH_CHECK(pi.numel() == pj.numel(), "pair_i/pair_j size mismatch");
  TORCH_CHECK(p14i.numel() == p14j.numel(), "pair14_i/pair14_j size mismatch");
  TORCH_CHECK(
      pi.numel() == cc.numel() && pi.numel() == a12.numel() && pi.numel() == b6.numel() &&
          pi.numel() == b10.numel(),
      "pair/coeff size mismatch");
  TORCH_CHECK(
      p14i.numel() == cc14.numel() && p14i.numel() == a1214.numel() && p14i.numel() == b614.numel() &&
          p14i.numel() == b1014.numel(),
      "pair14/coeff size mismatch");

  const int64_t natom = c.size(0);
  auto e_coul = torch::zeros({}, c.options());
  auto e_lj = torch::zeros({}, c.options());
  auto e_coul14 = torch::zeros({}, c.options());
  auto e_lj14 = torch::zeros({}, c.options());
  auto force = torch::zeros_like(c);

  AT_DISPATCH_FLOATING_TYPES(c.scalar_type(), "nonbonded_energy_force_preparam_cpu", [&] {
    std::vector<scalar_t> force_acc(static_cast<size_t>(natom) * 3, static_cast<scalar_t>(0));
    double ec = 0.0;
    double el = 0.0;
    double ec14 = 0.0;
    double el14 = 0.0;

    accumulate_pairs_preparam_cpu<scalar_t>(
        c.data_ptr<scalar_t>(),
        pi.data_ptr<int64_t>(),
        pj.data_ptr<int64_t>(),
        cc.data_ptr<scalar_t>(),
        a12.data_ptr<scalar_t>(),
        b6.data_ptr<scalar_t>(),
        b10.data_ptr<scalar_t>(),
        pi.numel(),
        natom,
        min_pairs_per_thread,
        force_acc,
        ec,
        el);

    accumulate_pairs_preparam_cpu<scalar_t>(
        c.data_ptr<scalar_t>(),
        p14i.data_ptr<int64_t>(),
        p14j.data_ptr<int64_t>(),
        cc14.data_ptr<scalar_t>(),
        a1214.data_ptr<scalar_t>(),
        b614.data_ptr<scalar_t>(),
        b1014.data_ptr<scalar_t>(),
        p14i.numel(),
        natom,
        min_pairs_per_thread,
        force_acc,
        ec14,
        el14);

    auto* force_ptr = force.data_ptr<scalar_t>();
    std::copy(force_acc.begin(), force_acc.end(), force_ptr);
    e_coul.fill_(static_cast<scalar_t>(ec));
    e_lj.fill_(static_cast<scalar_t>(el));
    e_coul14.fill_(static_cast<scalar_t>(ec14));
    e_lj14.fill_(static_cast<scalar_t>(el14));
  });

  return {e_coul, e_lj, e_coul14, e_lj14, force};
}

}  // namespace

std::vector<torch::Tensor> nonbonded_energy_force(
    const torch::Tensor& coords,
    const torch::Tensor& charge,
    const torch::Tensor& atom_type,
    const torch::Tensor& lj_acoef,
    const torch::Tensor& lj_bcoef,
    const torch::Tensor& hb_acoef,
    const torch::Tensor& hb_bcoef,
    const torch::Tensor& nb_index,
    const torch::Tensor& pair_i,
    const torch::Tensor& pair_j,
    const torch::Tensor& pair14_i,
    const torch::Tensor& pair14_j,
    const torch::Tensor& pair14_inv_scee,
    const torch::Tensor& pair14_inv_scnb,
    int64_t chunk_size,
    bool cpu_fast,
    int64_t min_pairs_per_thread) {
  // CPU-specialized path: hand-written loop kernel with OpenMP threading.
  if (coords.device().is_cpu() && cpu_fast) {
    return nonbonded_energy_force_cpu(
        coords,
        charge,
        atom_type,
        lj_acoef,
        lj_bcoef,
        hb_acoef,
        hb_bcoef,
        nb_index,
        pair_i,
        pair_j,
        pair14_i,
        pair14_j,
        pair14_inv_scee,
        pair14_inv_scnb,
        min_pairs_per_thread);
  }
  // ATen path for autograd-friendly CPU execution.
  return nonbonded_energy_force_aten(
      coords,
      charge,
      atom_type,
      lj_acoef,
      lj_bcoef,
      hb_acoef,
      hb_bcoef,
      nb_index,
      pair_i,
      pair_j,
      pair14_i,
      pair14_j,
      pair14_inv_scee,
      pair14_inv_scnb,
      chunk_size);
}

PYBIND11_MODULE(TORCH_EXTENSION_NAME, m) {
  m.doc() = "hessian_ff native nonbonded extension";
  m.def(
      "nonbonded_energy_force",
      &nonbonded_energy_force,
      "Compute nonbonded energies and forces (CPU-only backend)",
      pybind11::arg("coords"),
      pybind11::arg("charge"),
      pybind11::arg("atom_type"),
      pybind11::arg("lj_acoef"),
      pybind11::arg("lj_bcoef"),
      pybind11::arg("hb_acoef"),
      pybind11::arg("hb_bcoef"),
      pybind11::arg("nb_index"),
      pybind11::arg("pair_i"),
      pybind11::arg("pair_j"),
      pybind11::arg("pair14_i"),
      pybind11::arg("pair14_j"),
      pybind11::arg("pair14_inv_scee"),
      pybind11::arg("pair14_inv_scnb"),
      pybind11::arg("chunk_size"),
      pybind11::arg("cpu_fast") = true,
      pybind11::arg("min_pairs_per_thread") = -1);
  m.def(
      "nonbonded_energy_force_preparam",
      &nonbonded_energy_force_preparam,
      "Compute nonbonded energies and forces from precomputed pair coefficients (CPU-only)",
      pybind11::arg("coords"),
      pybind11::arg("pair_i"),
      pybind11::arg("pair_j"),
      pybind11::arg("pair_coul_coeff"),
      pybind11::arg("pair_a12_coeff"),
      pybind11::arg("pair_b6_coeff"),
      pybind11::arg("pair_b10_coeff"),
      pybind11::arg("pair14_i"),
      pybind11::arg("pair14_j"),
      pybind11::arg("pair14_coul_coeff"),
      pybind11::arg("pair14_a12_coeff"),
      pybind11::arg("pair14_b6_coeff"),
      pybind11::arg("pair14_b10_coeff"),
      pybind11::arg("chunk_size"),
      pybind11::arg("cpu_fast") = true,
      pybind11::arg("min_pairs_per_thread") = -1);
}
