#include <torch/extension.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

// References:
// - hessian_ff term definitions and analytical force equations:
//   - hessian_ff/terms/bond.py
//   - hessian_ff/terms/angle.py
//   - hessian_ff/terms/dihedral.py
//   - hessian_ff/terms/cmap.py
// - OpenMM conventions used for compatibility:
//   - Periodic torsion / dihedral sign conventions:
//     https://github.com/openmm/openmm/blob/master/platforms/reference/src/SimTKReference/ReferenceProperDihedralBond.cpp
//   - CMAP bicubic map conventions:
//     https://github.com/openmm/openmm/blob/master/openmmapi/src/CMAPTorsionForceImpl.cpp

namespace {

constexpr double kTwoPi = 6.283185307179586476925286766559;

inline double clamp(double x, double lo, double hi) {
  return std::max(lo, std::min(hi, x));
}

inline double wrap_0_2pi(double x) {
  double y = std::fmod(x, kTwoPi);
  if (y < 0.0) {
    y += kTwoPi;
  }
  return y;
}

inline void cross3(
    double ax,
    double ay,
    double az,
    double bx,
    double by,
    double bz,
    double& cx,
    double& cy,
    double& cz) {
  cx = ay * bz - az * by;
  cy = az * bx - ax * bz;
  cz = ax * by - ay * bx;
}

inline double dot3(
    double ax,
    double ay,
    double az,
    double bx,
    double by,
    double bz) {
  return ax * bx + ay * by + az * bz;
}

template <typename scalar_t>
inline void add_force(scalar_t* force, int64_t atom, double fx, double fy, double fz) {
  force[3 * atom + 0] += static_cast<scalar_t>(fx);
  force[3 * atom + 1] += static_cast<scalar_t>(fy);
  force[3 * atom + 2] += static_cast<scalar_t>(fz);
}

template <typename scalar_t>
inline double dihedral_angle_from_points(
    const scalar_t* coords,
    int64_t i,
    int64_t j,
    int64_t k,
    int64_t l,
    double* out_cross12,
    double* out_cross23,
    double* out_v1,
    double* out_v2,
    double* out_v3) {
  // OpenMM-compatible signed dihedral convention used in hessian_ff/terms/dihedral.py.
  const double p0x = static_cast<double>(coords[3 * i + 0]);
  const double p0y = static_cast<double>(coords[3 * i + 1]);
  const double p0z = static_cast<double>(coords[3 * i + 2]);
  const double p1x = static_cast<double>(coords[3 * j + 0]);
  const double p1y = static_cast<double>(coords[3 * j + 1]);
  const double p1z = static_cast<double>(coords[3 * j + 2]);
  const double p2x = static_cast<double>(coords[3 * k + 0]);
  const double p2y = static_cast<double>(coords[3 * k + 1]);
  const double p2z = static_cast<double>(coords[3 * k + 2]);
  const double p3x = static_cast<double>(coords[3 * l + 0]);
  const double p3y = static_cast<double>(coords[3 * l + 1]);
  const double p3z = static_cast<double>(coords[3 * l + 2]);

  const double v1x = p0x - p1x;
  const double v1y = p0y - p1y;
  const double v1z = p0z - p1z;
  const double v2x = p2x - p1x;
  const double v2y = p2y - p1y;
  const double v2z = p2z - p1z;
  const double v3x = p2x - p3x;
  const double v3y = p2y - p3y;
  const double v3z = p2z - p3z;

  double c12x, c12y, c12z;
  double c23x, c23y, c23z;
  cross3(v1x, v1y, v1z, v2x, v2y, v2z, c12x, c12y, c12z);
  cross3(v2x, v2y, v2z, v3x, v3y, v3z, c23x, c23y, c23z);

  const double n12 = std::max(std::sqrt(dot3(c12x, c12y, c12z, c12x, c12y, c12z)), 1.0e-12);
  const double n23 = std::max(std::sqrt(dot3(c23x, c23y, c23z, c23x, c23y, c23z)), 1.0e-12);
  const double cos_phi = clamp(dot3(c12x, c12y, c12z, c23x, c23y, c23z) / (n12 * n23), -1.0, 1.0);
  double phi = std::acos(cos_phi);

  const double sign_probe = dot3(v1x, v1y, v1z, c23x, c23y, c23z);
  if (sign_probe < 0.0) {
    phi = -phi;
  }

  if (out_cross12 != nullptr) {
    out_cross12[0] = c12x;
    out_cross12[1] = c12y;
    out_cross12[2] = c12z;
  }
  if (out_cross23 != nullptr) {
    out_cross23[0] = c23x;
    out_cross23[1] = c23y;
    out_cross23[2] = c23z;
  }
  if (out_v1 != nullptr) {
    out_v1[0] = v1x;
    out_v1[1] = v1y;
    out_v1[2] = v1z;
  }
  if (out_v2 != nullptr) {
    out_v2[0] = v2x;
    out_v2[1] = v2y;
    out_v2[2] = v2z;
  }
  if (out_v3 != nullptr) {
    out_v3[0] = v3x;
    out_v3[1] = v3y;
    out_v3[2] = v3z;
  }
  return phi;
}

template <typename scalar_t>
inline void accumulate_dihedral_force(
    const scalar_t* coords,
    int64_t i,
    int64_t j,
    int64_t k,
    int64_t l,
    double dE_dphi,
    scalar_t* force) {
  double cross12[3], cross23[3], v1[3], v2[3], v3[3];
  (void)dihedral_angle_from_points(
      coords,
      i,
      j,
      k,
      l,
      cross12,
      cross23,
      v1,
      v2,
      v3);

  const double norm_cross12_sq = std::max(dot3(cross12[0], cross12[1], cross12[2], cross12[0], cross12[1], cross12[2]), 1.0e-24);
  const double norm_cross23_sq = std::max(dot3(cross23[0], cross23[1], cross23[2], cross23[0], cross23[1], cross23[2]), 1.0e-24);
  const double norm_v2_sq = std::max(dot3(v2[0], v2[1], v2[2], v2[0], v2[1], v2[2]), 1.0e-24);
  const double norm_v2 = std::max(std::sqrt(norm_v2_sq), 1.0e-12);

  const double f0 = (-dE_dphi * norm_v2) / norm_cross12_sq;
  const double f3 = (dE_dphi * norm_v2) / norm_cross23_sq;
  const double f1 = dot3(v1[0], v1[1], v1[2], v2[0], v2[1], v2[2]) / norm_v2_sq;
  const double f2 = dot3(v3[0], v3[1], v3[2], v2[0], v2[1], v2[2]) / norm_v2_sq;

  const double ff0x = f0 * cross12[0];
  const double ff0y = f0 * cross12[1];
  const double ff0z = f0 * cross12[2];
  const double ff3x = f3 * cross23[0];
  const double ff3y = f3 * cross23[1];
  const double ff3z = f3 * cross23[2];

  const double sx = f1 * ff0x - f2 * ff3x;
  const double sy = f1 * ff0y - f2 * ff3y;
  const double sz = f1 * ff0z - f2 * ff3z;

  const double ff1x = ff0x - sx;
  const double ff1y = ff0y - sy;
  const double ff1z = ff0z - sz;
  const double ff2x = ff3x + sx;
  const double ff2y = ff3y + sy;
  const double ff2z = ff3z + sz;

  add_force(force, i, ff0x, ff0y, ff0z);
  add_force(force, j, -ff1x, -ff1y, -ff1z);
  add_force(force, k, -ff2x, -ff2y, -ff2z);
  add_force(force, l, ff3x, ff3y, ff3z);
}

std::vector<torch::Tensor> bonded_energy_force_cpu(
    const torch::Tensor& coords,
    const torch::Tensor& bond_i,
    const torch::Tensor& bond_j,
    const torch::Tensor& bond_k,
    const torch::Tensor& bond_r0,
    const torch::Tensor& angle_i,
    const torch::Tensor& angle_j,
    const torch::Tensor& angle_k,
    const torch::Tensor& angle_k0,
    const torch::Tensor& angle_t0,
    const torch::Tensor& dihed_i,
    const torch::Tensor& dihed_j,
    const torch::Tensor& dihed_k,
    const torch::Tensor& dihed_l,
    const torch::Tensor& dihed_force,
    const torch::Tensor& dihed_period,
    const torch::Tensor& dihed_phase,
    const torch::Tensor& cmap_type,
    const torch::Tensor& cmap_i,
    const torch::Tensor& cmap_j,
    const torch::Tensor& cmap_k,
    const torch::Tensor& cmap_l,
    const torch::Tensor& cmap_m,
    const torch::Tensor& cmap_size,
    const torch::Tensor& cmap_delta,
    const torch::Tensor& cmap_offset,
    const torch::Tensor& cmap_coeff) {
  auto c = coords.contiguous();
  TORCH_CHECK(c.device().is_cpu(), "bonded_energy_force_cpu expects CPU coords");
  TORCH_CHECK(c.dim() == 2 && c.size(1) == 3, "coords must be [N,3]");

  auto bi = bond_i.contiguous();
  auto bj = bond_j.contiguous();
  auto bk = bond_k.contiguous();
  auto br0 = bond_r0.contiguous();

  auto ai = angle_i.contiguous();
  auto aj = angle_j.contiguous();
  auto ak = angle_k.contiguous();
  auto ak0 = angle_k0.contiguous();
  auto at0 = angle_t0.contiguous();

  auto di = dihed_i.contiguous();
  auto dj = dihed_j.contiguous();
  auto dk = dihed_k.contiguous();
  auto dl = dihed_l.contiguous();
  auto df = dihed_force.contiguous();
  auto dp = dihed_period.contiguous();
  auto dph = dihed_phase.contiguous();

  auto ct = cmap_type.contiguous();
  auto ci = cmap_i.contiguous();
  auto cj = cmap_j.contiguous();
  auto ck = cmap_k.contiguous();
  auto cl = cmap_l.contiguous();
  auto cm = cmap_m.contiguous();
  auto csize = cmap_size.contiguous();
  auto cdelta = cmap_delta.contiguous();
  auto coff = cmap_offset.contiguous();
  auto ccoef = cmap_coeff.contiguous();

  TORCH_CHECK(bi.scalar_type() == torch::kInt64, "bond_i must be int64");
  TORCH_CHECK(ai.scalar_type() == torch::kInt64, "angle_i must be int64");
  TORCH_CHECK(di.scalar_type() == torch::kInt64, "dihed_i must be int64");
  TORCH_CHECK(ct.scalar_type() == torch::kInt64, "cmap_type must be int64");
  TORCH_CHECK(csize.scalar_type() == torch::kInt64, "cmap_size must be int64");
  TORCH_CHECK(coff.scalar_type() == torch::kInt64, "cmap_offset must be int64");

  auto force = torch::zeros_like(c);
  auto e_bond = torch::zeros({}, c.options());
  auto e_angle = torch::zeros({}, c.options());
  auto e_dihed = torch::zeros({}, c.options());
  auto e_cmap = torch::zeros({}, c.options());

  AT_DISPATCH_FLOATING_TYPES(c.scalar_type(), "bonded_energy_force_cpu", [&] {
    const scalar_t* xyz = c.data_ptr<scalar_t>();
    scalar_t* f = force.data_ptr<scalar_t>();
    const int64_t natom = c.size(0);
    const int64_t stride = 3 * natom;
    const int64_t nbond = bi.numel();
    const int64_t nangle = ai.numel();
    const int64_t ndihed = di.numel();
    const int64_t ncmap = ct.numel();

    const int64_t* bi_ptr = bi.data_ptr<int64_t>();
    const int64_t* bj_ptr = bj.data_ptr<int64_t>();
    const scalar_t* bk_ptr = bk.data_ptr<scalar_t>();
    const scalar_t* br0_ptr = br0.data_ptr<scalar_t>();

    const int64_t* ai_ptr = ai.data_ptr<int64_t>();
    const int64_t* aj_ptr = aj.data_ptr<int64_t>();
    const int64_t* ak_ptr = ak.data_ptr<int64_t>();
    const scalar_t* ak0_ptr = ak0.data_ptr<scalar_t>();
    const scalar_t* at0_ptr = at0.data_ptr<scalar_t>();

    const int64_t* di_ptr = di.data_ptr<int64_t>();
    const int64_t* dj_ptr = dj.data_ptr<int64_t>();
    const int64_t* dk_ptr = dk.data_ptr<int64_t>();
    const int64_t* dl_ptr = dl.data_ptr<int64_t>();
    const scalar_t* df_ptr = df.data_ptr<scalar_t>();
    const scalar_t* dp_ptr = dp.data_ptr<scalar_t>();
    const scalar_t* dph_ptr = dph.data_ptr<scalar_t>();

    const int64_t* ct_ptr = ct.data_ptr<int64_t>();
    const int64_t* ci_ptr = ci.data_ptr<int64_t>();
    const int64_t* cj_ptr = cj.data_ptr<int64_t>();
    const int64_t* ck_ptr = ck.data_ptr<int64_t>();
    const int64_t* cl_ptr = cl.data_ptr<int64_t>();
    const int64_t* cm_ptr = cm.data_ptr<int64_t>();
    const int64_t* csize_ptr = csize.data_ptr<int64_t>();
    const scalar_t* cdelta_ptr = cdelta.data_ptr<scalar_t>();
    const int64_t* coff_ptr = coff.data_ptr<int64_t>();
    const scalar_t* ccoef_ptr = ccoef.data_ptr<scalar_t>();

    int nthreads = 1;
#ifdef _OPENMP
    nthreads = std::max(1, omp_get_max_threads());
#endif
    std::vector<scalar_t> force_tls(
        static_cast<size_t>(nthreads) * static_cast<size_t>(stride),
        static_cast<scalar_t>(0));
    std::vector<double> eb_tls(static_cast<size_t>(nthreads), 0.0);
    std::vector<double> ea_tls(static_cast<size_t>(nthreads), 0.0);
    std::vector<double> ed_tls(static_cast<size_t>(nthreads), 0.0);
    std::vector<double> ec_tls(static_cast<size_t>(nthreads), 0.0);

    // Bond
#ifdef _OPENMP
#pragma omp parallel num_threads(nthreads)
#endif
    {
      int tid = 0;
#ifdef _OPENMP
      tid = omp_get_thread_num();
#endif
      scalar_t* fl = force_tls.data() + static_cast<size_t>(tid) * static_cast<size_t>(stride);
      double eb_local = 0.0;
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
      for (int64_t p = 0; p < nbond; ++p) {
        const int64_t i = bi_ptr[p];
        const int64_t j = bj_ptr[p];
        const double k = static_cast<double>(bk_ptr[p]);
        const double r0 = static_cast<double>(br0_ptr[p]);

        const double dx = static_cast<double>(xyz[3 * j + 0] - xyz[3 * i + 0]);
        const double dy = static_cast<double>(xyz[3 * j + 1] - xyz[3 * i + 1]);
        const double dz = static_cast<double>(xyz[3 * j + 2] - xyz[3 * i + 2]);
        const double r2 = std::max(dx * dx + dy * dy + dz * dz, 1.0e-24);
        const double inv_r = 1.0 / std::sqrt(r2);
        const double r = r2 * inv_r;
        const double dr = r - r0;
        eb_local += k * dr * dr;
        const double fs = 2.0 * k * dr * inv_r;
        add_force(fl, i, fs * dx, fs * dy, fs * dz);
        add_force(fl, j, -fs * dx, -fs * dy, -fs * dz);
      }
      eb_tls[static_cast<size_t>(tid)] += eb_local;
    }

    // Angle
#ifdef _OPENMP
#pragma omp parallel num_threads(nthreads)
#endif
    {
      int tid = 0;
#ifdef _OPENMP
      tid = omp_get_thread_num();
#endif
      scalar_t* fl = force_tls.data() + static_cast<size_t>(tid) * static_cast<size_t>(stride);
      double ea_local = 0.0;
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
      for (int64_t p = 0; p < nangle; ++p) {
        const int64_t i = ai_ptr[p];
        const int64_t j = aj_ptr[p];
        const int64_t k = ak_ptr[p];
        const double ktheta = static_cast<double>(ak0_ptr[p]);
        const double theta0 = static_cast<double>(at0_ptr[p]);

        const double d0x = static_cast<double>(xyz[3 * j + 0] - xyz[3 * i + 0]);
        const double d0y = static_cast<double>(xyz[3 * j + 1] - xyz[3 * i + 1]);
        const double d0z = static_cast<double>(xyz[3 * j + 2] - xyz[3 * i + 2]);
        const double d1x = static_cast<double>(xyz[3 * j + 0] - xyz[3 * k + 0]);
        const double d1y = static_cast<double>(xyz[3 * j + 1] - xyz[3 * k + 1]);
        const double d1z = static_cast<double>(xyz[3 * j + 2] - xyz[3 * k + 2]);

        double px, py, pz;
        cross3(d0x, d0y, d0z, d1x, d1y, d1z, px, py, pz);

        const double r20 = std::max(dot3(d0x, d0y, d0z, d0x, d0y, d0z), 1.0e-24);
        const double r21 = std::max(dot3(d1x, d1y, d1z, d1x, d1y, d1z), 1.0e-24);
        const double rp = std::max(std::sqrt(dot3(px, py, pz, px, py, pz)), 1.0e-12);
        const double dot = dot3(d0x, d0y, d0z, d1x, d1y, d1z);
        const double cos_theta = clamp(dot / std::sqrt(r20 * r21), -1.0, 1.0);
        const double theta = std::acos(cos_theta);

        const double dtheta = theta - theta0;
        ea_local += ktheta * dtheta * dtheta;
        const double dE_dtheta = 2.0 * ktheta * dtheta;

        const double term_i = dE_dtheta / (r20 * rp);
        const double term_k = -dE_dtheta / (r21 * rp);

        double cix, ciy, ciz;
        double ckx, cky, ckz;
        cross3(d0x, d0y, d0z, px, py, pz, cix, ciy, ciz);
        cross3(d1x, d1y, d1z, px, py, pz, ckx, cky, ckz);

        const double fix = cix * term_i;
        const double fiy = ciy * term_i;
        const double fiz = ciz * term_i;
        const double fkx = ckx * term_k;
        const double fky = cky * term_k;
        const double fkz = ckz * term_k;
        const double fjx = -(fix + fkx);
        const double fjy = -(fiy + fky);
        const double fjz = -(fiz + fkz);

        add_force(fl, i, fix, fiy, fiz);
        add_force(fl, j, fjx, fjy, fjz);
        add_force(fl, k, fkx, fky, fkz);
      }
      ea_tls[static_cast<size_t>(tid)] += ea_local;
    }

    // Dihedral
#ifdef _OPENMP
#pragma omp parallel num_threads(nthreads)
#endif
    {
      int tid = 0;
#ifdef _OPENMP
      tid = omp_get_thread_num();
#endif
      scalar_t* fl = force_tls.data() + static_cast<size_t>(tid) * static_cast<size_t>(stride);
      double ed_local = 0.0;
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
      for (int64_t p = 0; p < ndihed; ++p) {
        const int64_t i = di_ptr[p];
        const int64_t j = dj_ptr[p];
        const int64_t k = dk_ptr[p];
        const int64_t l = dl_ptr[p];
        const double kf = static_cast<double>(df_ptr[p]);
        const double n = std::abs(static_cast<double>(dp_ptr[p]));
        const double phase = static_cast<double>(dph_ptr[p]);

        const double phi = dihedral_angle_from_points<scalar_t>(
            xyz, i, j, k, l, nullptr, nullptr, nullptr, nullptr, nullptr);
        const double delta = n * phi - phase;
        ed_local += kf * (1.0 + std::cos(delta));
        const double dE_dphi = -kf * n * std::sin(delta);
        accumulate_dihedral_force<scalar_t>(xyz, i, j, k, l, dE_dphi, fl);
      }
      ed_tls[static_cast<size_t>(tid)] += ed_local;
    }

    // CMAP
#ifdef _OPENMP
#pragma omp parallel num_threads(nthreads)
#endif
    {
      int tid = 0;
#ifdef _OPENMP
      tid = omp_get_thread_num();
#endif
      scalar_t* fl = force_tls.data() + static_cast<size_t>(tid) * static_cast<size_t>(stride);
      double ec_local = 0.0;
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
      for (int64_t p = 0; p < ncmap; ++p) {
        const int64_t tmap = ct_ptr[p];
        const int64_t i = ci_ptr[p];
        const int64_t j = cj_ptr[p];
        const int64_t k = ck_ptr[p];
        const int64_t l = cl_ptr[p];
        const int64_t m = cm_ptr[p];

        const double phi = dihedral_angle_from_points<scalar_t>(
            xyz, i, j, k, l, nullptr, nullptr, nullptr, nullptr, nullptr);
        const double psi = dihedral_angle_from_points<scalar_t>(
            xyz, j, k, l, m, nullptr, nullptr, nullptr, nullptr, nullptr);
        const double delta = static_cast<double>(cdelta_ptr[tmap]);
        const int64_t size = csize_ptr[tmap];

        const double ang_a = wrap_0_2pi(phi + kTwoPi);
        const double ang_b = wrap_0_2pi(psi + kTwoPi);
        const double u = ang_a / delta;
        const double v = ang_b / delta;
        int64_t su = static_cast<int64_t>(std::floor(u));
        int64_t sv = static_cast<int64_t>(std::floor(v));
        su = std::min(su, size - 1);
        sv = std::min(sv, size - 1);
        const double da = u - static_cast<double>(su);
        const double db = v - static_cast<double>(sv);

        const int64_t patch = su + size * sv;
        const int64_t coeff_row = coff_ptr[tmap] + patch;
        const scalar_t* coeff = ccoef_ptr + 16 * coeff_row;

        const double db2 = db * db;
        const double da2 = da * da;
        const double da3 = da2 * da;

        double ppoly[4];
        double pd[4];
        for (int r = 0; r < 4; ++r) {
          const double c0 = static_cast<double>(coeff[4 * r + 0]);
          const double c1 = static_cast<double>(coeff[4 * r + 1]);
          const double c2 = static_cast<double>(coeff[4 * r + 2]);
          const double c3 = static_cast<double>(coeff[4 * r + 3]);
          ppoly[r] = c0 + c1 * db + c2 * db2 + c3 * db2 * db;
          pd[r] = c1 + 2.0 * c2 * db + 3.0 * c3 * db2;
        }

        const double e_val = ppoly[0] + ppoly[1] * da + ppoly[2] * da2 + ppoly[3] * da3;
        ec_local += e_val;

        const double dE_dda = ppoly[1] + 2.0 * ppoly[2] * da + 3.0 * ppoly[3] * da2;
        const double dE_ddb = pd[0] + pd[1] * da + pd[2] * da2 + pd[3] * da3;
        const double dE_dphi = dE_dda / delta;
        const double dE_dpsi = dE_ddb / delta;

        accumulate_dihedral_force<scalar_t>(xyz, i, j, k, l, dE_dphi, fl);
        accumulate_dihedral_force<scalar_t>(xyz, j, k, l, m, dE_dpsi, fl);
      }
      ec_tls[static_cast<size_t>(tid)] += ec_local;
    }

    double eb = 0.0;
    double ea = 0.0;
    double ed = 0.0;
    double ec = 0.0;
    for (int tid = 0; tid < nthreads; ++tid) {
      eb += eb_tls[static_cast<size_t>(tid)];
      ea += ea_tls[static_cast<size_t>(tid)];
      ed += ed_tls[static_cast<size_t>(tid)];
      ec += ec_tls[static_cast<size_t>(tid)];
    }

    for (int64_t idx = 0; idx < stride; ++idx) {
      double acc = 0.0;
      for (int tid = 0; tid < nthreads; ++tid) {
        acc += static_cast<double>(
            force_tls[static_cast<size_t>(tid) * static_cast<size_t>(stride) + static_cast<size_t>(idx)]);
      }
      f[idx] = static_cast<scalar_t>(acc);
    }

    e_bond.fill_(static_cast<scalar_t>(eb));
    e_angle.fill_(static_cast<scalar_t>(ea));
    e_dihed.fill_(static_cast<scalar_t>(ed));
    e_cmap.fill_(static_cast<scalar_t>(ec));
  });

  return {e_bond, e_angle, e_dihed, e_cmap, force};
}

}  // namespace

std::vector<torch::Tensor> bonded_energy_force(
    const torch::Tensor& coords,
    const torch::Tensor& bond_i,
    const torch::Tensor& bond_j,
    const torch::Tensor& bond_k,
    const torch::Tensor& bond_r0,
    const torch::Tensor& angle_i,
    const torch::Tensor& angle_j,
    const torch::Tensor& angle_k,
    const torch::Tensor& angle_k0,
    const torch::Tensor& angle_t0,
    const torch::Tensor& dihed_i,
    const torch::Tensor& dihed_j,
    const torch::Tensor& dihed_k,
    const torch::Tensor& dihed_l,
    const torch::Tensor& dihed_force,
    const torch::Tensor& dihed_period,
    const torch::Tensor& dihed_phase,
    const torch::Tensor& cmap_type,
    const torch::Tensor& cmap_i,
    const torch::Tensor& cmap_j,
    const torch::Tensor& cmap_k,
    const torch::Tensor& cmap_l,
    const torch::Tensor& cmap_m,
    const torch::Tensor& cmap_size,
    const torch::Tensor& cmap_delta,
    const torch::Tensor& cmap_offset,
    const torch::Tensor& cmap_coeff) {
  if (!coords.device().is_cpu()) {
    TORCH_CHECK(false, "bonded_energy_force currently supports CPU tensors only");
  }
  return bonded_energy_force_cpu(
      coords,
      bond_i,
      bond_j,
      bond_k,
      bond_r0,
      angle_i,
      angle_j,
      angle_k,
      angle_k0,
      angle_t0,
      dihed_i,
      dihed_j,
      dihed_k,
      dihed_l,
      dihed_force,
      dihed_period,
      dihed_phase,
      cmap_type,
      cmap_i,
      cmap_j,
      cmap_k,
      cmap_l,
      cmap_m,
      cmap_size,
      cmap_delta,
      cmap_offset,
      cmap_coeff);
}

PYBIND11_MODULE(TORCH_EXTENSION_NAME, m) {
  m.doc() = "hessian_ff native bonded extension";
  m.def(
      "bonded_energy_force",
      &bonded_energy_force,
      "Compute bonded energy/force (bond/angle/dihedral/cmap)");
}
