#include <torch/extension.h>

#include <vector>

using torch::indexing::Slice;

namespace {

constexpr double kCoulomb = 332.0637132991921;

void scatter_pair_blocks_flat(
    torch::Tensor& h_flat,
    const torch::Tensor& row_atoms,
    const torch::Tensor& col_atoms,
    const torch::Tensor& blocks,
    int64_t ndof) {
  if (row_atoms.numel() == 0) {
    return;
  }
  auto offs = torch::arange(3, row_atoms.options());
  auto rows = row_atoms.unsqueeze(1) * 3 + offs.unsqueeze(0);  // [P,3]
  auto cols = col_atoms.unsqueeze(1) * 3 + offs.unsqueeze(0);  // [P,3]
  auto lin = rows.unsqueeze(2) * ndof + cols.unsqueeze(1);     // [P,3,3]
  h_flat.index_add_(0, lin.reshape({-1}), blocks.reshape({-1}));
}

torch::Tensor central_pair_blocks(
    const torch::Tensor& rij,
    const torch::Tensor& dV_dr,
    const torch::Tensor& d2V_dr2) {
  auto r2 = (rij * rij).sum(-1).clamp_min(1.0e-24);
  auto inv_r = torch::rsqrt(r2);
  auto u = rij * inv_r.unsqueeze(-1);
  auto uu = u.unsqueeze(-1) * u.unsqueeze(-2);
  auto eye = torch::eye(3, rij.options()).unsqueeze(0);
  auto b = dV_dr * inv_r;
  auto a = d2V_dr2 - b;
  return a.unsqueeze(-1).unsqueeze(-1) * uu + b.unsqueeze(-1).unsqueeze(-1) * eye;
}

torch::Tensor scatter_pair_blocks_impl(
    const torch::Tensor& blocks,
    const torch::Tensor& ii,
    const torch::Tensor& jj,
    const torch::Tensor& active_map,
    int64_t ndof) {
  auto li = active_map.index_select(0, ii);
  auto lj = active_map.index_select(0, jj);
  auto li_valid = li >= 0;
  auto lj_valid = lj >= 0;
  auto both = li_valid & lj_valid;

  auto h2 = torch::zeros({ndof, ndof}, blocks.options());
  auto h_flat = h2.reshape({ndof * ndof});

  if (li_valid.any().item<bool>()) {
    auto li_idx = li.index({li_valid});
    auto blk = blocks.index({li_valid});
    scatter_pair_blocks_flat(h_flat, li_idx, li_idx, blk, ndof);
  }
  if (lj_valid.any().item<bool>()) {
    auto lj_idx = lj.index({lj_valid});
    auto blk = blocks.index({lj_valid});
    scatter_pair_blocks_flat(h_flat, lj_idx, lj_idx, blk, ndof);
  }
  if (both.any().item<bool>()) {
    auto li_idx = li.index({both});
    auto lj_idx = lj.index({both});
    auto blk = -blocks.index({both});
    scatter_pair_blocks_flat(h_flat, li_idx, lj_idx, blk, ndof);
    scatter_pair_blocks_flat(h_flat, lj_idx, li_idx, blk, ndof);
  }
  return h2;
}

}  // namespace

torch::Tensor bond_hessian(
    const torch::Tensor& coords,
    const torch::Tensor& bond_i,
    const torch::Tensor& bond_j,
    const torch::Tensor& bond_k,
    const torch::Tensor& bond_r0,
    const torch::Tensor& active_map,
    int64_t ndof) {
  if (bond_i.numel() == 0) {
    return torch::zeros({ndof, ndof}, coords.options());
  }

  auto li = active_map.index_select(0, bond_i);
  auto lj = active_map.index_select(0, bond_j);
  auto keep = (li >= 0) | (lj >= 0);
  if (!keep.any().item<bool>()) {
    return torch::zeros({ndof, ndof}, coords.options());
  }

  auto ii = bond_i.index({keep});
  auto jj = bond_j.index({keep});
  auto kk = bond_k.index({keep});
  auto r0 = bond_r0.index({keep});

  auto rij = coords.index_select(0, jj) - coords.index_select(0, ii);
  auto r2 = (rij * rij).sum(-1).clamp_min(1.0e-24);
  auto inv_r = torch::rsqrt(r2);
  auto r = r2 * inv_r;
  auto dr = r - r0;
  auto dV_dr = 2.0 * kk * dr;
  auto d2V_dr2 = 2.0 * kk;
  auto blocks = central_pair_blocks(rij, dV_dr, d2V_dr2);
  return scatter_pair_blocks_impl(blocks, ii, jj, active_map, ndof);
}

torch::Tensor nonbonded_hessian(
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
    const c10::optional<torch::Tensor>& inv_scee_opt,
    const c10::optional<torch::Tensor>& inv_scnb_opt,
    const torch::Tensor& active_map,
    int64_t ndof) {
  if (pair_i.numel() == 0) {
    return torch::zeros({ndof, ndof}, coords.options());
  }

  auto li = active_map.index_select(0, pair_i);
  auto lj = active_map.index_select(0, pair_j);
  auto keep = (li >= 0) | (lj >= 0);
  if (!keep.any().item<bool>()) {
    return torch::zeros({ndof, ndof}, coords.options());
  }

  auto ii = pair_i.index({keep});
  auto jj = pair_j.index({keep});

  c10::optional<torch::Tensor> inv_scee = c10::nullopt;
  c10::optional<torch::Tensor> inv_scnb = c10::nullopt;
  if (inv_scee_opt.has_value()) {
    inv_scee = inv_scee_opt.value().index({keep});
  }
  if (inv_scnb_opt.has_value()) {
    inv_scnb = inv_scnb_opt.value().index({keep});
  }

  auto rij = coords.index_select(0, jj) - coords.index_select(0, ii);
  auto r2 = (rij * rij).sum(-1).clamp_min(1.0e-24);
  auto inv_r = torch::rsqrt(r2);
  auto inv_r2 = inv_r * inv_r;

  auto qq = charge.index_select(0, ii) * charge.index_select(0, jj);
  auto c = qq * static_cast<double>(kCoulomb);
  if (inv_scee.has_value()) {
    c = c * inv_scee.value();
  }
  auto dV_dr = -c * inv_r2;
  auto d2V_dr2 = 2.0 * c * inv_r2 * inv_r;

  auto ti = atom_type.index_select(0, ii);
  auto tj = atom_type.index_select(0, jj);
  auto raw = nb_index.index({ti, tj});

  auto dV_dr_nb = torch::zeros_like(dV_dr);
  auto d2V_dr2_nb = torch::zeros_like(d2V_dr2);

  auto lj_mask = raw > 0;
  if (lj_mask.any().item<bool>()) {
    auto lj_idx = raw.index({lj_mask}) - 1;
    auto a = lj_acoef.index_select(0, lj_idx);
    auto b = lj_bcoef.index_select(0, lj_idx);
    auto x = inv_r.index({lj_mask});
    auto x2 = x * x;
    auto x4 = x2 * x2;
    auto x6 = x4 * x2;
    auto x7 = x6 * x;
    auto x8 = x7 * x;
    auto x13 = x7 * x6;
    auto x14 = x13 * x;
    dV_dr_nb.index_put_({lj_mask}, -12.0 * a * x13 + 6.0 * b * x7);
    d2V_dr2_nb.index_put_({lj_mask}, 156.0 * a * x14 - 42.0 * b * x8);
  }

  auto hb_mask = raw < 0;
  if (hb_mask.any().item<bool>() && hb_acoef.numel() > 0 && hb_bcoef.numel() > 0) {
    auto hb_idx_full = (-raw.index({hb_mask})) - 1;
    auto valid = (hb_idx_full >= 0) & (hb_idx_full < hb_acoef.numel()) & (hb_idx_full < hb_bcoef.numel());
    if (valid.any().item<bool>()) {
      auto hb_pos = torch::nonzero(hb_mask).reshape(-1).index({valid});
      auto hb_idx = hb_idx_full.index({valid});
      auto a = hb_acoef.index_select(0, hb_idx);
      auto b = hb_bcoef.index_select(0, hb_idx);
      auto x = inv_r.index({hb_pos});
      auto x2 = x * x;
      auto x4 = x2 * x2;
      auto x6 = x4 * x2;
      auto x10 = x6 * x4;
      auto x11 = x10 * x;
      auto x12 = x6 * x6;
      auto x13 = x12 * x;
      auto x14 = x13 * x;
      dV_dr_nb.index_put_({hb_pos}, -12.0 * a * x13 + 10.0 * b * x11);
      d2V_dr2_nb.index_put_({hb_pos}, 156.0 * a * x14 - 110.0 * b * x12);
    }
  }

  if (inv_scnb.has_value()) {
    dV_dr_nb = dV_dr_nb * inv_scnb.value();
    d2V_dr2_nb = d2V_dr2_nb * inv_scnb.value();
  }

  dV_dr = dV_dr + dV_dr_nb;
  d2V_dr2 = d2V_dr2 + d2V_dr2_nb;

  auto blocks = central_pair_blocks(rij, dV_dr, d2V_dr2);
  return scatter_pair_blocks_impl(blocks, ii, jj, active_map, ndof);
}

torch::Tensor scatter_local_hessian(
    const torch::Tensor& local_h,
    const torch::Tensor& dof_idx,
    int64_t ndof) {
  auto h2 = torch::zeros({ndof, ndof}, local_h.options());
  if (local_h.numel() == 0 || dof_idx.numel() == 0) {
    return h2;
  }
  TORCH_CHECK(dof_idx.scalar_type() == torch::kInt64, "dof_idx must be int64");
  TORCH_CHECK(local_h.dim() == 3, "local_h must be [T,m,m]");
  TORCH_CHECK(dof_idx.dim() == 2, "dof_idx must be [T,m]");
  TORCH_CHECK(local_h.size(0) == dof_idx.size(0), "local_h/dof_idx term count mismatch");
  TORCH_CHECK(local_h.size(1) == dof_idx.size(1), "local_h/dof_idx local dof mismatch");
  TORCH_CHECK(local_h.size(1) == local_h.size(2), "local_h must be square on local dof axes");

  const int64_t m = dof_idx.size(1);
  auto rows = dof_idx.unsqueeze(2).expand({dof_idx.size(0), m, m});
  auto cols = dof_idx.unsqueeze(1).expand({dof_idx.size(0), m, m});
  auto valid = (rows >= 0) & (cols >= 0);
  if (!valid.any().item<bool>()) {
    return h2;
  }

  auto lin = (rows * ndof + cols).masked_select(valid);
  auto vals = local_h.masked_select(valid);
  auto h_flat = h2.reshape({ndof * ndof});
  h_flat.index_add_(0, lin, vals);
  return h2;
}

PYBIND11_MODULE(TORCH_EXTENSION_NAME, m) {
  m.doc() = "hessian_ff analytical Hessian helper extension";
  m.def("bond_hessian", &bond_hessian, "Bond analytical Hessian");
  m.def("nonbonded_hessian", &nonbonded_hessian, "Nonbonded analytical Hessian");
  m.def("scatter_local_hessian", &scatter_local_hessian, "Scatter local term Hessian blocks");
}
