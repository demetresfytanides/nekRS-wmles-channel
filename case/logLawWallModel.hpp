// logLawWallModel.hpp
//
// Smooth-wall log-law wall model, evaluated per boundary node. Host-only,
// header-only, no nekRS/MPI/device dependency (unit test: tests/test_wall_model.cpp).
//
// The mesh boundary sits a distance zInput above a "virtual wall"; the layer
// between the virtual wall and the mesh boundary is not simulated. At every
// boundary node i:
//
//   u_t   = u - (u.n) n                   tangential velocity at the node
//   u_tau : solves |u_t| = u_tau * U+(zInput*u_tau/nu),  U+(y+) = ln(y+)/kappa + B
//   tr    = -u_tau^2 * u_t / |u_t|        traction vector (zero when |u_t| == 0)
//
// tr is a global-frame vector; the device BC (channel.oudf) projects it onto
// the kernel's own tangent basis, so nothing here assumes which axis is the
// wall normal (n comes from the mesh geometric factors).
//
// Traction block handed to the device: three component blocks of length
// fieldOffset, entry [c*fieldOffset + idxVol], read by the BC kernel at
// bc->idxVol. Entries at non-boundary nodes are 0 and never read.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

// Inside a nekRS .udf, dfloat comes from udf.hpp (double in the default
// build). Standalone builds (the unit test) define WM_STANDALONE first.
#ifdef WM_STANDALONE
using dfloat = double;
#endif

namespace wm {

struct LogLawParams {
  dfloat kappa = dfloat(0.41);
  dfloat B = dfloat(5.2);

  void validate() const
  {
    if (!(kappa > dfloat(0))) {
      throw std::runtime_error("wm::LogLawParams: kappa must be > 0, got " + std::to_string(kappa));
    }
  }
};

struct LogLaw {
  LogLawParams params;

  void setup(const LogLawParams &p)
  {
    p.validate();
    params = p;
  }

  // U+ = ln(y+)/kappa + B (log region only; no sublayer blending).
  dfloat uPlusOf(dfloat yPlus) const
  {
    if (!(yPlus > dfloat(0))) {
      throw std::runtime_error("wm::LogLaw::uPlusOf: yPlus must be > 0, got " + std::to_string(yPlus));
    }
    return (dfloat(1) / params.kappa) * std::log(yPlus) + params.B;
  }

  // Solve Umag = uTau * U+(y*uTau/nu) for uTau by Newton iteration.
  dfloat invertUTau(dfloat Umag, dfloat y, dfloat nu, int maxIter = 50, dfloat tol = dfloat(1e-10)) const
  {
    if (!(y > dfloat(0))) {
      throw std::runtime_error("wm::LogLaw::invertUTau: y must be > 0, got " + std::to_string(y));
    }
    if (!(nu > dfloat(0))) {
      throw std::runtime_error("wm::LogLaw::invertUTau: nu must be > 0, got " + std::to_string(nu));
    }
    if (!(Umag >= dfloat(0))) {
      throw std::runtime_error("wm::LogLaw::invertUTau: Umag must be >= 0, got " + std::to_string(Umag));
    }
    if (Umag == dfloat(0)) {
      return dfloat(0);
    }

    // Start from the viscous-sublayer estimate u_tau = sqrt(U nu / y).
    dfloat uTau = std::sqrt(Umag * nu / y);
    if (!(uTau > dfloat(0))) {
      uTau = dfloat(1e-6) * std::max(Umag, dfloat(1));
    }

    for (int it = 0; it < maxIter; ++it) {
      const dfloat yPlus = uTau * y / nu;
      const dfloat uPlus = uPlusOf(yPlus);
      const dfloat F = uTau * uPlus - Umag;
      const dfloat dF = uPlus + dfloat(1) / params.kappa;
      const dfloat step = F / dF;
      dfloat uNew = uTau - step;
      if (!(uNew > dfloat(0))) {
        uNew = uTau * dfloat(0.5);
      }
      const dfloat denom = std::max(uNew, dfloat(1e-12));
      const bool converged = std::fabs(uNew - uTau) <= tol * denom;
      uTau = uNew;
      if (converged) {
        break;
      }
    }
    return uTau;
  }
};

// One boundary node: its volume-array index and the unit surface normal there.
struct BoundaryNode {
  std::size_t idxVol = 0;
  dfloat nx = 0, ny = 0, nz = 0;
};

struct NodeResult {
  dfloat uTau = 0;
  dfloat uTangMag = 0;
  dfloat tr[3] = {0, 0, 0};
};

inline NodeResult nodeTraction(const LogLaw &law, dfloat zInput, dfloat nu, dfloat ux, dfloat uy, dfloat uz,
                               dfloat nx, dfloat ny, dfloat nz)
{
  const dfloat nmag = std::sqrt(nx * nx + ny * ny + nz * nz);
  if (!(nmag > dfloat(0))) {
    throw std::runtime_error("wm::nodeTraction: zero-length surface normal");
  }
  nx /= nmag;
  ny /= nmag;
  nz /= nmag;
  const dfloat un = ux * nx + uy * ny + uz * nz;
  const dfloat tx = ux - un * nx, ty = uy - un * ny, tz = uz - un * nz;
  NodeResult r;
  r.uTangMag = std::sqrt(tx * tx + ty * ty + tz * tz);
  r.uTau = law.invertUTau(r.uTangMag, zInput, nu);
  if (r.uTangMag > dfloat(0)) {
    const dfloat s = -r.uTau * r.uTau / r.uTangMag;
    r.tr[0] = s * tx;
    r.tr[1] = s * ty;
    r.tr[2] = s * tz;
  }
  return r;
}

// Fill the 3*fieldOffset traction block (and optional per-node uTau / |u_t|,
// one entry per node in `nodes` order) from the volume velocity arrays.
inline void fillTraction(const LogLaw &law, dfloat zInput, dfloat nu, const std::vector<dfloat> &ux,
                         const std::vector<dfloat> &uy, const std::vector<dfloat> &uz,
                         const std::vector<BoundaryNode> &nodes, std::size_t fieldOffset,
                         std::vector<dfloat> &trBlock, std::vector<dfloat> *uTauOut = nullptr,
                         std::vector<dfloat> *uTangOut = nullptr)
{
  trBlock.assign(3 * fieldOffset, dfloat(0));
  if (uTauOut) {
    uTauOut->assign(nodes.size(), dfloat(0));
  }
  if (uTangOut) {
    uTangOut->assign(nodes.size(), dfloat(0));
  }
  for (std::size_t k = 0; k < nodes.size(); ++k) {
    const BoundaryNode &nd = nodes[k];
    if (nd.idxVol >= fieldOffset || nd.idxVol >= ux.size() || nd.idxVol >= uy.size() || nd.idxVol >= uz.size()) {
      throw std::runtime_error("wm::fillTraction: node index outside the velocity / traction arrays");
    }
    const NodeResult r =
        nodeTraction(law, zInput, nu, ux[nd.idxVol], uy[nd.idxVol], uz[nd.idxVol], nd.nx, nd.ny, nd.nz);
    for (int c = 0; c < 3; ++c) {
      trBlock[c * fieldOffset + nd.idxVol] = r.tr[c];
    }
    if (uTauOut) {
      (*uTauOut)[k] = r.uTau;
    }
    if (uTangOut) {
      (*uTangOut)[k] = r.uTangMag;
    }
  }
}

// Sample statistics for the per-wall diagnostics (linear-interpolated
// percentiles, population standard deviation).
struct SampleStats {
  std::size_t n = 0;
  dfloat mean = 0, stddev = 0, min = 0, max = 0, p5 = 0, p95 = 0;
};

inline dfloat sampleQuantile(const std::vector<dfloat> &sorted, dfloat q)
{
  if (sorted.empty()) {
    return dfloat(0);
  }
  const dfloat pos = q * dfloat(sorted.size() - 1);
  const std::size_t lo = static_cast<std::size_t>(std::floor(pos));
  const std::size_t hi = std::min(lo + 1, sorted.size() - 1);
  const dfloat w = pos - dfloat(lo);
  return (dfloat(1) - w) * sorted[lo] + w * sorted[hi];
}

inline SampleStats sampleStats(std::vector<dfloat> v)
{
  SampleStats s;
  s.n = v.size();
  if (v.empty()) {
    return s;
  }
  std::sort(v.begin(), v.end());
  dfloat sum = 0;
  for (dfloat x : v) {
    sum += x;
  }
  s.mean = sum / dfloat(v.size());
  dfloat var = 0;
  for (dfloat x : v) {
    var += (x - s.mean) * (x - s.mean);
  }
  s.stddev = std::sqrt(var / dfloat(v.size()));
  s.min = v.front();
  s.max = v.back();
  s.p5 = sampleQuantile(v, dfloat(0.05));
  s.p95 = sampleQuantile(v, dfloat(0.95));
  return s;
}

} // namespace wm
