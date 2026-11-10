// Unit test for case/logLawWallModel.hpp (no nekRS needed).
//   g++ -std=c++17 -O2 -I case tests/test_wall_model.cpp -o test_wall_model && ./test_wall_model
#define WM_STANDALONE
#include "logLawWallModel.hpp"

#include <cstdio>
#include <cstdlib>

static int fails = 0;
static void check(bool ok, const char *what, double got, double want)
{
  std::printf("%s  %-55s got %.15g  want %.15g\n", ok ? "PASS" : "FAIL", what, got, want);
  if (!ok) {
    ++fails;
  }
}

int main()
{
  wm::LogLaw law;
  law.setup(wm::LogLawParams{0.41, 5.2});
  const double nu = 1.0 / 5200.0, z = 100.0 / 5200.0;

  // 1. Inversion round trip: U = u_tau * (ln(z u_tau/nu)/kappa + B) -> u_tau.
  for (double ut : {0.05, 0.4, 0.71, 1.0, 2.5}) {
    const double U = ut * law.uPlusOf(z * ut / nu);
    const double got = law.invertUTau(U, z, nu);
    char msg[96];
    std::snprintf(msg, sizeof msg, "invertUTau round trip, u_tau=%.2f", ut);
    check(std::fabs(got - ut) <= 1e-12 * ut, msg, got, ut);
  }
  // 2. Reference value: u_tau = 1 at y+ = 100 needs U = ln(100)/0.41 + 5.2.
  check(std::fabs(law.uPlusOf(100.0) - (std::log(100.0) / 0.41 + 5.2)) < 1e-14, "U+(100)", law.uPlusOf(100.0),
        std::log(100.0) / 0.41 + 5.2);

  // 3. Traction: magnitude u_tau^2, opposite the tangential velocity, normal
  //    velocity ignored, for a tilted (non-axis) normal.
  const double n[3] = {0.3, 0.9, 0.1};
  const double nm = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
  const double u[3] = {14.0, 5.0, -3.0};
  const auto r = wm::nodeTraction(law, z, nu, u[0], u[1], u[2], n[0], n[1], n[2]);
  const double un = (u[0] * n[0] + u[1] * n[1] + u[2] * n[2]) / nm;
  double t[3], tm = 0, trm = 0, trn = 0, dot = 0;
  for (int c = 0; c < 3; ++c) {
    t[c] = u[c] - un * n[c] / nm;
    tm += t[c] * t[c];
    trm += r.tr[c] * r.tr[c];
    trn += r.tr[c] * n[c] / nm;
  }
  tm = std::sqrt(tm);
  trm = std::sqrt(trm);
  for (int c = 0; c < 3; ++c) {
    dot += r.tr[c] * t[c];
  }
  check(std::fabs(r.uTangMag - tm) < 1e-12 * tm, "|u_t| removes the normal component", r.uTangMag, tm);
  check(std::fabs(trm - r.uTau * r.uTau) < 1e-12, "|tr| = u_tau^2", trm, r.uTau * r.uTau);
  check(std::fabs(trn) < 1e-12, "tr is tangential (tr.n = 0)", trn, 0.0);
  check(std::fabs(dot / (trm * tm) + 1.0) < 1e-12, "tr opposite u_t (cosine = -1)", dot / (trm * tm), -1.0);
  check(std::fabs(r.uTau - law.invertUTau(tm, z, nu)) < 1e-12, "u_tau from |u_t| only", r.uTau,
        law.invertUTau(tm, z, nu));

  // 4. Zero velocity -> zero traction.
  const auto r0 = wm::nodeTraction(law, z, nu, 0, 0, 0, 0, 1, 0);
  check(r0.uTau == 0 && r0.tr[0] == 0 && r0.tr[1] == 0 && r0.tr[2] == 0, "zero velocity -> zero traction", r0.uTau, 0);

  // 5. Block layout: entry [c*fieldOffset + idxVol].
  std::vector<double> ux(10, 0), uy(10, 0), uz(10, 0), blk;
  ux[7] = 16.0;
  std::vector<wm::BoundaryNode> nodes(1);
  nodes[0].idxVol = 7;
  nodes[0].ny = -1;
  wm::fillTraction(law, z, nu, ux, uy, uz, nodes, 10, blk);
  const double ut16 = law.invertUTau(16.0, z, nu);
  check(blk.size() == 30 && std::fabs(blk[7] + ut16 * ut16) < 1e-14 && blk[17] == 0 && blk[27] == 0,
        "fillTraction block layout", blk[7], -ut16 * ut16);

  std::printf("%s (%d failures)\n", fails ? "FAILED" : "ALL PASSED", fails);
  return fails ? EXIT_FAILURE : EXIT_SUCCESS;
}
