# nekRS wall-modelled channel at Re_τ = 5200

This is a wall-modelled large-eddy simulation (WMLES) of turbulent plane channel flow at friction Reynolds number Re_τ = 5200, written for [nekRS](https://github.com/Nek5000/nekRS) v26.0.

- The channel is periodic in the streamwise and spanwise directions and is driven by a constant mean pressure gradient.
  - So the turbulence sustains itself. There is no inflow, no recycling and no external forcing other than the pressure gradient.
- The near-wall layer is not resolved.
  - The mesh stops at y⁺ = 100 above two *virtual walls*.
  - At each boundary node, a log-law wall model supplies the wall shear stress as a traction boundary condition.

This case was taken out of a larger sediment-transport project. Everything related to sediment, beds, the moving mesh and the case generator has been removed. What is left is the flow solver setup that ran there as job 9045014 (see [Status](#6-status-what-has-been-run-and-what-it-showed)).

---

## 1. What is simulated

| Quantity | Value |
|---|---|
| Flow | incompressible, constant density and viscosity |
| Domain | Lx × Ly × Lz = 4 × (2 − 2z₁) × 2, in units of the half-height h |
| Virtual walls | y = 0 and y = 2 (not part of the mesh) |
| Mesh boundaries | y = z₁ and y = 2 − z₁, z₁ = 100/5200 = 0.019231 (y⁺ = 100) |
| Periodicity | x and z |
| Forcing | constant G = −dp/dx = ρu_τ²/h = 1 |
| Re_τ | h u_τ / ν = 5200 (ν = 1/5200) |
| Elements | 16 × 20 × 8 = 2560 hexahedra, polynomial order 7 (8 GLL points per direction) |
| Mean spacing | Δx⁺ ≈ 186, Δz⁺ ≈ 186 (element size / 7) |
| Wall-normal spacing | first element 0.0483h (251 wall units); midplane element 0.170h (887) |
| Wall-normal stretching | geometric, ratio 1.15 per element, mirror-symmetric (`scripts/wall_stretch.py`) |
| Time step | Δt = 1×10⁻⁴ h/u_τ, 15 000 steps (t = 1.5 h/u_τ) |

**Units.** Lengths are in h (measured to the virtual walls), velocities in the nominal u_τ and density ρ = 1. With these units ν = 1/Re_τ, and the forcing −dp/dx = 1 fixes the virtual-wall friction velocity at 1 in a statistically steady state.

**Why "self-sustained".**
- In a periodic channel with a fixed pressure gradient, turbulence feeds on the mean shear that the forcing maintains. No inflow turbulence has to be generated.
- The mean flow adjusts until the wall stress balances the pressure gradient.
- Unlike a constant-flow-rate controller, the pressure gradient fixes u_τ (and so Re_τ), not the bulk velocity.

---

## 2. How the problem is set up

### 2.1 Forcing
`userf()` in `case/channel.udf` adds the constant body force −dp/dx = ρ·1²/h = 1 to the x-momentum equation every step, by filling the explicit forcing array `nrs->fluid->o_EXT`. Do not set `constFlowRate` in the `.par` file; the UDF stops with an error if it is set.

### 2.2 Boundary conditions
- **x and z:** periodic.
- **The two mesh boundaries (y = z₁ and y = 2 − z₁):** `zeroDirichletN/udfNeumann`.
  - The wall-normal velocity is zero (no penetration).
  - The two tangential components get a prescribed traction (shear stress), supplied per node by the wall model (§3).
  - Pressure is homogeneous Neumann.
- **Boundary ids:** assigned in `case/channel.usr`. Both walls carry the code `W` in `channel.box` and are told apart by face number:
  - id 1 = top boundary, y = 1.98;
  - id 2 = bottom boundary, y = 0.019.

### 2.3 Initial condition (fresh start only; skipped on restart)
- **Mean profile.** Measured from the nearest virtual wall, with u_τ = 1:
  - U⁺ = y⁺ below y⁺ = 11.53;
  - U⁺ = ln(y⁺)/κ + B above that.
- **Perturbation.** One divergence-free roll/streak pair, **u′ = ∇ × A**, with:
  - A = (g(y) sin βz, 0, g(y) sin αx), g(y) = [(y − y_lo)(y_hi − y)]²;
  - α = 2π/Lx, β = 2π/Lz;
  - amplitude = `ic_amp` × the centreline velocity (`ic_amp` = 0.024228).
- **Properties of the perturbation:**
  - It is exactly divergence-free.
  - It vanishes, with its derivative, at the virtual walls.
- **Symmetry caveat.** This IC is exactly symmetric under (x, y, z) → (x + Lx/2, 2 − y, z + Lz/2), v → −v. The solver keeps that symmetry until round-off breaks it.
  - In job 9045014 the two walls were mirror images of each other to output precision until t ≈ 0.8. They only began to evolve independently around t ≈ 1.1.
  - If you need the walls to be independent early, add a non-symmetric perturbation.

### 2.4 Numerics (from `case/channel.par`)
- **Time stepping:** `tombo2` (second order) with characteristic subcycling (`advectionSubCyclingSteps = 8`).
- **Subgrid model:** none explicit. The high-pass-filter relaxation `regularization = hpfrt + nModes=1 + scalingCoeff=10` removes energy from the highest polynomial mode and acts as the subgrid dissipation.
- **Solver tolerances:** pressure 1e-6, velocity 1e-8.

---

## 3. How the shear velocity is computed

The wall model runs on the host every time step (`wallModelStep()` in `case/channel.udf`). The maths is in `case/logLawWallModel.hpp`. The sequence is:

**Step 0, setup (once, `wallModelSetup()`).** For every face node on boundary ids 1 and 2, store:
- its index in the volume arrays (`vmapM`);
- its unit outward normal **n** (from the mesh surface geometric factors `sgeo`);
- its surface quadrature weight times Jacobian (`WsJ`);
- its coordinates.

So the wall normal is taken from the mesh. The code never assumes it is ±y.

**Step 1, sample the velocity.** Copy the current velocity field (the solution at the end of the previous step) from the device to the host.

**Step 2, tangential velocity at the boundary node itself.** At each boundary node:

$$\mathbf{u}_t = \mathbf{u} - (\mathbf{u}\cdot\mathbf{n})\,\mathbf{n}, \qquad U_t = \lvert \mathbf{u}_t \rvert .$$

- The input is the node's own instantaneous velocity.
- There is no plane average, no sampling at a different height and no time filtering.
- So u_τ varies from node to node and from step to step.

**Step 3, invert the log law.** The node sits at height z₁ = 100/5200 above its virtual wall (`wm_z_input` in the `.par`). Solve the smooth-wall log law

$$U_t = u_\tau \, U^+\!\left(\frac{z_1 u_\tau}{\nu}\right), \qquad U^+(y^+) = \frac{1}{\kappa}\ln y^+ + B, \qquad \kappa = 0.41,\; B = 5.2,$$

for u_τ by Newton iteration (`wm::LogLaw::invertUTau`):

| Newton step | Formula |
|---|---|
| initial guess | $u_\tau^{0} = \sqrt{U_t\,\nu / z_1}$ (the viscous-sublayer estimate) |
| residual | $F(u_\tau) = u_\tau\, U^+(z_1 u_\tau/\nu) - U_t$ |
| derivative | $F'(u_\tau) = U^+(z_1 u_\tau/\nu) + 1/\kappa$ |
| update | $u_\tau \leftarrow u_\tau - F/F'$; if that would give $u_\tau \le 0$, use $u_\tau \leftarrow u_\tau/2$ instead |
| stop | relative change $\le 10^{-10}$, or 50 iterations |

- If U_t = 0, then u_τ = 0.
- Only the log branch is used. With z₁⁺ = 100 nominal, the boundary node is always in the log layer unless u_τ falls below about 0.12. Job 9045014 did have individual nodes down to u_τ ≈ 0.005; the log law is still applied there.
- The unit test (`tests/test_wall_model.cpp`) checks that the inversion recovers u_τ to 1e-12 for u_τ from 0.05 to 2.5.

**Step 4, traction.** The wall shear stress is applied opposite the local tangential velocity:

$$\boldsymbol{\tau} = -\,u_\tau^{2}\,\frac{\mathbf{u}_t}{U_t} \qquad (\text{kinematic},\ \rho = 1;\ \lvert\boldsymbol{\tau}\rvert = u_\tau^{2}).$$

- The three global components of τ are written into `o_usrwrk` as three blocks of length `fieldOffset`, at each node's volume index.
- On the device, `udfNeumann` in `case/channel.oudf` projects τ onto the kernel's local tangent vectors t₁ and t₂ (`tr1 = τ·t₁`, `tr2 = τ·t₂`).
- In the first 4 steps the block is read back from the device and compared with the host copy. The run stops if they differ.

**Time lag.** The traction used during step n+1 comes from the velocity at the end of step n (explicit, one step behind). At step 0 it comes from the initial condition.

### 3.1 The u_τ values printed in the log
Every `diag_every` steps (default 100), each wall prints a `[logbnd-wall]` line. It contains two different u_τ measures.

- **`uTau_mean` (also `_std`, `_min`, `_max`, `_p5`, `_p95`):** statistics over the per-node u_τ from step 3.
  - Nodes shared between elements are counted once, matched by (x, z).
  - The periodic image lines (x = 0 and x = Lx, z = 0 and z = Lz) have different coordinates, so they are counted twice. This leaves 6441 nodes per wall instead of 6272.
  - So `uTau_mean` is an **unweighted** node mean, slightly biased towards the domain edges.
- **`Fx`, `Fy`, `Fz`, `A`:** the delivered traction integrated over the wall with the surface quadrature weights, Σ WsJ·τ, and the wall area A = Lx·Lz = 8.
  - **u_τ,F = √(−Fx/A)** is the area-weighted friction velocity actually delivered to the flow.
  - `scripts/extract_log.py` writes it as `utauF1` and `utauF2`.

**What value to expect.** The simulated region ends at y = z₁, not at the virtual wall. In a statistically steady state the traction at the mesh boundary balances the pressure gradient over the simulated height only:

$$2\,\tau_b = G\,(2 - 2z_1) \;\Rightarrow\; \tau_b = 1 - z_1 = 0.981, \qquad u_{\tau,F} = \sqrt{1 - z_1} = 0.990,$$

with G = −dp/dx = 1 the imposed pressure gradient.

The full-height value u_τ = 1 belongs to the virtual wall.

---

## 4. Repository layout

```
case/
  channel.par            nekRS parameter file (numerics, BC types, [CASEDATA])
  channel.udf            forcing, initial condition, wall model glue, diagnostics
  channel.oudf           device BC: projects the per-node traction onto t1/t2
  channel.usr            assigns boundary ids 1 (top) and 2 (bottom)
  channel.box            genbox mesh definition
  logLawWallModel.hpp    log law, u_tau inversion, per-node traction (host, header-only)
  udf.cmake              puts case/ on the UDF include path
scripts/
  make_mesh.sh           genbox + genmap -> channel.re2 / channel.ma2
  wall_stretch.py        wall-normal element boundaries for channel.box
  run_aurora.sh          PBS script for ALCF Aurora (build, then run, 12 ranks on 1 node)
  extract_log.py         log -> CSV time series (U_b, u_tau per wall, momentum budget)
tests/
  test_wall_model.cpp    unit test of logLawWallModel.hpp (no nekRS needed)
```

---

## 5. How to run

**Requirements:**
- nekRS v26.0;
- Nek5000's `genbox` and `genmap` (from a Nek5000 build: `cd Nek5000/tools && ./maketools genbox genmap`);
- Python 3 with numpy, for the helper scripts;
- a C++17 compiler for the unit test.

```bash
# 0. unit test of the wall model (seconds, any machine)
g++ -std=c++17 -O2 -I case tests/test_wall_model.cpp -o test_wall_model && ./test_wall_model

# 1. make a run directory and the mesh
mkdir run && cp case/* run/ && cd run
NEK_TOOLS=/path/to/Nek5000/bin ../scripts/make_mesh.sh

# 2a. Aurora: edit the account, then submit from the run directory
qsub -A <project> ../scripts/run_aurora.sh

# 2b. elsewhere (example: CUDA, 4 ranks)
mpirun -np 4 $NEKRS_HOME/bin/nekrs --setup channel --backend cuda

# 3. post-process the log
python3 ../scripts/extract_log.py <logfile> --out timeseries.csv
```

**Settings and outputs:**
- **Tunable parameters** are in `[CASEDATA]` in `channel.par`: κ, B, z₁, the virtual-wall positions, the IC amplitude and the diagnostics interval.
- **If you change the mesh extent,** keep `wm_z_input`, `vwall_lo` and `vwall_hi` consistent with it. The UDF checks that the mesh boundaries lie exactly z₁ inside the virtual walls.
- **Field files** (`channel0.f*`) are written every `writeInterval` = 1000 steps.
- **Restarts:** set `startFrom` in `[GENERAL]`. The IC is skipped when a restart file is given.

---

## 6. Status: what has been run and what it showed

**This setup is not yet a validated WMLES.** It runs stably, but in the one long run so far the wall stress and bulk velocity had not settled.

### 6.1 Job 9045014 (sediment-project copy of this case)
The run used the same mesh, parameters and equations as this repo: 12 ranks on Aurora, 15 000 steps, t = 0 → 1.5.

- **Stability.** No NaN or abort. CFL ranged 0.19–0.87, max |u| 26–77.
- **u_τ.**
  - The mean of both walls falls from 1.0 to about 0.40 at t = 0.7.
  - It rises during t ≈ 0.9–1.1 and stays at 0.69–0.72 from t = 1.1 to 1.5 (u_τ,F 0.665–0.697).
  - The steady-state expectation is u_τ,F = 0.990 (§3.1). The delivered traction is about half of what the pressure gradient supplies.
- **Bulk velocity.** Because of that, the bulk velocity was still rising at t = 1.5: U_b over the simulated region was 24.89, with dU_b/dt = +0.54.
  - Lee & Moser (2015) DNS gives 24.31 over the same y range and 24.10 over the full channel.
  - The momentum budget H dU_b/dt = G H + (F₁ + F₂)/A closes to within 3% of G H at every printed time. So the imbalance is physical and not an accounting error.
- **Mean profile at t = 1.4–1.5.** U(y) − U(boundary) against the log law, with y⁺ measured from the mesh boundary, agrees within 0.43–0.63 U⁺ (max |difference| per wall) for y⁺ = 35–370. At t = 1.3 it was below the log law by up to 1.4.
- **Wall symmetry.** See the caveat in §2.3. The walls had evolved independently for only about 0.4 time units by t = 1.5.

These numbers come from the sediment project's run record (`validation_state/records/V6a/logbnd-2026-11-06/` there).

### 6.2 Equivalence of this repo with that run
- **Mesh.** `scripts/make_mesh.sh` produces a `channel.re2` that is byte-identical to the 9045014 mesh, except the 128 bottom-boundary face codes (`mv` → `W`).
- **200-step run on 12 ranks.** Used to check that the trimmed UDF reproduces the original diagnostics (Aurora job 9064489, build and run exit codes 0).
  - **Result: reproduced.** All wall-model lines at steps 0–3, 100 and 200 match 9045014. This covers setup, grid, sync, `Ubulk_res`, `maxU`, and every u_τ statistic and velocity mean.
  - **Only two differences:**
    - `Fx` differs in its last printed digit, at most 2e-10 (about 6e-11 relative).
    - `Fz` is printed as `-0` in some places and `0` in others.
  - Both differences are consistent with round-off in the parallel sum. Example at step 200, wall 1:
    - `Ubulk_res = 23.82730375`
    - `uTau_mean = 0.66648223`
    - `Fx = -3.5533269963` (original run: `-3.5533269961`)

---

## 7. What was removed from the original project

The original case was generated by the sediment project's case generator. The UDF there was 1840 lines, with many alternative code paths behind `[CASEDATA]` switches. This repo keeps only the path that job 9045014 actually ran.

**Removed:**
- **Sediment-transport and moving-mesh items:**
  - bed tags and bed-shape and gravity settings;
  - the dependency on the sediment library (`src/`).
- **Wall-model variants not used in the logbnd run:**
  - plane-averaged sampling at an interior exchange height;
  - the u_τ time filter (EMA);
  - the prescribed-constant-traction test;
  - the flux diagnostic;
  - time averaging (`tavg`);
  - the symmetric-IC probe;
  - generator bookkeeping.

**Unchanged:**
- the equations;
- the parameter values;
- the IC formula;
- the per-node traction rule;
- the format of the diagnostic lines.
