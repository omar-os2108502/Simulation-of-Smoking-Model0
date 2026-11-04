# Basic Smoking Dynamics Model (MATLAB)

MATLAB implementation of a four-compartment ODE model of smoking dynamics, based on the basic model of Sharomi & Gumel (2008). The code simulates how a population moves between potential smokers, active smokers, temporary quitters, and permanent quitters, and shows how the smoking reproduction number **R₀** decides whether smoking dies out or persists.

## Model

The total population is held constant at `N = P + S + Qt + Qp`.

| Compartment | Meaning |
|---|---|
| `P`  | Potential smokers |
| `S`  | Active smokers |
| `Qt` | Temporary quitters |
| `Qp` | Permanent quitters |

```
dP/dt  = μN − β·P·S/N − μP
dS/dt  = β·P·S/N + α·Qt − (μ + γ)·S
dQt/dt = γ(1 − r)·S − (μ + α)·Qt
dQp/dt = γ·r·S − μ·Qp
```

| Parameter | Description | Value used |
|---|---|---|
| `β` | Smoking transmission (recruitment) rate | 0.03 / 0.30 |
| `μ` | Natural mortality / turnover rate (per year) | 0.05 (= 1/20) |
| `γ` | Quitting rate (per year) | 0.30 |
| `α` | Relapse rate of temporary quitters (per year) | 0.30 |
| `r` | Fraction of quitters who quit permanently | 0.50 |

### Reproduction number

```
R₀ = β(μ + α) / [ μ(μ + α) + γ(αr + μ) ]
```

- **R₀ < 1**: the smoking-free equilibrium `E₀ = (N, 0, 0, 0)` is stable and smoking dies out.
- **R₀ > 1**: `E₀` is unstable and smoking persists at an endemic equilibrium.

With the parameters above:

| Scenario | β | R₀ | Outcome |
|---|---|---|---|
| 1 | 0.03 | ≈ 0.1355 | Eradication |
| 2 | 0.30 | ≈ 1.3548 | Persistence |

## Repository contents

| File | Description |
|---|---|
| `run_model0_simulation.m` | Self-contained, full analysis. Runs both scenarios for 150 years, computes R₀, builds the Jacobian at the smoking-free equilibrium (and at the numerically found endemic equilibrium for Scenario 2), prints the eigenvalues, and plots the trajectories. |
| `simulate_basic_smoking_model.m` | Lighter driver that simulates both scenarios for 60 years and plots the results. Uses the two helper files below. |
| `basic_smoking_ode.m` | Right-hand side of the ODE system, for use with `ode45`. |
| `calculate_Ro.m` | Computes R₀ from `β, μ, α, γ, r`. |

## Usage

Requires MATLAB (no extra toolboxes; only `ode45`, `eig`, and basic plotting). GNU Octave should also work, though this has not been tested.

Full analysis with stability output:

```matlab
run_model0_simulation
```

This prints R₀, the eigenvalues at each equilibrium, and the endemic equilibrium `E* = (P, S, Qt, Qp)`, and opens a figure with the eradication and persistence trajectories side by side.

Quick simulation only:

```matlab
simulate_basic_smoking_model
```

Using the helpers directly:

```matlab
R0 = calculate_Ro(0.30, 0.05, 0.30, 0.30, 0.50);

y0 = [1600; 50; 10; 0];   % [P; S; Qt; Qp]
N  = sum(y0);
[t, y] = ode45(@(t,y) basic_smoking_ode(t, y, 0.30, 0.05, 0.30, 0.30, 0.50, N), [0 60], y0);
plot(t, y);
```

## Initial conditions

`P(0) = 1600`, `S(0) = 50`, `Qt(0) = 10`, `Qp(0) = 0`, giving `N = 1660`.

## Reference

Sharomi, O. & Gumel, A. B. (2008). Curtailing smoking dynamics: A mathematical modeling approach. *Applied Mathematics and Computation*, 195(2), 475–499.

## License

Add a license of your choice (e.g. MIT) here.