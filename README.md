# The Numerical Integrator Sandbox

**[Click Here to Launch the Live Interactive Dashboard](https://numerical-integrator-sandbox.streamlit.app)**

An interactive, educational web application designed to visually demonstrate and compare how different numerical integration schemes handle a classic 2D planar two-body gravitational orbit. 

This sandbox isolates the mathematical behavior of three foundational solvers (**Forward Euler**, **Classical RK4**, and **Symplectic Leapfrog**) by executing them simultaneously under identical initial physics conditions and timestep constraints.

## Why This Tool Exists

In introductory physics and computational astrophysics, understanding why certain algorithms fail while others succeed can be highly abstract. This application turns mathematical truncation errors, step-size instabilities, and energy conservation violations into dramatic, real-time visual events. It is a plug-and-play classroom resource for educators to help students build genuine intuition for numerical analysis.

---

## Core Features

* **Simultaneous Multi-Method Tracking:** Trajectories for all three solvers are rendered on a unified, interactive Plotly canvas to immediately highlight path divergences caused strictly by algorithmic architecture.
* **Relative Energy Error Analysis:** Tracks specific orbital energy ($E = K + U$) on a logarithmic scale over time (`abs((E_t - E_0) / E_0)`), explicitly plotting the exact physical drift signature of each solver.
* **Interactive Physics Sandbox:** Dynamic sidebar sliders allow users to adjust the initial position, tangential velocity, total simulation timeframe, and—most importantly—the numerical timestep size ($dt$).
* **Real-Time Stability Warnings:** Features automatic UI alerts when the selected timestep exceeds stability thresholds, signaling to students precisely where first-order approximations collapse.

---

## The Physics & Integrators Under the Hood

The simulation models a bound two-body planetary orbit around a central fixed stellar mass located at the origin under normalized units where $GM = 1.0$. The initial velocity is configured as purely tangential ($v_x = 0$).

### 1. Forward Euler (1st-Order Explicit)
* **The Math:** `state_next = state + dt * derivatives(state)`
* **Behavior:** Simple and computationally cheap, but accumulates local truncation errors aggressively. Lacks symplectic structure, causing it to systematically inject artificial energy into the system, forcing orbits to rapidly spiral outward.

### 2. Classical RK4 (4th-Order Multi-Stage)
* **The Math:** Evaluates derivatives at four distinct intermediate stages (`k1` through `k4`) across each step, combining them via fixed Taylor-series weights.
* **Behavior:** Highly precise over short-to-medium horizons due to its 4th-order global accuracy. However, because it is non-symplectic, it will still experience a slow, secular energy drift over long integration spans.

### 3. Symplectic Leapfrog (2nd-Order Verlet)
* **The Math:** Staggers position and velocity updates by alternating half-steps (kick-drift-kick).
* **Behavior:** Though globally only 2nd-order accurate, its *symplectic* architecture forces it to exactly conserve a "shadow" Hamiltonian near the true one. Instead of drifting endlessly, its energy error stays bounded within a tight, oscillating window indefinitely.

---

## Installation & Local Execution

This dashboard is engineered as a single, self-contained Python script to make deployment seamless.

### Prerequisites
Ensure you have Python 3 installed along with the required scientific computing libraries:
```bash
pip install streamlit numpy plotly
```

### Running the Application
Clone this repository, navigate to the directory, and run the following command to spin up the dashboard locally:
```bash
streamlit run numerical_integrator_sandbox.py
```
___

## Classroom Application Ideas

**The Energy Conservation Challenge:** Have students look at the Energy Drift Analysis tab at a very small timestep (e.g., dt=0.005) and observe how RK4 achieves an initially lower error than Leapfrog, but Leapfrog stays bounded while RK4 secularly drifts over time.

**The Chaos Trigger:** Instruct students to gradually drag the dt slider above 0.1. Watch the Forward Euler trajectory immediately break down into an unstable outward spiral while Leapfrog maintains an eccentric loop.
