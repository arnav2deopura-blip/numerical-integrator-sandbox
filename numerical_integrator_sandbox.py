import streamlit as st
import numpy as np
import plotly.graph_objects as go

st.set_page_config(page_title="The Numerical Integrator Sandbox", page_icon="🛰️", layout="wide", initial_sidebar_state="expanded")

GM = 1.0
EPS = 1e-12

def compute_derivatives(state):
    x, y, vx, vy = state
    r = np.sqrt(x * x + y * y)
    r3 = max(r * r * r, EPS)
    ax = -GM * x / r3
    ay = -GM * y / r3
    return np.array([vx, vy, ax, ay])

def compute_energy(state):
    x, y, vx, vy = state
    r = max(np.sqrt(x * x + y * y), EPS)
    v2 = vx * vx + vy * vy
    return 0.5 * v2 - GM / r

def euler_step(state, dt):
    return state + dt * compute_derivatives(state)

def rk4_step(state, dt):
    k1 = compute_derivatives(state)
    k2 = compute_derivatives(state + 0.5 * dt * k1)
    k3 = compute_derivatives(state + 0.5 * dt * k2)
    k4 = compute_derivatives(state + dt * k3)
    return state + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)

def leapfrog_step(state, dt):
    x, y, vx, vy = state
    r = np.sqrt(x * x + y * y)
    r3 = max(r * r * r, EPS)
    ax = -GM * x / r3
    ay = -GM * y / r3
    vx_half = vx + 0.5 * dt * ax
    vy_half = vy + 0.5 * dt * ay
    x_new = x + dt * vx_half
    y_new = y + dt * vy_half
    r_new = np.sqrt(x_new * x_new + y_new * y_new)
    r3_new = max(r_new * r_new * r_new, EPS)
    ax_new = -GM * x_new / r3_new
    ay_new = -GM * y_new / r3_new
    vx_new = vx_half + 0.5 * dt * ax_new
    vy_new = vy_half + 0.5 * dt * ay_new
    return np.array([x_new, y_new, vx_new, vy_new])

STEP_FUNCTIONS = {"Forward Euler": euler_step, "Classical RK4": rk4_step, "Leapfrog": leapfrog_step}

@st.cache_data(show_spinner=False)
def run_simulation(method, x0, vy0, t_total, dt):
    n_steps = max(int(t_total / dt), 1)
    state = np.array([x0, 0.0, 0.0, vy0])
    trajectory = np.zeros((n_steps + 1, 4))
    energies = np.zeros(n_steps + 1)
    trajectory[0] = state
    energies[0] = compute_energy(state)
    step_fn = STEP_FUNCTIONS[method]
    for i in range(1, n_steps + 1):
        state = step_fn(state, dt)
        trajectory[i] = state
        energies[i] = compute_energy(state)
    times = np.linspace(0.0, n_steps * dt, n_steps + 1)
    return times, trajectory, energies

COLORS = {"Forward Euler": "#EF553B", "Classical RK4": "#636EFA", "Leapfrog": "#00CC96"}
METHODS = ["Forward Euler", "Classical RK4", "Leapfrog"]

st.title("The Numerical Integrator Sandbox")
st.markdown(
    "An interactive comparison of three numerical integration schemes applied to a classic two-body "
    "gravitational orbit (central mass at the origin, **GM = 1**). Adjust the initial conditions and, most "
    "importantly, the **timestep size**, to see how each method's mathematical structure shapes its accuracy "
    "and long-term stability."
)

with st.sidebar:
    st.header("Initial Conditions")
    x0 = st.slider("Initial Position X", min_value=0.5, max_value=2.0, value=1.0, step=0.05)
    vy0 = st.slider("Initial Velocity Y", min_value=0.3, max_value=1.5, value=1.0, step=0.05)
    t_total = st.slider("Total Simulation Time", min_value=5.0, max_value=100.0, value=20.0, step=1.0)
    st.markdown("---")
    st.subheader("Timestep Size (dt)")
    dt = st.slider("dt (integration step)", min_value=0.001, max_value=0.5, value=0.01, step=0.005)
    st.caption(f"This will require **{int(t_total / dt):,}** integration steps per method.")
    st.markdown("---")
    st.caption("Central mass fixed at the origin with GM = 1. Initial velocity is purely tangential (vx₀ = 0).")

if dt > 0.1:
    st.warning(
        f"⚠️ **Numerical Instability Warning (dt = {dt:.3f}):** Forward Euler is a first-order method, so its "
        "local error scales with dt² and its global error scales with dt. At this timestep, Euler is injecting "
        "artificial energy into the system at every step, causing its orbit to spiral outward instead of closing "
        "on itself. Watch the Trajectory plot for the outward spiral and the Energy Drift plot for its runaway "
        "error curve. Try lowering dt below 0.1 to see Euler recover.",
    )

results = {}
with st.spinner("Integrating orbits with all three methods..."):
    for method in METHODS:
        results[method] = run_simulation(method, x0, vy0, t_total, dt)

tab1, tab2 = st.tabs(["Trajectory Comparison", "Energy Drift Analysis"])

with tab1:
    st.subheader("Orbital Trajectories")
    st.markdown(
        "All three integrators start from the exact same initial state and are advanced using the exact same "
        "timestep, so any divergence between the paths below is caused entirely by the numerical method itself, "
        "not by the physics. A perfectly integrated orbit under an inverse-square force is a closed ellipse "
        "(or, for high enough velocity, an open hyperbola/parabola)."
    )
    fig_traj = go.Figure()
    fig_traj.add_trace(
        go.Scatter(
            x=[0], y=[0], mode="markers", marker=dict(size=20, color="#FFD700", symbol="star",
            line=dict(color="white", width=1)), name="Central Mass",
        )
    )
    for method in METHODS:
        _, traj, _ = results[method]
        fig_traj.add_trace(
            go.Scatter(
                x=traj[:, 0], y=traj[:, 1], mode="lines", name=method,
                line=dict(color=COLORS[method], width=2.2),
            )
        )
    fig_traj.update_layout(
        template="plotly_dark", height=650, margin=dict(l=10, r=10, t=40, b=10),
        xaxis_title="x", yaxis_title="y",
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
    )
    fig_traj.update_yaxes(scaleanchor="x", scaleratio=1)
    st.plotly_chart(fig_traj, use_container_width=True)

with tab2:
    st.subheader("Relative Energy Error Over Time")
    st.markdown(
        "Specific orbital energy (kinetic + potential) is a conserved quantity in the true two-body problem. "
        "The plot below tracks `abs((E(t) - E(0)) / E(0))` on a **logarithmic scale**, which reveals three "
        "very different error signatures:\n\n"
        "- **Forward Euler** (1st order, non-symplectic): error grows roughly linearly with time and blows up "
        "quickly as dt increases, since the method systematically injects energy into the orbit.\n"
        "- **Classical RK4** (4th order, non-symplectic): highly accurate locally, so its error starts far lower "
        "than Euler's, but because it is still not symplectic, its error can drift slowly and secularly over "
        "very long integrations.\n"
        "- **Leapfrog** (2nd order, symplectic): its error does not grow over time. Instead it oscillates within "
        "a small, bounded band, because the method exactly conserves a nearby 'shadow' Hamiltonian rather than "
        "the true energy."
    )
    fig_energy = go.Figure()
    for method in METHODS:
        times, _, energies = results[method]
        e0 = energies[0]
        rel_err = np.abs((energies - e0) / e0) + 1e-16
        fig_energy.add_trace(
            go.Scatter(x=times, y=rel_err, mode="lines", name=method, line=dict(color=COLORS[method], width=2.2))
        )
    fig_energy.update_layout(
        template="plotly_dark", height=550, margin=dict(l=10, r=10, t=40, b=10),
        xaxis_title="Time", yaxis_title="Relative Energy Error |ΔE / E₀|",
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
    )
    fig_energy.update_yaxes(type="log")
    st.plotly_chart(fig_energy, use_container_width=True)

st.markdown("---")
with st.expander("About the three integration methods"):
    st.markdown(
        "**Forward Euler** advances the state using only the derivative evaluated at the current point: "
        "`state_next = state + dt * f(state)`. It is simple and cheap, but its first-order local truncation "
        "error and lack of symplectic structure make it a poor choice for long-term orbital integration.\n\n"
        "**Classical RK4** evaluates the derivative at four stages within each step (`k1` through `k4`) and "
        "combines them with fixed weights to cancel out lower-order error terms, achieving 4th-order global "
        "accuracy. It is far more accurate per step than Euler but is still not symplectic.\n\n"
        "**Leapfrog** (kick-drift-kick / velocity Verlet) staggers velocity and position updates by a half "
        "step. Despite being only 2nd-order accurate, its symplectic structure means it exactly conserves a "
        "Hamiltonian that stays close to the true one, so orbital energy oscillates but never systematically "
        "drifts, even over extremely long simulations."
    )