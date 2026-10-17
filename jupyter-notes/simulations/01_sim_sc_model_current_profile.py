from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "battery-model"))

from sc_model import SCModel
import numpy as np
import matplotlib.pyplot as plt

model = SCModel(capacity_ah=5.0, r0=0.01, r1=0.02, tau1=10, gamma=0.1, m=0.05, m0=0.02)

# Definition der Simulationszeiten
dt = 1.0
t_end = 180.0

time = np.arange(0, t_end + dt, dt)
#print(time)
current = np.zeros(len(time))
current[(time >= 20) & (time < 60)] = 2.0
current[(time >= 100) & (time < 140)] = -2.0
print(current)

# Initial states
x = np.array([
    0.8,   # SOC
    0.0,   # RC current
    0.0    # dynamic hysteresis
])

s = 0.0

# Storage for results
soc_store = np.zeros(len(time))
ir_store = np.zeros(len(time))
hyst_store = np.zeros(len(time))
s_store = np.zeros(len(time))
voltage_store = np.zeros(len(time))

# Store initial states
soc_store[0] = x[0]
ir_store[0] = x[1]
hyst_store[0] = x[2]
s_store[0] = s

################################
# Start of Simulation
################################
for k in range(len(time)-1):

    i_k = current[k]

    # Instantaneous hysteresis state
    s = model.update_sign_state(i_k, s)

    # Output from current state
    voltage_store[k] = model.h(x, i_k, s)

    # Continuous-state propagation
    x = model.f(x, i_k, dt)

    # Store states
    soc_store[k + 1] = x[0]
    ir_store[k + 1] = x[1]
    hyst_store[k + 1] = x[2]
    s_store[k + 1] = s

s = model.update_sign_state(current[-1], s)
voltage_store[-1] = model.h(x, current[-1], s)
s_store[-1] = s

fig, axes = plt.subplots(3, 1, sharex=True, figsize=(10, 8))

axes[0].plot(time, ir_store)
axes[0].set_ylabel("RC current [A]")
axes[0].grid()

axes[1].plot(time, hyst_store)
axes[1].set_ylabel("Dynamic hysteresis h [-]")
axes[1].grid()

axes[2].plot(time, s_store)
axes[2].set_ylabel("Instantaneous hysteresis s [-]")
axes[2].set_xlabel("Time [s]")
axes[2].grid()

fig.tight_layout()
plt.show()

fig, axes = plt.subplots(6, 1, sharex=True, figsize=(10, 8))

axes[0].plot(time, current)
axes[0].set_ylabel("Input current [A]")
axes[0].grid()

axes[1].plot(time, voltage_store)
axes[1].set_ylabel("Terminal voltage [V]")
axes[1].grid()

axes[2].plot(time, 0.05 * hyst_store)
axes[2].set_ylabel("Dynamic hysteresis * m [-]")
axes[2].set_xlabel("Time [s]")
axes[2].grid()

axes[3].plot(time, 0.02 * s_store)
axes[3].set_ylabel("Instantaneous hysteresis m * s [-]")
axes[3].set_xlabel("Time [s]")
axes[3].grid()

axes[4].plot(time,-1 *  ir_store * 0.02)
axes[4].set_ylabel("RC voltage  [V]")
axes[4].set_xlabel("Time [s]")
axes[4].grid()

axes[5].plot(time, current * 0.01)
axes[5].set_ylabel("Instantaneous voltage drop [V]")
axes[5].set_xlabel("Time [s]")
axes[5].grid()

fig.tight_layout()
plt.show()