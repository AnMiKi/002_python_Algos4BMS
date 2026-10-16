######################################################################
## This file is for testing, probing, plotting, playing around
######################################################################

from sc_model import SCModel
import numpy as np
import matplotlib.pyplot as plt

model = SCModel(capacity_ah=5.0, r0=0.01, r1=0.02, tau1=10)

# Definition der Simulationszeiten
dt = 1.0
t_end = 100.0

# Anfangszustand
x = np.array([
    0.8, #SoC
    0.0 # RC current
])

# Erzeugen der Zeitachse
time = np.arange(0, t_end + dt, dt)

# Erzeugen des Stromprofils
current = np.zeros(len(time))
pulse = (time >= 20) & (time < 60)
current[pulse] = 2.0

# Speicher allokieren
soc_store = np.zeros(len(time))
ir_store = np.zeros(len(time))
voltage_store = np.zeros(len(time))
ocv_store = np.zeros(len(time))
r0_store = np.zeros(len(time))
rc_store = np.zeros(len(time))

# Anfangszustand speichern
soc_store[0] = x[0]
ir_store[0] = x[1]

##################################
# Start of Simulation
##################################
for k in range(len(time) - 1):

    # Current at this time step
    i_k = current[k]

    # 1. Calculate terminal voltage from current state
    voltage_store[k] = model.h(x, i_k)
    ocv_store[k] = model.ocv(x[0])
    r0_store[k] = model.r0 * i_k
    rc_store[k] = model.r1 * x[1]

    # 2. Propagate state to next sample
    x = model.f(x, i_k, dt)

    # 3. Store new state
    soc_store[k + 1] = x[0]
    ir_store[k + 1] = x[1]

voltage_store[-1] = model.h(x, current[-1])
ocv_store[-1] = model.ocv(x[0])
r0_store[-1] = model.r0 * current[-1]
rc_store[-1] = model.r1 * x[1]

plt.figure(1)

plt.plot(time, current)

plt.xlabel("Time [s]")
plt.ylabel("Current [A]")
plt.grid()

plt.show()

plt.figure(2)

plt.plot(time, ir_store)

plt.xlabel("Time [s]")
plt.ylabel("RC current [A]")
plt.grid()

plt.show()

plt.figure(3)

plt.plot(time, voltage_store)

plt.xlabel("Time [s]")
plt.ylabel("Terminal voltage [V]")
plt.grid()

plt.show()

plt.figure()

plt.plot(time, voltage_store, label="Terminal voltage")
plt.plot(time, ocv_store, label="OCV")

plt.xlabel("Time [s]")
plt.ylabel("Voltage [V]")
plt.legend()
plt.grid()

plt.show()