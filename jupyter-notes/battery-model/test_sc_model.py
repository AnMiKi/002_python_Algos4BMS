######################################################################
## This file is part of the Battery Model project.
## It is dedicted to test the class implemnted for the self correcting
## cell model implemendted in sc_model.py
######################################################################

from sc_model import SCModel
import numpy as np

# #############################
# # Test of SoC calculation
# #############################
# # init an instance
# model = SCModel(capacity_ah = 5.0, r0 = 10e-3, r1 = 10e-2, tau1 = 10)

# # Test Zelle entladen (positives Vorzeichen)
# x_next = model.f(x = [0.5, 0.0], current = 2.0, dt = 10)

# print(x_next[1])

# assert abs(x_next[1] - 1.264241) < 1e-4

# # Test der Relaxation
# x_next = model.f(x = [0.8 , 2.0], current = 0.0, dt = 10.0)

# print(x_next[1])

# assert abs(x_next[1] - 0.7357588823428847) < 1e-12

# #############################
# # Test of OCV calculation
# #############################

# # Test @ 100% SoC
# ocv = model.ocv(x = 1.0, current = 0.0)

# print(ocv)

# assert abs(ocv - 4.2) < 1e-12

# # Test @ 0% SoC
# ocv = model.ocv(x = 0.0, current = 0.0)

# print(ocv)

# assert abs(ocv - 3.0) < 1e-12

# #############################
# # Test of the output equation
# #############################

# # No current
# assert abs(model.h([0.5, 0.0], 0.0) - 3.6) < 1e-4

# # Discharge: terminal voltage must drop
# assert abs(model.h([0.5,0.0], 10.0) - 3.5) < 1e-4

# # Charge: terminal voltage must rise
# assert abs(model.h([0.5,0.0], -10.0) - 3.7) < 1e-4

#############################
# Test of the hysteresis (dynamic part)
#############################

# model = SCModel(capacity_ah = 5, r0 = 10e-3, r1 = 10e-2, tau1 = 10, gamma = 900, m = 30e-3);

# x = np.array([
#     0.8,
#     0.0,
#     0.0
# ])

# # UT-1 -> I = 2A => Entladen
# x_next = model.f(x, 2, 10);

# print(x_next[2])

# # UT-2 -> I = -2A => Laden
# x_next = model.f(x, -2, 10);

# print(x_next[2])

# # UT-3 -> I = 0A, dt = 1000s

# x = np.array([
#     0.8,
#     0.0,
#     0.7
# ])

# x_next = model.f(x, 0, 1000);

# print(x_next[2])


#############################
# Test of the hysteresis (instant part)
#############################
# model = SCModel(capacity_ah = 5, r0 = 10e-3, r1 = 10e-2, tau1 = 10, gamma = 900, m = 30e-3, m0 = 1);

# # Test des Vorzeichens in Entladerichtung
# s = model.update_sign_state(current = 2.0, s_prev = -1)
# print(s)
# assert(np.abs(s - 1)) < 1e-12

# # Test des Vorzeichens in Laderichtung
# s = model.update_sign_state(current = -2.0, s_prev = -1)
# print(s)
# assert(np.abs(s + 1)) < 1e-12

# # Test von sher kleinen Strömen
# # Test des Vorzeichens in Laderichtung
# s = model.update_sign_state(current = 0.01, s_prev = -1)
# print(s)
# assert(np.abs(s + 1)) < 1e-12 # soll auf virherigem Wert bleiben

#############################
# Test of the effective current calculation
#############################
model = SCModel(capacity_ah = 1, r0 = 10e-3, r1 = 10e-2, tau1 = 10, gamma = 900, m = 30e-3, m0 = 1, eta = 0.9);

# discharge
x = np.array([0.5, 0.0, 0.0])
x_next = model.f(x, current=0.1, dt=3600)

assert abs(x_next[0] - 0.4) < 1e-12

x = np.array([0.5, 0.0, 0.0])
x_next = model.f(x, current=-0.1, dt=3600)

assert abs(x_next[0] - 0.59) < 1e-12

