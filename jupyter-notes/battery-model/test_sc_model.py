######################################################################
## This file is part of the Battery Model project.
## It is dedicted to test the class implemnted for the self correcting
## cell model implemendted in sc_model.py
######################################################################

from sc_model import SCModel


#############################
# Test of SoC calculation
#############################
# init an instance
model = SCModel(capacity_ah = 5.0, r0 = 10e-3)

# Test Zelle entladen (positives Vorzeichen)
z_next = model.f(x = 0.8, current = 1.0, dt = 3600.0)

print(z_next)

assert abs(z_next - 0.6) < 1e-12

# Test Zelle laden (negative Vorzeichen)
z_next = model.f(x = 0.8, current = -1.0, dt = 3600.0)

print(z_next)

assert abs(z_next - 1.0) < 1e-12

#############################
# Test of OCV calculation
#############################

# Test @ 100% SoC
ocv = model.ocv(x = 1.0)

print(ocv)

assert abs(ocv - 4.2) < 1e-12

# Test @ 0% SoC
ocv = model.ocv(x = 0.0)

print(ocv)

assert abs(ocv - 3.0) < 1e-12

#############################
# Test of the output equation
#############################

# Test beim Entladen der Zelle
v = model.h(x = 0.5833, current = 10)

print(v)

assert abs(v - 3.6) < 1e-4

# Test beim Laden der Zelle
v = model.h(x = 0.5833, current = -10)

print(v)

assert abs(v - 3.8) < 1e-4