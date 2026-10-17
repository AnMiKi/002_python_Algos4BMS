######################################################################
## This file is part of the Battery Model project.
## It implements a class for the enhanced self-correcting cell model.
######################################################################

import numpy as np

class SCModel:
    
    '''Self Correcting Cell Model'''

    def __init__(self, capacity_ah, r0, r1, tau1, gamma, m, m0, eta_charge):
        # Def of initial model params
        self.capacity_ah = capacity_ah
        self.r0 = r0
        self.r1 = r1
        self.tau1 = tau1
        self.gamma = gamma
        self.m = m
        self.m0 = m0
        self.eta_charge = eta_charge

    # TBD equation
    def f(self, x, current, dt):

        z = x[0]
        i_r = x[1]
        hyst = x[2]

        i_eff = self.effective_current(current)

        # 1. SOC state: Coulomb counting
        z_next = z - (i_eff * dt) / (3600 * self.capacity_ah)

        # 2. RC polarization state
        a = np.exp(-dt/self.tau1) 
        i_r_next = a * i_r + (1 - a) * i_eff 

        # 3. Dynamic hysteresis state
        a_h = np.exp((-np.abs(self.gamma * i_eff * dt)/(3600 * self.capacity_ah)))
        hyst_next = a_h * hyst - (1 - a_h) * np.sign(current)

        x_next = np.array([
            z_next,
            i_r_next,
            hyst_next
        ])

        return x_next

    # Simple OCV model
    def ocv(self, soc):
        ocv = 3.0 + 1.2 * soc
        return ocv

    # Output equation
    def h(self, x, current, s):

        z = x[0]
        i_r = x[1]
        hyst = x[2]

        v = self.ocv(z)  
        v = v + (self.m * hyst) + (self.m0 * s) - (self.r0 * current) - (self.r1 * i_r)

        return v

    # instante Hysterese
    def update_sign_state(self, current, s_prev):

        s = np.sign(current)

        if np.abs(current) < self.capacity_ah / 100:
            s = s_prev

        return s

    # Effective current for charge/discharge efficiency
    def effective_current(self, current):
        if current >= 0:
            return current
        else:
            return current * self.eta_charge