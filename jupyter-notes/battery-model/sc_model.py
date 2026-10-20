######################################################################
## This file is part of the Battery Model project.
## It implements a class for the enhanced self-correcting cell model.
######################################################################

import numpy as np

class SCModel:
    
    '''Self Correcting Cell Model'''

    def __init__(self,
                capacity_ah,
                r0,
                r,
                tau,
                gamma,
                m,
                m0,
                eta_charge, 
                soc_grid,
                ocv0_values,
                ocv_rel_values):
        
        # Def of initial model params
        self.capacity_ah = capacity_ah
        self.r0 = r0
        self.r = np.asarray(r, dtype=float)
        self.tau = np.asarray(tau, dtype=float)
        if len(self.r) != len(self.tau):
            raise ValueError("Length of r and tau must be the same.")
        self.gamma = gamma
        self.m = m
        self.m0 = m0
        self.eta_charge = eta_charge
        self.soc_grid = np.asarray(soc_grid, dtype= float)
        if np.any(np.diff(self.soc_grid) <= 0):
            raise ValueError("soc_grid must be strictly increasing.")
        if len(self.soc_grid) < 2:
            raise ValueError("At least two OCV/SOC points are required.")
        self.ocv0_values = np.asarray(ocv0_values, dtype=float)
        if len(self.soc_grid) != len(self.ocv0_values):
            raise ValueError("Length of soc_grid and ocv0_values must be the same.")
        self.ocv_rel_values = np.asarray(ocv_rel_values, dtype=float)
        if len(self.soc_grid) != len(self.ocv_rel_values):
            raise ValueError("Length of soc_grid and ocv_rel_values must be the same.")
        
    # State equation
    def f(self, x, current, dt):

        z = x[0] # SoC
        i_r = x[1:-1] # RC currents 
        hyst = x[-1] # dynamic hysteresis

        i_eff = self.effective_current(current)

        # 1. SOC state: Coulomb counting
        z_next = z - (i_eff * dt) / (3600 * self.capacity_ah)

        # 2. RC polarization state
        a = np.exp(-dt/self.tau) 
        i_r_next = a * i_r + (1 - a) * i_eff 

        # 3. Dynamic hysteresis state
        a_h = np.exp((-np.abs(self.gamma * i_eff * dt)/(3600 * self.capacity_ah)))
        hyst_next = a_h * hyst - (1 - a_h) * np.sign(current)

        x_next = np.concatenate((
            [z_next],
            i_r_next,
            [hyst_next]
        ))

        return x_next

    # OCV calculation
    def ocv(self, soc, temperature):
        
        ocv0 = self._interp_extrap(soc, self.ocv0_values)
        ocv_rel = self._interp_extrap(soc, self.ocv_rel_values)

        ocv = ocv0 + temperature * ocv_rel

        return ocv

    # Output equation
    def h(self, x, current, s, temperature):

        z = x[0]
        i_r = x[1:-1]
        hyst = x[-1]

        v = self.ocv(z, temperature)  
        v = v + (self.m * hyst) + (self.m0 * s) - (self.r0 * current) - np.dot(self.r ,i_r)

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

    # Calculation of the slope of the OCV curve at a given SoC
    def docv_dsoc(self, soc, temperature):
        slope_ocv0 = self._slope(soc, self.ocv0_values)
        slope_ocv_rel = self._slope(soc, self.ocv_rel_values)
        slope = slope_ocv0 + temperature * slope_ocv_rel

        return slope

    def _interp_extrap(self, soc, values):

        if soc < self.soc_grid[0]:
        # lower extrapolation
            slope_low = (values[1] - values[0])/(self.soc_grid[1] - self.soc_grid[0])
            ocv = values[0] + slope_low * (soc - self.soc_grid[0])  
        elif soc > self.soc_grid[-1]:
        # upper extrapolation
            slope_high = (values[-1] - values[-2])/(self.soc_grid[-1] - self.soc_grid[-2])
            ocv = values[-1] + slope_high * (soc - self.soc_grid[-1]) 
        else:
        # interpolation 
            ocv = np.interp(soc, self.soc_grid, values)
        return ocv

    def _slope(self, soc, values):
        if soc <= self.soc_grid[0]:
            # lower extrapolation
            slope_low = (values[1] - values[0])/(self.soc_grid[1] - self.soc_grid[0])
            return slope_low
        elif soc >= self.soc_grid[-1]:
            # upper extrapolation
            slope_high = (values[-1] - values[-2])/(self.soc_grid[-1] - self.soc_grid[-2])
            return slope_high
        else:
            # interpolation
            idx = np.searchsorted(self.soc_grid, soc) - 1
            slope = (values[idx + 1] - values[idx]) / (self.soc_grid[idx + 1] - self.soc_grid[idx])
            return slope