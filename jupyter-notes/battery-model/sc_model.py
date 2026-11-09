######################################################################
## This file is part of the Battery Model project.
## It implements a class for the enhanced self-correcting cell model.
######################################################################

import numpy as np

class SCModel:
    """Self Correcting Cell Model."""

    def __init__(
        self,
        capacity_values,
        r0_values,
        r_values,
        tau_values,
        gamma,
        m,
        m0,
        eta_charge,
        soc_grid,
        ocv0_values,
        ocv_rel_values,
        temperature_grid,
    ):

        """Initialize the Enhanced Self-Correcting cell model."""

        # ------------------------------------------------------------------
        # Temperature grid
        # ------------------------------------------------------------------
        self.temperature_grid = np.asarray(temperature_grid, dtype=float)

        if self.temperature_grid.ndim != 1:
            raise ValueError("temperature_grid must be a 1D array.")

        if np.any(np.diff(self.temperature_grid) <= 0):
            raise ValueError("temperature_grid must be strictly increasing.")

        n_temperatures = len(self.temperature_grid)

        # ------------------------------------------------------------------
        # Temperature-dependent scalar parameters
        # ------------------------------------------------------------------
        self.capacity_values = np.asarray(capacity_values, dtype=float)
        self.r0_values = np.asarray(r0_values, dtype=float)

        if self.capacity_values.ndim != 1:
            raise ValueError("capacity_values must be a 1D array.")

        if self.r0_values.ndim != 1:
            raise ValueError("r0_values must be a 1D array.")

        if len(self.capacity_values) != n_temperatures:
            raise ValueError(
                "Length of capacity_values must match temperature_grid."
            )

        if len(self.r0_values) != n_temperatures:
            raise ValueError(
                "Length of r0_values must match temperature_grid."
            )

        # ------------------------------------------------------------------
        # Temperature-dependent RC parameters
        # Shape: (number of temperatures, number of RC branches)
        # ------------------------------------------------------------------
        self.r_values = np.asarray(r_values, dtype=float)
        self.tau_values = np.asarray(tau_values, dtype=float)

        if self.r_values.ndim != 2:
            raise ValueError(
                "r_values must be a 2D array with one row per temperature."
            )

        if self.tau_values.ndim != 2:
            raise ValueError(
                "tau_values must be a 2D array with one row per temperature."
            )

        if self.r_values.shape != self.tau_values.shape:
            raise ValueError(
                "r_values and tau_values must have the same shape."
            )

        if self.r_values.shape[0] != n_temperatures:
            raise ValueError(
                "Number of temperature rows in r_values and tau_values "
                "must match temperature_grid."
            )

        # ------------------------------------------------------------------
        # ESC parameters
        # ------------------------------------------------------------------
        self.gamma = gamma
        self.m = m
        self.m0 = m0
        self.eta_charge = eta_charge

        # ------------------------------------------------------------------
        # OCV-SOC lookup tables
        # ------------------------------------------------------------------
        self.soc_grid = np.asarray(soc_grid, dtype=float)
        self.ocv0_values = np.asarray(ocv0_values, dtype=float)
        self.ocv_rel_values = np.asarray(ocv_rel_values, dtype=float)

        if self.soc_grid.ndim != 1:
            raise ValueError("soc_grid must be a 1D array.")

        if len(self.soc_grid) < 2:
            raise ValueError(
                "At least two OCV/SOC support points are required."
            )

        if np.any(np.diff(self.soc_grid) <= 0):
            raise ValueError("soc_grid must be strictly increasing.")

        if self.ocv0_values.ndim != 1:
            raise ValueError("ocv0_values must be a 1D array.")

        if self.ocv_rel_values.ndim != 1:
            raise ValueError("ocv_rel_values must be a 1D array.")

        if len(self.ocv0_values) != len(self.soc_grid):
            raise ValueError(
                "Length of ocv0_values must match soc_grid."
            )

        if len(self.ocv_rel_values) != len(self.soc_grid):
            raise ValueError(
                "Length of ocv_rel_values must match soc_grid."
            )
        
    # State equation
    def f(self, x, current, dt, temperature):

        capacity_ah = self._interp_temperature(temperature, self.capacity_values) # Interpolate capacity based on temperature

        x = self._validate_state(x)  # Validate the state vector

        # Temperaturabhängige Parameter bestimmen
        tau = self._interp_temperature_vector(temperature, self.tau_values)

        z = x[0] # SoC
        i_r = x[1:-1] # RC currents 
        hyst = x[-1] # dynamic hysteresis

        i_eff = self.effective_current(current)

        # 1. SOC state: Coulomb counting
        z_next = z - (i_eff * dt) / (3600 * capacity_ah)

        # 2. RC polarization state
        a = np.exp(-dt/tau) 
        i_r_next = a * i_r + (1 - a) * i_eff 

        # 3. Dynamic hysteresis state
        a_h = np.exp((-np.abs(self.gamma * i_eff * dt)/(3600 * capacity_ah)))
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

        x = self._validate_state(x)  # Validate the state vector

        z = x[0]
        i_r = x[1:-1]
        hyst = x[-1]

        # Temperaturabhängiger Innenwiderstand
        r0  = self._interp_temperature(temperature, self.r0_values)

        # Temperaturabhängiger RC-Widerstand
        r = self._interp_temperature_vector(temperature, self.r_values)

        v = self.ocv(z, temperature)  
        v = v + (self.m * hyst) + (self.m0 * s) - (r0 * current) - np.dot(r ,i_r)

        return v

    # instante Hysterese
    def update_sign_state(self, current, s_prev, temperature):

        capacity_ah = self._interp_temperature(temperature, self.capacity_values) # Interpolate capacity based on temperature

        s = np.sign(current)

        if np.abs(current) < capacity_ah / 100:
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

    def _interp_temperature(self, temperature, values):
        # Temperatur begrenzen auf min oder max der Temperaturwerte
        temperature_clamped =np.clip(temperature,
                                     self.temperature_grid[0],
                                     self.temperature_grid[-1])
        value = np.interp(temperature_clamped,
                          self.temperature_grid,
                          values)

        return value

    # Interpolation of a vector of values based on temperature
    def _interp_temperature_vector(self, temperature, values):
        temperature_clamped =np.clip(temperature,
                                     self.temperature_grid[0],
                                     self.temperature_grid[-1])

        n_rc = values.shape[1] # Number of RC elements
        result = np.zeros(n_rc) # Initialize result array
        for j in range(n_rc):
            result[j] = np.interp(
            temperature_clamped,
            self.temperature_grid,
            values[:, j])
        return result # Return the interpolated values

    def _validate_state(self, x):

        x = np.asarray(x, dtype=float) # Ensure x is a numpy array of floats

        n_rc = self.r_values.shape[1] # Number of RC elements
        expected_length = 1 + n_rc + 1 # SOC + RC currents + dynamic hysteresis

        if x.ndim != 1:
            raise ValueError(f"State vector x must be a 1D array of length {expected_length}.")
        if len(x) != expected_length:
            raise ValueError(f"State vector x must have length {expected_length}, but got length {len(x)}.")

        return x