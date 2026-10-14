######################################################################
## This file is part of the Battery Model project.
## It implements a class for the advanced self correcting cell model.
######################################################################

import numpy as np

class SCModel:
    
    '''Self Correcting Cell Model'''

    def __init__(self, capacity_ah, r0):
        # Def of initial model params
        self.capacity_ah = capacity_ah
        self.r0 = r0

    # TBD equation
    def f(self, x, current, dt):
        x_next = x - (current * dt) / (3600 * self.capacity_ah)
        return x_next

    # Simple OCV model
    def ocv(self, x):
        ocv = 3.0 + 1.2 * x
        return ocv

    # Output equation
    def h(self, x, current):
        v = self.ocv(x) # - (self.r0 * current)
        return v
