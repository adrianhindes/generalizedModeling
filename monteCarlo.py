# -*- coding: utf-8 -*-
"""
Created on Mon May 20 16:35:53 2019

@author: hindesa
Monte Carlo Primer
"""

import numpy as np
import random
import math
from matplotlib import pyplot as plt

def get_rand(a, b):
    r = np.abs(a-b)
    c = random.uniform(0,1)
    return min([a,b]) +r*c