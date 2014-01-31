#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Gravitational constant
G = 6.67e-11
PC_in_m = 3.1e16
Msun = 1.99e30

__C_AU              = 1.495978707e+11        # m
__C_GRAVIT          = 6.6742e-11             # m³/kg/s²
__C_SOLMASS         = 1.9891e+30             # kg
__C_SOLMASS_AU3DAY2 = 2.9591220828559115e-04 # au³/day²


# mass constants - base unit: au³/day²
CONV_M_SOLMASS = (1./__C_SOLMASS_AU3DAY2)

# length constants - base unit au
CONV_X_PARSEC = 4.848136813e-6
CONV_X_KM     = 149597870.7

# time constants - base unit day
CONV_T_YEAR   = 0.2737909330e-2
CONV_T_SECOND = 86400.
CONVINT_X_AU_SQRT      = 32.  # sqrt(1000)
CONVINT_M_AU3DAY2_SQRT = 100.

CONVINT_X_AU      = (CONVINT_X_AU_SQRT * CONVINT_X_AU_SQRT)
CONVINT_M_AU3DAY2 = (CONVINT_M_AU3DAY2_SQRT * CONVINT_M_AU3DAY2_SQRT)
CONVINT_T_DAY     = (CONVINT_X_AU_SQRT * CONVINT_X_AU_SQRT * CONVINT_X_AU_SQRT / CONVINT_M_AU3DAY2_SQRT)

CONV_M = CONV_M_SOLMASS * CONVINT_M_AU3DAY2
CONV_X = CONV_X_PARSEC  * CONVINT_X_AU
CONV_T = CONV_T_YEAR    * CONVINT_T_DAY
CONV_V = CONV_X / CONV_T
CONV_E = (CONV_M * CONV_M / CONV_X)
