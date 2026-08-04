### Basic packages import

import numpy as np

### src imports

import Cosmology as cosmo
import Constants_Units as cu

### defining dN_dt

def dN_dt_func(SFRD, z_arr, n_gamma_val, f_esc_val):
    return SFRD * (1 + z_arr) ** 3 * n_gamma_val * f_esc_val / (cu.m_p * cu.kgToSM)

### Parameters in the reionization equation

X_H = 0.75 #Fraction of hydrogen
alpha_B = 2.59 * 10 ** (-13) #Case B recombination coefficient at T=3*10^4 K in cm^3 * s^-1

alpha_B_unit = alpha_B * cu.cmToMpc ** 3 / cu.sToYr #Coefficient in Mpc^3 * yr^-1

def n_H(z): #Mean proper number denisty of hydrogen atoms in Mpc^-3
    return X_H * cosmo.Om_b * cosmo.rho_c * (1 + z) ** 3 / (cu.m_p * cu.kgToSM)

def dt_dz(z): #In yr
    return -1 / ((1 + z) * cosmo.H(z))

def Test_Clump(z):
    return 3