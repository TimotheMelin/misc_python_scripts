from numpy import *
import os, re, sys
from IPython.display import display


# # CBS 3pts: E_n = E_CBS + Ae(-(n-1)) + Be(-(n-2)^2), n is zeta of basis set
def extrapolation_3pts(input_list,elem,method,verbose=None):
     e_vals = array(input_list.loc['{}.dz'.format(method):'{}.qz'.format(method),elem])

     E_CBS = ((1+exp(2))*e_vals[0] - (exp(1) + exp(3) + exp(5))*e_vals[1] + exp(6)*e_vals[2])/((exp(1)-1)*(exp(5) - exp(2) - 1))

     return E_CBS



# En = E_CBS + A/((n+0.5)**4)
def extrapolation_2pts_4(input_list,elem,calc,basis):
    if 'dz' in basis[0] and 'tz' in basis[1]:
        n1 = 2
        n2 = 3
    if 'tz' in basis[0] and 'qz' in basis[1]:
        n1 = 3
        n2 = 4

    E1 = input_list.loc['{}.{}'.format(calc,basis[0]),elem]
    E2 = input_list.loc['{}.{}'.format(calc,basis[1]),elem]

    E_CBS = ((E2/((0.5+n1)**4)) - (E1/(0.5+n2)**4))/((1/(0.5+n1)**4)-(1/(0.5+n2)**4))

    return E_CBS

# En = E_CBS + A/n**3
def extrapolation_2pts_3(input_list,elem,calc,basis):
    if 'dz' in basis[0] and 'tz' in basis[1]:
        n1 = 2
        n2 = 3
    if 'tz' in basis[0] and 'qz' in basis[1]:
        n1 = 3
        n2 = 4

    E1 = input_list.loc['{}.{}'.format(calc,basis[0]),elem]
    E2 = input_list.loc['{}.{}'.format(calc,basis[1]),elem]

#     A = ((E2 - E1)*n1**3*n2**3)/((n1**3-n2**3))
    E_CBS = ((E1*n1**3) - (E2*n2**3))/((n1**3)-n2**3)
    return E_CBS
