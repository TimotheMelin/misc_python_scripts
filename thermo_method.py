from numpy import *
import os, re, sys
from IPython.display import display
import math

import extrapolation as extrap
import misceleneous as mis

# f-ccCA version using UCCSD(T)
def f_ccCA(input_list,elem,SO=False,verbose=False):

    E_CBS_mp2 = extrap.extrapolation_3pts(input_list,elem,'corre')
    E_CBS_hf = extrap.extrapolation_2pts_hf(input_list,elem,'hf')
    E_cv = input_list.loc['mp2.pcvdz',elem] - input_list.loc['mp2.dz',elem]
    E_cc =  input_list.loc['ccsdt.tz',elem] - input_list.loc['mp2.tz',elem]
    E_ccCA = E_CBS_mp2 + E_CBS_hf + E_cv + E_cc

    if verbose:
        so = input_list.loc['SO',elem]
        print('E_CBS_MP2',E_CBS_mp2)
        print('E_CBS_HF',E_CBS_hf)
        print('E_CV',E_cv)
        print('E_CC',E_cc)
        print('E_ccCA',E_ccCA)
        print('E_SO',so)
        print('!!! E_SO must be in hartree !!!')
        print('E_ccCA + SO',E_ccCA+so)

    if SO == True:
        so = input_list.loc['SO',elem]
        return E_ccCA + so
    else:
        return E_ccCA

#Main ccCA function for organic molecules. SO is directly taken from experiment with the J average.
def ccCA(input_list,elem,verbose=False):

    so = input_list.loc['SO',elem]
    E_CBS_mp2 = extrap.extrapolation_3pts(input_list,elem,'corre')
    E_CBS_hf = extrap.extrapolation_2pts_hf(input_list,elem,'hf')
    E_cc = input_list.loc['ccsdt.dz',elem] - input_list.loc['mp2.dz',elem]
    E_cv = input_list.loc['mp2.pcvtz',elem] - input_list.loc['mp2.atz',elem]
    E_sr = input_list.loc['mp2.tz.dk',elem] - input_list.loc['mp2.tz',elem]

    E_ccCA = E_CBS_hf + E_CBS_mp2 + E_cc + E_cv + E_sr + so

    if verbose:
        print('E_CBS_MP2',E_CBS_mp2)
        print('E_CBS_HF',E_CBS_hf)
        print('E_CV',E_cv)
        print('E_CC',E_cc)
        print('E_SR',E_sr)
        print('E_SO',so)
        print('!!! E_SO must be in hartree and be calculated before using this function and added in the excelsheet !!!')
        print('E_ccCA + SO',E_ccCA)

    else:
        return E_ccCA


# DLPNO-ccCA for ogranic molecules
def DLPNO_ccCA(input_list,elem,verbose=False):

    so = input_list.loc['SO',elem]
    E_CBS_mp2 = extrap.extrapolation_3pts(input_list,elem,'DLPNO-corre')
    E_CBS_hf = extrap.extrapolation_2pts_hf(input_list,elem,'DLPNO-hf')
    E_cc = input_list.loc['DLPNO-ccsdt.dz',elem] - input_list.loc['DLPNO-mp2.dz',elem]
    E_cv = input_list.loc['DLPNO-mp2.pcvtz',elem] - input_list.loc['DLPNO-mp2.atz',elem]
    E_sr = input_list.loc['DLPNO-mp2.tz.dk',elem] - input_list.loc['DLPNO-mp2.tz',elem]

    E_ccCA = E_CBS_hf + E_CBS_mp2 + E_cc + E_cv + E_sr + so

    if verbose:
        print('E_CBS_MP2',E_CBS_mp2)
        print('E_CBS_HF',E_CBS_hf)
        print('E_CV',E_cv)
        print('E_CC',E_cc)
        print('E_SR',E_sr)
        print('E_SO',so)
        print('!!! E_SO must be in hartree and be calculated before using this function and added in the excelsheet !!!')
        print('E_ccCA + SO',E_ccCA)

    else:
        return E_ccCA

# rp-ccCA for 5d compounds. CBS only at Mp2 energies and So only for atoms
def rp_ccCA(input_list,elem,verbose=False):

    E_cc = input_list.loc['ccsdt.tz', elem] - input_list.loc['mp2.tz', elem]
    E_cv = input_list.loc['ccsdt.pcvdz', elem] - input_list.loc['ccsdt.dz', elem]
    E_CBS_mp2 = extrap.extrapolation_3pts(input_list,elem,'mp2')
    so = input_list.loc['SO',elem]

    E_ccCA = E_CBS_mp2 + E_cc + E_cv + so

    if verbose:
        print('E_CBS_MP2',E_CBS_mp2)
        print('E_CV',E_cv)
        print('E_CC',E_cc)
        print('E_SO',so)
        print('!!! E_SO must be in hartree and be calculated before using this function and added in the excelsheet !!!')
        print('E_ccCA + SO',E_ccCA)

    else:
        return E_ccCA

