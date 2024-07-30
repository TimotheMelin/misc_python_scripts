from numpy import *
import os, re, sys
from IPython.display import display
import math
from scipy.constants import physical_constants
from scipy.constants import calorie, N_A

# Define strandard units
invcm2kcal = (physical_constants['inverse meter-joule relationship'][0]*N_A/10)/calorie
hartree2kcal = physical_constants['hartree-joule relationship'][0]*N_A/1000/calorie
hartree2kj = 2625.5002


# Create 2 arrays with list of atoms and number of appearance
def create_atom_cmpd_list(usr_inp):
    int_set = []
    elem_list = []
    a = re.findall(r'[A-Z][a-z]*|\d+|\(|\)', usr_inp)
    for p in range(1,len(a)):
        if a[p-1].isalpha() and a[p].isalpha():
            if p == len(a)-1:
                int_set.append(1.0)
                int_set.append(1.0)
                elem_list.append(str(a[p-1]))
                elem_list.append(str(a[p]))
            else:
                int_set.append(1.0)
                elem_list.append(str(a[p-1]))

        elif a[p-1].isalpha() and a[p].isdigit():
            int_set.append(float(a[p]))
            elem_list.append(str(a[p-1]))
    if a[-2].isdigit() and a[-1].isalpha():
        int_set.append(1)
        elem_list.append(str(a[-1]))

    return array(elem_list), array(int_set)

#def ordered_elemnames(exp_set):
#    M_list = []
#    for i in exp_set.axes[0]:
#
#        if re.search('[2-9]',i):
#            e = re.search('[2-9]',i).group()
#            e = re.sub('\d','$_{}$'.format(e),i)
#            M_list.append(e)
#        else:
#            M_list.append(i)
#    return  M_list

# Create list of molecules name for plotting
def ordered_elemnames(molecules):
    list_M = []
    for j in range(len(molecules)):
        atoms_list, number_atoms = create_atom_cmpd_list(molecules[j])
        name = ''
        for i in range(len(atoms_list)):
            a = '{}'.format(atoms_list[i])
            if number_atoms[i] == 1.0:
                n = ''
            else:
                n = '$_{{{}}}$'.format(int(number_atoms[i]))
            name = name + a + n
        list_M.append(name)
    return list_M

# RMSD function
def RMSD(target,values):
    total = 0
    for i in range(len(values)):
        total += sqrt((values[i]-target)**2)
    return mean(total)

# J average to calculate atomic SO. Data from NIST
def SO_atoms(Js_array,Energy_Js):
    return sum((2*Js_array + 1)*Energy_Js)/sum(2*Js_array + 1)*invcm2kcal


# Calculate Molecule(SO) based on experimental data
def calc_so_exp(elem):
    atom_set, int_list = create_atom_cmpd_list(elem)
    input_list = atoms_input
    SO = 0
    for i in range(len(atom_set)):
        SO += atoms_input.loc['SO_kcal', atom_set[i]]*int_list[i]

    return SO/hartree2kcal


