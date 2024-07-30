import os, re, string
import sys, argparse
from numpy import *
import pandas as pd
top_dir = os.getcwd()
import h5py

def get_atom(usr_inp):
    int_set = []
    elem_list = []
    a = re.findall(r'[A-Z][a-z]*|\d+|\(|\)', usr_inp)
    if len(a) == 1:
        elem_list = a
        int_set = [1]
    else:
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


def get_geom(input,list_nbr):
    cmp=input
    Geom = []
    os.chdir(top_dir + '/geom')
    if len(list_nbr) >= 2:
        with open('{}.geom.txt'.format(cmp)) as f:
            for line in f:
                Geom.append(line)
    else:
        Geom.append('{} 0.0 0.0 0.0'.format(cmp))
    os.chdir('..')
    return Geom

def get_symmetry(input):
    if input == 'HgO':
        return 'c2v'

def Data():
    os.chdir(top_dir)
    mult_info = pd.read_excel('multiplicty.xlsx')
    mult = mult_info.loc[:,'Mult']
    froz_info = pd.read_excel('Frozen_orbital.xlsx')
    froz = froz_info.loc[:,'Freeze']

    maindb = pd.HDFStore('input_gendb.h5','a')
    maindb['multiplicity'] = mult
    maindb['Freeze'] = froz
    maindb.close()


def get_multi(input):
    cmp=input
    mult_data = pd.read_hdf('input_gendb.h5','multiplicity')
    print(type(mult_data[cmp]))
    if mult_data[cmp] == 1:
        return 'singlet'
    elif mult_data[cmp] == 2:
        return 'doublet'
    elif mult_data[cmp] == 3:
        return 'triplet' 
    elif mult_data[cmp] == 4:
        return 'quartet'
    elif mult_data[cmp] == 5:
        return 'quintet'
    elif mult_data[cmp] == 6:
        return 'sextet'
    elif mult_data[cmp] == 7:
        return 'septet'
    else:
        return 'octet'

def get_freeze(list_atom,list_nbr):
    somme_freeze = 0
    freeze_data = pd.read_hdf('input_gendb.h5','Freeze')
    for i in range(0,len(list_atom)):
        somme_freeze = somme_freeze + int(freeze_data[list_atom[i]]) * int(list_nbr[i])
    return str(somme_freeze)

def write_preambule(compound,method,name_base,base,memory1,list_atom,count,ecp):
    Main_file=[]

    if count==0:
        Main_file.append('start {}.{}\n'.format(compound,method))
        Main_file.append('title "{} {}"\n'.format(compound,method))
        if method=='cbs':
            Main_file.append('memory {} m \n'.format(memory1))
        else:
            Main_file.append('memory stack 2400 mb head 200 mb global 1200 mb\n')
        Main_file.append('echo\n\n')
        Main_file.append('geometry\n')
        Main_file.append('symmetry {}\n'.format(get_symmetry(compound)))
        GEOM = get_geom(compound,list_nbr)
        for j in range(0,len(GEOM)):
            Main_file.append(GEOM[j])
        Main_file.append('\n')
    #----------Basis-----------#
    if method=='cbs':
        Main_file.append('basis {} spherical\n'.format(name_base))
        for i in range(0,len(list_atom)):
            if list_atom[i] == 'Hg' or list_atom[i] == 'W' or list_atom[i] == 'Hf' or list_atom[i] == 'Ta' or list_atom[i] == 'Br' or list_atom[i] == 'I':
                Main_file.append(' {} library aug-cc-pv{}z-pp\n'.format(list_atom[i],base))
            else:
                Main_file.append(' {} library aug-cc-pv{}z\n'.format(list_atom[i],base))
        Main_file.append('end\n\n')
    
    elif method=='cc':
        Main_file.append('basis spherical\n')
        for i in range(0,len(list_atom)):
            if list_atom[i] == 'Hg' or list_atom[i] == 'W' or list_atom[i] == 'Hf' or list_atom[i] == 'Ta' or list_atom[i] == 'Br' or list_atom[i] == 'I':
                Main_file.append(' {} library aug-cc-pvtz-pp\n'.format(list_atom[i]))
            else:
                Main_file.append(' {} library aug-cc-pvtz\n'.format(list_atom[i]))
        Main_file.append('end\n\n')
    
    else:
        Main_file.append('basis spherical print\n')
        for i in range(0,len(list_atom)):
            if list_atom[i] == 'Hg' or list_atom[i] == 'W' or list_atom[i] == 'Hf' or list_atom[i] == 'Ta' or list_atom[i] == 'Br' or list_atom[i] == 'I':
                Main_file.append('{} library cc-pwCVDZ-PP\n'.format(list_atom[i]))
                Main_file.append('{} library aug-cc-pvdz-pp_diffuse\n'.format(list_atom[i]))
            else:
                Main_file.append('{} library aug-cc-pwCVDZ\n'.format(list_atom[i]))
        Main_file.append('end\n\n')

    #----------ECP-------------#
    if ecp=='yes':
        Main_file.append('ecp\n')
        for i in range(0,len(list_atom)):
            if list_atom[i] == 'Hg' or list_atom[i] == 'W' or list_atom[i] == 'Hf' or list_atom[i] == 'Ta' or list_atom[i] == 'Br' or list_atom[i] == 'I':
                Main_file.append(' {} library "stuttgart-koeln_mcdhf_rsc_ecp"\n'.format(list_atom[i]))
        Main_file.append('end\n\n')

    return Main_file

def write_calculation(compound,method,name_base,list_atom,list_nbr):
    Main_file=[]

    if method=='cbs':
        Main_file.append('set "ao basis" {}\n\n'.format(name_base))
        Main_file.append('scf\n')
        Main_file.append(' {}\n'.format(get_multi(compound)))
        Main_file.append(' rohf\n')
        if name_base == 'double':
            Main_file.append(' vectors input atomic output {}.{}.movecs\n'.format(compound,name_base))
        elif name_base=='triple':
            Main_file.append(' vectors input project double {}.double.movecs output {}.{}.movecs\n'.format(compound,compound,name_base))
        else:
            Main_file.append(' vectors input project triple {}.triple.movecs output {}.{}.movecs\n'.format(compound,compound,name_base))
        Main_file.append('end\n\n')
        Main_file.append('tce\n')
        Main_file.append(' scf\n')
        Main_file.append(' mp2\n')
        Main_file.append(' io ga\n')
        Main_file.append(' freeze {}\n'.format(get_freeze(list_atom,list_nbr)))
        Main_file.append('end\n\n')
        Main_file.append('task tce energy\n\n')
    else:
        Main_file.append('scf\n')
        Main_file.append(' {}\n'.format(get_multi(compound)))
        Main_file.append(' rohf\n')
        Main_file.append(' vectors input atomic output {}.{}.movecs\n'.format(compound,method))
        Main_file.append('end\n\n')
        
 
    return Main_file
    

#---------------------------- Get the data---------------------------------#
parser = argparse.ArgumentParser(description='''
Hi I make the generator input file for rp_ccCA for the 5d metal only. To be able to run it you need:
- A folder name 'basis' where are all of your basis for each atom, e.g W (W-aug-cc-pcvdz-PP.txt, W-aug-cc-pVDZ-PP.txt, W-aug-cc-pVQZ-PP.txt, W-aug-cc-pVTZ-PP.txt, W-cc-pVTZ-PP.txt) And your pseudo (W-ECP60MDF.txt). You need 60 for the metal, 10 for Br and 28 for I.
- A folder with all the geometry named 'geom' e.g WCl6 (WCl6.geom.txt)
- An excel file with the multiplicty for each compound (multiplicty.xlsx) and one with the charge of the atom where your removed the number of electrons in the ECP (Chrage.xlsx)

It is not optimized but it works. If you want to add something or improve the code, you are free to do.
''',
epilog='\n Usage: input-gen-molpro.py -c WCl6 \n')
parser.add_argument('-c','--compound', help="Name of compound", required=True)
args = vars(parser.parse_args())
if args['compound']:
        compound= args['compound']
else:
        compound = input("Enter a compound: ")
new_dir = '{}'.format(compound)
if os.path.isdir(new_dir) is False:
        os.mkdir(new_dir)
os.chdir(new_dir)

#-----------------Initialization data------------------#

list_atom=get_atom(compound)[0]
list_nbr=get_atom(compound)[1]
Data()


#-----------------Write the cbs file------------------#
file_out = open (top_dir + '/{}/{}.cbs.nw'.format(compound,compound), 'w', newline='\n')
File=[]

Test = write_preambule(compound,'cbs','double','d',800,list_atom,0,'no')

for i in range(0,len(Test)):
    File.append(Test[i])
File.append('\n')

Test = write_preambule(compound,'cbs','triple','t',200,list_atom,1,'no')

for i in range(0,len(Test)):
    File.append(Test[i])
File.append('\n')

Test = write_preambule(compound,'cbs','quad','q',200,list_atom,1,'yes')

for i in range(0,len(Test)):
    File.append(Test[i])



Test = write_calculation(compound,'cbs','double',list_atom,list_nbr)
for i in range(0,len(Test)):
    File.append(Test[i])
File.append('\n')

Test = write_calculation(compound,'cbs','triple',list_atom,list_nbr)
for i in range(0,len(Test)):
    File.append(Test[i])
File.append('\n')

Test = write_calculation(compound,'cbs','quad',list_atom,list_nbr)
for i in range(0,len(Test)):
    File.append(Test[i])
File.append('\n')



for k in range(0, len(File)):
        file_out.write(str(File[k]))

file_out.close()



#-----------------Write the cc file-------------------#

file_out = open (top_dir + '/{}/{}.cc.nw'.format(compound,compound), 'w', newline='\n')
File=[]

Test = write_preambule(compound,'cc','double','d',200,list_atom,0,'yes')

for i in range(0,len(Test)):
    File.append(Test[i])
File.append('\n')

Test =  write_calculation(compound,'cc','double',list_atom,list_nbr)
for i in range(0,len(Test)):
    File.append(Test[i])
File.append('\n')

File.append('tce\n')
File.append(' scf\n')
File.append(' io ga\n')
File.append(' mbpt2\n')
File.append(' freeze {}\n'.format(get_freeze(list_atom,list_nbr)))
File.append('end\n\n')
File.append('task tce energy\n\n')

File.append('tce\n')
File.append(' scf\n')
File.append(' io ga\n')
File.append(' ccsd(t)\n')
File.append(' freeze {}\n'.format(get_freeze(list_atom,list_nbr)))
File.append(' tilesize 6\n')
File.append(' 2eorb\n')
File.append(' 2emet 15\n')
File.append(' tilesize 6\n')
File.append('end\n\n')
File.append('task tce energy\n')


for k in range(0, len(File)):
        file_out.write(str(File[k]))

file_out.close()


#-----------------Write the cv file-------------------#

file_out = open (top_dir + '/{}/{}.cv.nw'.format(compound,compound), 'w', newline='\n')
File=[]

Test = write_preambule(compound,'cv','double','d',200,list_atom,0,'yes')

for i in range(0,len(Test)):
    File.append(Test[i])
File.append('\n')

Test =  write_calculation(compound,'cv','double',list_atom,list_nbr)
for i in range(0,len(Test)):
    File.append(Test[i])
File.append('\n')

File.append('tce\n')
File.append(' scf\n')
File.append(' io ga\n')
File.append(' ccsd(t)\n')
File.append(' freeze {}\n'.format(get_freeze(list_atom,list_nbr)))
File.append(' 2eorb\n')
File.append(' 2emet 15\n')
File.append(' tilesize 6\n')
File.append('end\n\n')
File.append('set tce:nts t\n')
File.append('task tce energy\n\n')

File.append('tce\n')
File.append(' scf\n')
File.append(' io ga\n')
File.append(' ccsd(t)\n')
File.append(' freeze 0\n')
File.append(' 2eorb\n')
File.append(' 2emet 15\n')
File.append(' tilesize 6\n')
File.append('end\n\n')
File.append('task tce energy\n')

for k in range(0, len(File)):
        file_out.write(str(File[k]))

file_out.close()
