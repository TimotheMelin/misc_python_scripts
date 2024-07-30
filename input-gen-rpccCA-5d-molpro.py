import os, re, string
import sys, argparse
from numpy import *
import pandas as pd
top_dir = os.getcwd()
import h5py


#def nbr_atom(input,list_nbr):
#    cmp=input
#    os.chdir(top_dir + '/geom')
#    NumberOfLine = 0
#    if len(list_nbr) >=2:
#        f = open('%s.geom.txt' % cmp ,'r')
#        NumberOfLine=0
#        for line in f:
#            NumberOfLine += 1
#    else:
#        NumberOfLine = 1
#    os.chdir('..')
#    return NumberOfLine

def nbr_atom(input,list_nbr):
    cmp = input
    os.chdir(top_dir + '/geom')
    if len(list_nbr) >= 2 or int(list_nbr[0]) >= 2:
        file = open('{}.geom.txt'.format(cmp),'r')
        num_lines = len(file.readlines())
    else:
        num_lines = 1
    os.chdir('..')
    return num_lines

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


def ECP(input):
    cmp=input
    ecp=[]
    os.chdir(top_dir + '/basis2')
    for i in range(0,len(elem_list_atom)):
        if elem_list_atom[i] == 'S' or elem_list_atom[i] =='F' or elem_list_atom[i] =='Cl' or elem_list_atom[i] =='O':
            Ecp='False'
        elif elem_list_atom[i] == 'Br':
            with open('{}-ECP10MDF.txt'.format(elem_list_atom[i])) as f:
                for line in f:
                    ecp.append(line)
        elif elem_list_atom[i] == 'I':
            with open('{}-ECP28MDF.txt'.format(elem_list_atom[i])) as f:
                for line in f:
                    ecp.append(line)
        else:
            with open('{}-ECP60MDF.txt'.format(elem_list_atom[i])) as f:
                for line in f:
                    ecp.append(line)
    os.chdir('..')
    return ecp

def get_basis(input,base_type):
    cmp=input
    base=base_type
    BASIS=[]
    os.chdir(top_dir + '/basis2')
    for i in range(0,len(elem_list_atom)):
        if elem_list_atom[i] == 'S' or elem_list_atom[i] == 'Cl' or elem_list_atom[i] == 'F' or elem_list_atom[i] == 'O':
            with open('{}-{}.txt'.format(elem_list_atom[i],base)) as f:
                for line in f:
                    BASIS.append(line)
        else:
            with open('{}-{}-PP.txt'.format(elem_list_atom[i],base)) as f:
                for line in f:
                    BASIS.append(line)
    os.chdir('..') 
    return BASIS  

def get_geom(input,list_nbr):
    cmp=input
    Geom = []
    os.chdir(top_dir + '/geom')
    if len(list_nbr) >= 2 or int(list_nbr[0]) >= 2:
        with open('{}.geom.txt'.format(cmp)) as f:
            for line in f:
                Geom.append(line)
    else:
        Geom.append('{} 0.0 0.0 0.0'.format(cmp))
    os.chdir('..')
    return Geom

def Data():
    os.chdir(top_dir)
    mult_info = pd.read_excel('multiplicty.xlsx')
    charge_info= pd.read_excel('Chrage.xlsx',index_col=0)
    mult = mult_info.loc[:,'Mult']
    charge = charge_info.loc[:,'charge']
    
    maindb = pd.HDFStore('input_gendb.h5','a')
    maindb['multiplicity'] = mult
    maindb['Charge']=charge
    maindb.close()

def get_multi(input):
    cmp=input
    mult_data = pd.read_hdf('input_gendb.h5','multiplicity')
    Multiplicity = str(int(mult_data[cmp])-1)
    return Multiplicity


def calculate_charge(input):
    cmp=input
    charge_data = pd.read_hdf('input_gendb.h5','Charge')
    charge_list=[]
    for i in range(0,len(list_atom)):
        if list_atom[i] not in charge_data.index:
            print('No charge define for this atom')
        else:
            charge_list.append(charge_data[list_atom[i]])
    CHARGE=0
    if len(charge_list) == len(list_nbr):
        for i in range(0,len(list_atom)):
            CHARGE=CHARGE+charge_list[i]*list_nbr[i]
    elif len(list_nbr) == 0:
        for i in range(0,len(list_atom)):
            CHARGE=CHARGE+charge_list[i]
    else:
        list_nbr2= [1]+list_nbr
        for i in range(0,len(list_atom)):
            CHARGE=CHARGE+charge_list[i]*list_nbr2[i]
    return CHARGE



def restricted_or_not(Charge):
    Value=Charge
    if Value%2==0:
        Coupled='ccsd(t)'
        Mp2='mp2'
    else:
        Coupled='rccsd(t)'
        Mp2='rmp2'
    return Coupled,Mp2

def write_file(compound,method,base,memory,list_nbr,count):
    Main_file=[]
    
   
    if count==0:
        Main_file.append('***,{}.{}\n'.format(compound,method))
        Main_file.append('memory,{},m \n\n'.format(memory))    
        Main_file.append('geomtyp=xyz\n')
        Main_file.append('geometry={\n')
        Main_file.append(nbr_atom(compound,list_nbr))
        Main_file.append('\n\n')
        GEOM = get_geom(compound,list_nbr)
        for j in range(0,len(GEOM)):
            Main_file.append(GEOM[j])

        Main_file.append('\n}\n\n')

    elem_list_atom=get_atom(compound)[0]
    Base='{}'.format(base)
    Main_file.append('basis={\n')
    Main_file.append('!')
    Main_file.append(Base)
    Basis_set = get_basis(compound,base)
    for j in range(0,len(Basis_set)):
        Main_file.append(Basis_set[j])
    Pseudo = ECP(compound)
    for j in range(0,len(Pseudo)):
        Main_file.append(Pseudo[j])

    return Main_file 


#Arguments 
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



#-----------------Function data -------------------#
elem_list_atom=get_atom(compound)[0]
list_atom=get_atom(compound)[0]
list_nbr=get_atom(compound)[1]
Data()
#-----------Write the cbs file---------------------#

file_out = open (top_dir + '/{}/{}.cbs.com'.format(compound,compound), 'w', newline='\n')
File=[]
Test = write_file(compound,'cbs','aug-cc-pVDZ',200,list_nbr,0)

for i in range(0,len(Test)):
    File.append(Test[i])

File.append('}\n\n')
File.append('{rhf;wf,,,')
File.append(get_multi(compound))
File.append(';save,2100.2}\n')
File.append('rmp2\n')
File.append('e1=energy\n\n')

Test = write_file(compound,'cbs','aug-cc-pVTZ',200,list_nbr,1)

for i in range(0,len(Test)):
    File.append(Test[i])

File.append('}\n\n')
File.append('{rhf;wf,,,')
File.append(get_multi(compound))
File.append(';start,2100.2;save,2101.2}\n')
File.append('rmp2\n')
File.append('e2=energy\n\n')

Test = write_file(compound,'cbs','aug-cc-pVQZ',200,list_nbr,1)

for i in range(0,len(Test)):
    File.append(Test[i])
File.append('}\n\n')
File.append('{rhf;wf,,,')
File.append(get_multi(compound))
File.append(';start,2101.2;save,2102.2}\n')
File.append('rmp2\n')
File.append('e3=energy\n\n')
File.append('---')

for k in range(0, len(File)):
        file_out.write(str(File[k]))

file_out.close()

#------------Write the pcv file part 1-------------#

file_out = open (top_dir + '/{}/{}.pcv1.com'.format(compound,compound), 'w', newline='\n')
File = []
File = write_file(compound,'pcv1','aug-cc-pcvdz',1000,list_nbr,0)

File.append('}\n\n')
File.append('{rhf;wf,,,')
File.append(get_multi(compound))
File.append(';start,2100.2}\n')
File.append('{rccsd(t);}\n')
File.append('e1=energy\n')
File.append('---')

for k in range(0, len(File)):
        file_out.write(str(File[k]))

file_out.close()


#------------Write the pcv file part 2-------------#

file_out = open (top_dir + '/{}/{}.pcv2.com'.format(compound,compound), 'w', newline='\n')
File = []
File = write_file(compound,'pcv2','aug-cc-pcvdz',1500,list_nbr,0)

File.append('}\n\n')
File.append('{rhf;wf,,,')
File.append(get_multi(compound))
File.append(';start,2100.2}\n')
File.append('{rccsd(t);core,0,0,0,0}\n')
File.append('e1=energy\n')
File.append('---')

for k in range(0, len(File)):
        file_out.write(str(File[k]))

file_out.close()
 

#---------------Write the cc file----------------#

file_out = open (top_dir + '/{}/{}.cc.com'.format(compound,compound), 'w', newline='\n')
File = []
File = write_file(compound,'cc','cc-pVTZ',200,list_nbr,0)
#'''
#}
#{rhf;wf,,,'
#;start
#'''.format()

File.append('}\n\n')
File.append('{rhf;wf,,,')
File.append(get_multi(compound))
File.append(';start,2100.2}\n')
File.append(restricted_or_not(calculate_charge(compound))[1])
File.append('\ne1=energy\n')
File.append(restricted_or_not(calculate_charge(compound))[0])
File.append('\ne2=energy\n\n')
File.append('---')

for k in range(0, len(File)):
	file_out.write(str(File[k]))

file_out.close()
