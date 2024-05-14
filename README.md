# Compare_two_diploids_with_two_tetraploids

# ***** THIS IS FOR NOT TOO LARGE FILES. SEE THE END FOR LARGE FILES*********
# Prepare input file
```
 module load tabix
 bgzip -c file.vcf > file.vcf.gz
 tabix -p vcf file.vcf.gz
```
#  Now use vcftools to make a tab delimited file:
```
 module load StdEnv/2020 vcftools/0.1.16
 zcat file.vcf.gz | vcf-to-tab > out.tab
```
Then used this python script to process tab file
# Run it like this

You can find example files for each of the commandas below

```bash
python3 compare_pops_final.py Testing_2.tab my_samp_list_with_sex_and_pop_2.txt comparison_list_2.txt summary.txt
```

```python
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  6 11:11:59 2024

@author: Tharindu
"""

import array
print('\n\n***************************\n')
print('Comparing population specific sites')
# set path to the current folder
print('Please check the working directory before running')
print('\n\n***************************\n')
import time
time.sleep(2.5)
import os as os

import inspect
import pandas as pd
import numpy as np
import csv
import ast
from numpy.random import randint
from itertools import chain
import sys as sys

def get_script_path():
    return os.path.dirname(os.path.realpath(sys.argv[0]))

# to get the current working directory
directory = "/Users/Tharindu/Library/CloudStorage/OneDrive-McMasterUniversity/for_lab_and_research/Tharindu_on_Mac/lab/python_projects/hetero_comparison_4_pops"
#directory = r"C:\Users\thari\OneDrive - McMaster University\for_lab_and_research\Tharindu_on_Mac\lab\python_projects\hetero_comparison_4_pops"
os. chdir(directory)
print('working in', directory)
print('\n\n\n\n\n')


inputfile = sys.argv[1]

print('printing first 10 lines of the file')


def open_file(inputfile):
    try:
        with open(inputfile, 'r') as file:
            content = file.read()
            # Do something with the file content (e.g., print it)
            return content
            print(content, 10)
    except FileNotFoundError:
        print(f"Error: File '{inputfile}' not found.")
    except Exception as e:
        print(f"An error occurred: {e}")


file_contents = open_file(inputfile)


#reading first line
tab_header=file_contents.split('\n', 1)[0]

# splitting 
samp_list=tab_header.split('\t')

print('your sample list is')
my_samp_list=samp_list[3:len(samp_list)]
print(my_samp_list)

file = open('my_samp_list.txt','w')
for item in my_samp_list:
	file.write(item+"\n")
file.close()

print('\n\n\n\n\n')

print('Your sample list was created as my_samp_list.txt. Please open with excel and add the corresponding genders as M or F in second column and population names in the third column***AND SAVE IT and enter the file name with full path')

file_with_sex_path = sys.argv[2]

# creating needed dataframes



file_with_sex=pd.read_table(file_with_sex_path,header=None)

file_with_sex.columns=['file','sex','pop']

print('Input was')

print(file_with_sex)

print('\n\n***************************\n')
print('Creating seperate dataframes for different populations')
print('\n\n***************************\n')
time.sleep(2.5)

# get a list of names
names=file_with_sex['pop'].unique().tolist()

combination_file=sys.argv[3]

#get the summary file name
summary_saving_name=sys.argv[4]

combination_list=pd.read_table(combination_file)
combination_list.columns=['p1','p2','p3','p4']
        
# Create df to store data
summary_df = pd.DataFrame(columns=['p1','p2','p3','p4','b_shared_with_c','b_shared_with_d'])

for combo in range(0,len(combination_list)) :
        
    
        p1=combination_list['p1'][combo]
        p2=combination_list['p2'][combo]
        p3=combination_list['p3'][combo]
        p4=combination_list['p4'][combo]
        
        if p1!=p2 and p1!=p3 and p1!=p4 and p2!=p3 and p2!=p4 and p3!=p4:
            
            print('Running iteration '+str(combo+1)+' out of total '+str(len(combination_list)))
        
            #avoid using same individual as different pops
            
            print('calculating similarities and differences for '+str(len(combination_list))+ ' combinations. \n\n *****This might take hours depending on the number of populations******')
            print(p1,p2,p3,p4)
            
            print('Reading input file')
            
            # read selected columns
            file_with_locs=pd.read_csv(inputfile, sep='\t',usecols=['#CHROM','POS',p1,p2,p3,p4])
            
            #rename
            file_with_locs.columns=['chr','pos','p1','p2','p3','p4']
            
            # print('saving temp csv file')
            # #remove headers
            # file_with_locs.to_csv('temp.tab', sep='\t', index=False)
            
            
            
            print('\n\n***************************\n')
            print('Removing empties')
            print('\n\n***************************\n')
            #time.sleep(2.5)
            file_with_locs = file_with_locs[file_with_locs.ne('./.').all(1)]
            file_with_locs = file_with_locs[file_with_locs.ne('*').all(1)]
            
            # ***** pop1******
            
            print('pop1')
            
            #extracting pop1
            pop1_all_sites=file_with_locs.iloc[:, [0,1,2]]
            
            #name
            #colname_pop1=str(pop1_all_sites.columns[2])
            
            #splitting it to 2 columns
            pop1_all_sites[['var1','var2']] = pop1_all_sites.p1.str.split('/',expand=True)
            
            # Delete heterozygous rows 
            pop1_homo_only = pop1_all_sites.drop(pop1_all_sites[pop1_all_sites['var1'] != pop1_all_sites['var2']].index)
            
            #select columns
            
            selected_pop1=pop1_homo_only.iloc[:, [0,1,3]]
            
            #saving an example file to proofread
            
            #file_with_locs.to_csv('file_with_locs.tab', sep='\t', index=False)
            #pop1_all_sites.to_csv('pop1_all_sites.tab', sep='\t', index=False)
            #selected_pop1.to_csv('selected_pop_1.tab', sep='\t', index=False)
            
            
            print('Removed '+str(pop1_all_sites.shape[0]-pop1_homo_only.shape[0])+' heterozygous sites')
            
            # ****************
            
            # ***** pop2******
            
            print('pop2')
            
            #extracting pop2
            pop2_all_sites=file_with_locs.iloc[:, [0,1,3]]
            
            #name
            #colname_pop2=str(pop2_all_sites.columns[2])
            
            #splitting it to 2 columns
            pop2_all_sites[['var1','var2']] = pop2_all_sites.p2.str.split('/',expand=True)
            
            # Delete heterozygous rows 
            pop2_homo_only = pop2_all_sites.drop(pop2_all_sites[pop2_all_sites['var1'] != pop2_all_sites['var2']].index)
            
            #select columns
            
            selected_pop2=pop2_homo_only.iloc[:, [0,1,3]]
            
            
            print('Removed '+str(pop2_all_sites.shape[0]-pop2_homo_only.shape[0])+' heterozygous sites')
            
            # ****************
            
            # ***** pop3******
            
            print('pop3')
            
            #extracting pop3
            pop3_all_sites=file_with_locs.iloc[:, [0,1,4]]
            
            #name
            #colname_pop3=str(pop3_all_sites.columns[2])
            
            #splitting it to 2 columns
            pop3_all_sites[['var1','var2']] = pop3_all_sites.p3.str.split('/',expand=True)
            
            # Delete heterozygous rows 
            #pop3_homo_only = pop3_all_sites.drop(pop3_all_sites[pop3_all_sites['var1'] != pop3_all_sites['var2']].index)
            
            #select columns
            
            selected_pop3=pop3_all_sites.iloc[:, [0,1,3,4]]
            
            #saving an example file to proofread
            
            #selected_pop3.to_csv('selected_pop3.tab', sep='\t', index=False)
            
            print('loaded pop 3')
            
            # ****************
            
            # ***** pop4******
            
            print('pop4')
            
            #extracting pop3
            pop4_all_sites=file_with_locs.iloc[:, [0,1,5]]
            
            #name
            #colname_pop3=str(pop3_all_sites.columns[2])
            
            #splitting it to 2 columns
            pop4_all_sites[['var1','var2']] = pop4_all_sites.p4.str.split('/',expand=True)
            
            # Delete heterozygous rows 
            #pop4_homo_only = pop4_all_sites.drop(pop4_all_sites[pop4_all_sites['var1'] != pop4_all_sites['var2']].index)
            
            #select columns
            
            selected_pop4=pop4_all_sites.iloc[:, [0,1,3,4]]
            
            
            print('loaded pop 4')
            
            # ****************
            
            # Combine identifiers to make it easier to compare - pop1***
            
            #first I have to convert values to string
            selected_pop1['chr'] = selected_pop1['chr'].astype(str)
            selected_pop1['pos'] = selected_pop1['pos'].astype(str)
            
            #combine
            selected_pop1['Chr-pos'] = selected_pop1['chr'].astype(str) + selected_pop1['pos']
            
            #drop extras
            selected_pop1.drop(['chr', 'pos'], axis=1)
            
            # change column order
            
            selected_pop1 = selected_pop1[['Chr-pos', 'var1']]
            
            #rename columns
            
            selected_pop1.columns=['Chr-pos', 'p1']
            
            
            # **********

            # ****************
            
            # Combine identifiers to make it easier to compare - pop2***
            
            #first I have to convert values to string
            selected_pop2['chr'] = selected_pop2['chr'].astype(str)
            selected_pop2['pos'] = selected_pop2['pos'].astype(str)
            
            #combine
            selected_pop2['Chr-pos'] = selected_pop2['chr'].astype(str) + selected_pop2['pos']
            
            #drop extras
            selected_pop2.drop(['chr', 'pos'], axis=1)
            
            # change column order
            
            selected_pop2 = selected_pop2[['Chr-pos', 'var1']]
            
            #rename columns
            
            selected_pop2.columns=['Chr-pos', 'p2']
            
            
            # **********

            # ****************
            
            # Combine identifiers to make it easier to compare - pop3***
            
            #first I have to convert values to string
            selected_pop3['chr'] = selected_pop3['chr'].astype(str)
            selected_pop3['pos'] = selected_pop3['pos'].astype(str)
            
            #combine
            selected_pop3['Chr-pos'] = selected_pop3['chr'].astype(str) + selected_pop3['pos']
            
            #drop extras
            selected_pop3.drop(['chr', 'pos'], axis=1)
            
            # change column order
            
            selected_pop3 = selected_pop3[['Chr-pos','var1','var2']]
            
            #rename columns
            
            selected_pop3.columns=['Chr-pos', 'p3_var1','p3_var2']
            
            
            # **********
            
            # ****************
            
            # Combine identifiers to make it easier to compare - pop4***
            
            #first I have to convert values to string
            selected_pop4['chr'] = selected_pop4['chr'].astype(str)
            selected_pop4['pos'] = selected_pop4['pos'].astype(str)
            
            #combine
            selected_pop4['Chr-pos'] = selected_pop4['chr'].astype(str) + selected_pop4['pos']
            
            #drop extras
            selected_pop4.drop(['chr', 'pos'], axis=1)
            
            # change column order
            
            selected_pop4 = selected_pop4[['Chr-pos','var1','var2']]
            
            #rename columns
            
            selected_pop4.columns=['Chr-pos', 'p4_var1','p4_var2']
            
            
            # **********
            
            # **********
            
            #combine all df s
            
            final_dataset=pd.concat([selected_pop1.set_index('Chr-pos'),selected_pop2.set_index('Chr-pos')], axis=1, join='inner').reset_index()
            
            final_dataset=pd.concat([final_dataset.set_index('Chr-pos'),selected_pop3.set_index('Chr-pos')], axis=1, join='inner').reset_index()
            
            final_dataset=pd.concat([final_dataset.set_index('Chr-pos'),selected_pop4.set_index('Chr-pos')], axis=1, join='inner').reset_index()
            # ******
            #remove deletions (*)
            
            final_dataset = final_dataset[final_dataset.ne('*').all(1)]
            
            #uncomment this if you wanna save csv to check manually
            
            #final_dataset.to_csv('finalized_data_for_counting.tsv', sep="\t")
            
            #first I have to convert and remove NAa
            
            print('\n*********************\n')
            print('Counting similarities and differences.....\n This may take a while \n Please wait...')
            print('\n*********************\n')
            
            final_dataset= final_dataset.replace('NaN', np.nan)
            
            # drop rows with more that 2 different locus types (one value goes to Chr-pos. So dropping values>3)
            
            #count uniques
            final_dataset['num_uniq'] = [len(set(v[pd.notna(v)].tolist())) for v in final_dataset.values]
            
            
            #drop
            final_dataset = final_dataset.drop(final_dataset[final_dataset.num_uniq>3].index)
                
            
            # ************
            
            b_eq_c_sub1=final_dataset.apply(lambda x: x.p2 == x.p3_var1 != x.p1 and x.p2 !=x.p3_var2 and x.p2 !=x.p4_var1 and x.p2 !=x.p4_var2 and x.p3_var2==x.p1, axis = 1).value_counts().get(True, 0)
            b_eq_c_sub2=final_dataset.apply(lambda x: x.p2 == x.p3_var2 != x.p1 and x.p2 !=x.p3_var1 and x.p2 !=x.p4_var1 and x.p2 !=x.p4_var2 and x.p3_var1==x.p1, axis = 1).value_counts().get(True, 0)
            
            b_shared_with_c=b_eq_c_sub1+b_eq_c_sub2
            
             
            # for readable copy*******
            
            # b_eq_c_sub1=final_dataset.apply(lambda x: x.p2 == x.p3_var1 != x.p1 !=x.p4, axis = 1)
            # b_eq_c_sub2=final_dataset.apply(lambda x: x.p2 == x.p3_var2 != x.p1 !=x.p4, axis = 1)
            
            # b_eq_c_sub1.to_csv('b_eq_c_sub1_db.tsv', sep="\t")
            # b_eq_c_sub2.to_csv('b_eq_c_sub2_db.tsv', sep="\t")
            
            # ************
            
            # ************
            
            b_eq_d_sub1=final_dataset.apply(lambda x: x.p2 == x.p4_var1 != x.p1 and x.p2 !=x.p4_var2 and x.p2 !=x.p3_var1 and x.p2 !=x.p3_var2 and x.p4_var2==x.p1, axis = 1).value_counts().get(True, 0)
            b_eq_d_sub2=final_dataset.apply(lambda x: x.p2 == x.p4_var2 != x.p1 and x.p2 !=x.p4_var1 and x.p2 !=x.p3_var1 and x.p2 !=x.p3_var2 and x.p4_var1==x.p1,axis = 1).value_counts().get(True, 0)
            
            b_shared_with_d=b_eq_d_sub1+b_eq_d_sub2
            
            
            # for readable copy*******
            
            # b_eq_d_sub1=final_dataset.apply(lambda x: x.p2 == x.p4_var1 != x.p1 !=x.p3, axis = 1)
            # b_eq_d_sub2=final_dataset.apply(lambda x: x.p2 == x.p4_var2 != x.p1 !=x.p3, axis = 1)
            
            # b_eq_d_sub1.to_csv('b_eq_d_sub1_db.tsv', sep="\t")
            # b_eq_d_sub2.to_csv('b_eq_d_sub2_db.tsv', sep="\t")
            
            # ************
            
            print("b_shared_with_c in "+str(b_shared_with_c)+'\n'+'b_shared_with_d in '+str(b_shared_with_d))
            
            to_write=p1,p2,p3,p4,b_shared_with_c,b_shared_with_d
            
            summary_df.loc[len(summary_df.index)] = to_write
            
            # saving as tsv file 
            summary_df.to_csv(summary_saving_name, sep="\t",index=False)

```

# testing_2.tab
```txt
#CHROM	POS	REF	F_Ghana_WZ_BJE4687_combined__sorted.bam	F_IvoryCoast_xen228_combined__sorted.bam	XT7_WY_no_adapt__sorted.bam	all_ROM19161_sorted.bam	all_calcaratus_sorted.bam	mello_GermSeq_sorted.bam
Chr10	707	A/A	A/A	A/A	A/A	A/A	A/A	A/A
Chr10	708	A/A	A/A	A/A	T/T	A/A	T/A	A/A
Chr10	709	A/A	A/A	A/A	T/T	A/A	T/A	A/A
Chr10	710	A/A	A/A	A/A	T/T	A/A	T/A	T/A
Chr10	711	A/A	T/T	A/A	T/T	A/A	T/A	A/A
Chr10	712	A/A	T/T	A/A	T/T	A/A	T/A	A/G
```

# my_samp_list_with_sex_and_pop_2.txt
```txt
F_Ghana_WZ_BJE4687_combined__sorted.bam	F	Tropicalis
F_IvoryCoast_xen228_combined__sorted.bam	F	other
XT7_WY_no_adapt__sorted.bam	M	other
all_ROM19161_sorted.bam	F	Liberia
all_calcaratus_sorted.bam	F	Calcaratus
mello_GermSeq_sorted.bam	F	other
```
# comparison_list_2.txt

```txt
1	2	3	4
F_Ghana_WZ_BJE4687_combined__sorted.bam	XT7_WY_no_adapt__sorted.bam	all_calcaratus_sorted.bam	mello_GermSeq_sorted.bam
```
# summary.txt is just a name for summary file

# **************** FOR LARGE FILES **********************

First I splitted tab file into smaller files

```bash
split --numeric-suffixes=1 -l 6 combined_chrs.tab "combined_chrs_part"
```
then added headers
```
for i in combined_chrs_part*; do echo -e "#CHROM\tPOS\tREF\tF_Ghana_WZ_BJE4687_combined__sorted.bam\tF_IvoryCoast_xen228_combined__sorted.bam\tF_Nigeria_EUA0331_combined__sorted.bam\tF_Nigeria_EUA0333_combined__sorted.bam\tF_SierraLeone_AMNH17272_combined__sorted.bam\tF_SierraLeone_AMNH17274_combined__sorted.bam\tJBL052_concatscafs_sorted.bam\tM_Ghana_WY_BJE4362_combined__sorted.bam\tM_Ghana_ZY_BJE4360_combined__sorted.bam\tM_Nigeria_EUA0334_combined__sorted.bam\tM_Nigeria_EUA0335_combined__sorted.bam\tM_SierraLeone_AMNH17271_combined__sorted.bam\tM_SierraLeone_AMNH17273_combined__sorted.bam\tXT10_WZ_no_adapt._sorted.bam\tXT11_WW_trim_no_adapt_scafconcat_sorted.bam\tXT1_ZY_no_adapt._sorted.bam\tXT7_WY_no_adapt__sorted.bam\tall_ROM19161_sorted.bam\tall_calcaratus_sorted.bam\tmello_GermSeq_sorted.bam" | cat - ${i} >/tmp/out && mv /tmp/out ${i} ;done
```
Then I had to rename first 9 files to get rid of '0' at the begining of combination files

```
mv file01 file1
```
Then I made a directory for outputs
```
mkdir outputs
```
I used lists of individuals to be used as each of the populations and python script below to come up with all possible combinations

example 
s_list_1.txt

```txt
F_Ghana_WZ_BJE4687_combined__sorted.bam
F_IvoryCoast_xen228_combined__sorted.bam
F_Nigeria_EUA0331_combined__sorted.bam
F_Nigeria_EUA0333_combined__sorted.bam
F_SierraLeone_AMNH17272_combined__sorted.bam
F_SierraLeone_AMNH17274_combined__sorted.bam
M_Ghana_WY_BJE4362_combined__sorted.bam
M_Ghana_ZY_BJE4360_combined__sorted.bam
M_Nigeria_EUA0334_combined__sorted.bam
M_Nigeria_EUA0335_combined__sorted.bam
M_SierraLeone_AMNH17271_combined__sorted.bam
M_SierraLeone_AMNH17273_combined__sorted.bam
XT10_WZ_no_adapt._sorted.bam
XT11_WW_trim_no_adapt_scafconcat_sorted.bam
XT1_ZY_no_adapt._sorted.bam
XT7_WY_no_adapt__sorted.bam
all_ROM19161_sorted.bam
```
combination_generator.py

```python
# -*- coding: utf-8 -*-
"""
Created on Wed May  8 11:48:45 2024

@author: Terry
"""

# get all the possible populatiopn combinations

import array
import os
import inspect
import os
import pandas as pd
import numpy as np
import csv
import ast
from numpy.random import randint
from itertools import chain
import sys
import itertools

# to get the current working directory
#directory = r"C:\Users\thari\OneDrive - McMaster University\for_lab_and_research\Tharindu_on_Mac\lab\python_projects\hetero_comparison_4_pops\pop_combos"
directory = "/Users/Tharindu/Library/CloudStorage/OneDrive-McMasterUniversity/for_lab_and_research/Tharindu_on_Mac/lab/python_projects/hetero_comparison_4_pops/pop_combos"
os. chdir(directory)
print('working in', directory)
print('\n\n\n\n\n')


# ask for population lists to use as pop1 , 2 and 3

s_list_1=input('\n Please enter the file name with first population list\n')
s_list_2=input('\n Please enter the file name with second population list\n')
s_list_3=input('\n Please enter the file name with third population list\n')
s_list_4=input('\n Please enter the file name with third population list\n')


file1 = open(s_list_1, "r")
list1=list(csv.reader(file1, delimiter="\t"))
list1=list(chain.from_iterable(list1))

file2 = open(s_list_2, "r")
list2=list(csv.reader(file2, delimiter="\t"))
list2=list(chain.from_iterable(list2))


file3 = open(s_list_3, "r")
list3=list(csv.reader(file3, delimiter="\t"))
list3=list(chain.from_iterable(list3))

file4 = open(s_list_4, "r")
list4=list(csv.reader(file4, delimiter="\t"))
list4=list(chain.from_iterable(list4))

pop_lists=[list1,list2,list3,list4]

#pop_lists=[my_samp_list,my_samp_list,my_samp_list]


#get all combinations
all_pop_combinations=list(itertools.product(*pop_lists))

combination_list=all_pop_combinations


# #****** This filtering is only if popa*popb is same as popb*popa
# #remove duplicates
# combination_list=list({*map(tuple, map(sorted, all_pop_combinations))})

# *********************************************

# #convert list to df
combination_list=pd.DataFrame(combination_list)

# #****** This filtering is only if popa*popb is same as popb*popa
# #drop cols rows with equal vals
# combination_list = combination_list.drop(combination_list[combination_list[0] == combination_list[1]].index)


# *********************************************
print('Saving a list of all possible combinations')


combination_list.to_csv('combination_list.tsv',sep='\t',index=False)
```
Then I had to split the resulting combination file into smaller ones
```bash
split --numeric-suffixes=1 -l 6 combination_list.tsv "combination"
```
and added headers
```
for i in combination*; do echo -e "0\t1\t2\t3" | cat - ${i} >/tmp/out && mv /tmp/out ${i} ;done
```

# *** And removed the extra header in Combination1

Then I wrote 11 different bash scripts to run all individual combinations with each of the tab files 

Example 

run_comparisons_part_1.sh

```bash
#!/bin/sh
#SBATCH --job-name=bwa_505
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=64gb
#SBATCH --output=comp_pop.%J.out
#SBATCH --error=comp_pop.%J.err
#SBATCH --account=def-ben
#SBATCH --array=1-49

#SBATCH --mail-user=premacht@mcmaster.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL

module load StdEnv/2020 python/3.8.2
virtualenv --no-download ~/ENV
source ~/ENV/bin/activate
pip install --no-index --upgrade pip
pip install pandas --no-index

python compare_pops_final.py ./../combined_chrs_more_efficient/combined_chrs_part01 ./../combined_chrs_more_efficient/my_sample_list_with_sex_and_pop.txt ./combos/combination${SLURM_ARRAY_TASK_ID} ./summary_part1_${SLURM_ARRAY_TASK_ID}
```
then ran them all 
``` bash
for i in $(seq 1 11); do sbatch run_comparisons_part_${i}.sh;done
```








