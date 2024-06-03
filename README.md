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
            
            # rearrange populations to make sure they are in the proper order
            file_with_locs = file_with_locs.reindex(['#CHROM','POS',p1,p2,p3,p4], axis=1)
            
            # save to check
            file_with_locs.to_csv('file_with_locs',sep='\t',index=False)
           
            
            #rename
            file_with_locs.columns=['chr','pos','p1','p2','p3','p4']
            
            # print('saving temp csv file')
            
            
            
            
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
            
            
            b_eq_d_sub1=final_dataset.apply(lambda x: x.p2 == x.p4_var1 != x.p1 and x.p2 !=x.p4_var2 and x.p2 !=x.p3_var1 and x.p2 !=x.p3_var2 and x.p4_var2==x.p1, axis = 1).value_counts().get(True, 0)
            b_eq_d_sub2=final_dataset.apply(lambda x: x.p2 == x.p4_var2 != x.p1 and x.p2 !=x.p4_var1 and x.p2 !=x.p3_var1 and x.p2 !=x.p3_var2 and x.p4_var1==x.p1,axis = 1).value_counts().get(True, 0)
            
            b_shared_with_d=b_eq_d_sub1+b_eq_d_sub2
            
            # extract selected sites for proof read
            
            # *******THIS IS ONLU FOR TESTING. COMMENT OUT WHEN RUNNING TO INCREASE EFFICIENCY*******
            
            # b_eq_c_sites_1=final_dataset.apply(lambda x: x.p2 == x.p3_var1 != x.p1 and x.p2 !=x.p3_var2 and x.p2 !=x.p4_var1 and x.p2 !=x.p4_var2 and x.p3_var2==x.p1, axis = 1)
            # b_eq_c_sites_2=final_dataset.apply(lambda x: x.p2 == x.p3_var2 != x.p1 and x.p2 !=x.p3_var1 and x.p2 !=x.p4_var1 and x.p2 !=x.p4_var2 and x.p3_var1==x.p1, axis = 1)
            
            # # use or logic to add two lists together
            # b_eq_c_all_sites=np.logical_or(b_eq_c_sites_1,b_eq_c_sites_2)
            
            # b_eq_d_sites_1=final_dataset.apply(lambda x: x.p2 == x.p4_var1 != x.p1 and x.p2 !=x.p4_var2 and x.p2 !=x.p3_var1 and x.p2 !=x.p3_var2 and x.p4_var2==x.p1, axis = 1)
            # b_eq_d_sites_2=final_dataset.apply(lambda x: x.p2 == x.p4_var2 != x.p1 and x.p2 !=x.p4_var1 and x.p2 !=x.p3_var1 and x.p2 !=x.p3_var2 and x.p4_var1==x.p1,axis = 1)
            
            # # use or logic to add two lists together
            # b_eq_d_all_sites=np.logical_or(b_eq_d_sites_1,b_eq_d_sites_2)
            
            # final_dataset["b_eq_c_all_sites"]=b_eq_c_all_sites
            # final_dataset["b_eq_d_all_sites"]=b_eq_d_all_sites
            
            # final_dataset_save_name=p1+"_and_"+p2+"full_final_dataset.txt"
            
            # # saving as tsv file 
            # final_dataset.to_csv(final_dataset_save_name, sep="\t",index=False)
            #***************************************************
            
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
This sript then produces a bunch of splitted output files.
I combined them using with the following command (part 11 is just a very small portion with nothing considerable. Therefore only using 1-10)
```bash
mkdir combined_summaries
for i in {1..10};do tail -q -n +2 summary_part${i}_* > combined_summaries/combined_summary_part_${i}.tab ; done
```
then downloaded summaries
```
rsync -axvH --no-g --no-p premacht@graham.computecanada.ca:/scratch/premacht/python_projects_2023/hetero_comparison_4_pops/combined_summaries .
```
Put them in the combined_summaries folder and plot the results with following R script
```R
library("rstudioapi") 
setwd(dirname(getActiveDocumentContext()$path))
library(forcats)
library(reshape)
library(tidyr)
library(ggplot2)
#install.packages("rlist")
library(rlist)
#install.packages("insight")
library(insight)
#print_color("ERROR", "red")
library(plyr)
library(dplyr)

# first I have to load all data from different parts and add them together

# **** Chnage the 'range_of_parts' to match the no of parts here

#part 11 didn't have anything and was left out
range_of_parts<-c(1:10)
for (i in range_of_parts) {
  assign(paste('my_dataset_',i,sep = ''),read.csv(paste('./combined_summaries/combined_summary_part_',i,'.tab',sep = ''),sep = '\t',header = FALSE,col.names = c('p1','p2','p3','p4','shared_with_cal','shared_with_mello')))
}

# loading a placeholder df as the sum df(to add all values together)

my_data<-my_dataset_3

#set all of its values to 0

my_data$shared_with_cal<-0
my_data$shared_with_mello<-0

#checking for the same order of columns between dfs
for (part in range_of_parts) {
  current_df_name<-paste('my_dataset_',part,sep = '')
  current_df<-get(current_df_name)
  
  checked_list<-purrr::map2_dbl(my_dataset_3$p1, current_df$p1, ~sum(.x == .y))
  allSame <- function(x) length(unique(x)) == 1
  is_similiar<-allSame(checked_list)
  if (is_similiar==TRUE) {
    print(paste('Datase 1 pop1 order is equal to Dataset ',part,'pop1. GOOD TO PROCEED'))
  }
  if(is_similiar==FALSE) {
    print_color(paste('Datase 1 pop1 order is NOT equal to Dataset ',part,'pop1. *****PLEASE CHECK AND RE ARRANGE*****'),'red')
    Sys.sleep(2)
  }
  
  checked_list<-purrr::map2_dbl(my_dataset_3$p2, current_df$p2, ~sum(.x == .y))
  allSame <- function(x) length(unique(x)) == 1
  is_similiar<-allSame(checked_list)
  if (is_similiar==TRUE) {
    print(paste('Datase 1 pop2 order is equal to Dataset ',part,'pop2. GOOD TO PROCEED'))
  }
  if(is_similiar==FALSE) {
    print_color(paste('Datase 1 pop2 order is NOT equal to Dataset ',part,'pop2. *****PLEASE CHECK AND RE ARRANGE*****'),'red')
    Sys.sleep(2)
  }
  
  checked_list<-purrr::map2_dbl(my_dataset_3$p3, current_df$p3, ~sum(.x == .y))
  allSame <- function(x) length(unique(x)) == 1
  is_similiar<-allSame(checked_list)
  if (is_similiar==TRUE) {
    print(paste('Datase 1 pop3 order is equal to Dataset ',part,'pop3. GOOD TO PROCEED'))
  }
  if(is_similiar==FALSE) {
    print_color(paste('Datase 1 pop3 order is NOT equal to Dataset ',part,'pop3. *****PLEASE CHECK AND RE ARRANGE*****'),'red')
    Sys.sleep(2)
  }
}



#sum data for final dataset

for (part in range_of_parts) {
  current_dataset_name<-paste('my_dataset_',part,sep = '')
  current_dataset<-get(current_dataset_name)
  
  my_data$shared_with_cal<-my_data$shared_with_cal+current_dataset$shared_with_cal
  my_data$shared_with_mello<-my_data$shared_with_mello+current_dataset$shared_with_mello
}

#saving a copy to proof read

write.table(my_data,file = './combined_dataframe_unarranged.tsv', sep='\t', col.names = NA)

# renaming all samples
my_data[my_data == "F_SierraLeone_AMNH17272_combined__sorted.bam"] <- "SL_F1"
my_data[my_data == "F_SierraLeone_AMNH17274_combined__sorted.bam"] <- "SL_F2"
my_data[my_data == "M_SierraLeone_AMNH17271_combined__sorted.bam"] <- "SL_M1"
my_data[my_data == "M_SierraLeone_AMNH17273_combined__sorted.bam"] <- "SL_M2"
my_data[my_data == "all_ROM19161_sorted.bam"] <- "LB_F1"
my_data[my_data == "F_IvoryCoast_xen228_combined__sorted.bam"] <- "IC_F1"
my_data[my_data == "F_Ghana_WZ_BJE4687_combined__sorted.bam"] <- "GH_F1"
my_data[my_data == "M_Ghana_WY_BJE4362_combined__sorted.bam"] <- "GH_M1"
my_data[my_data == "M_Ghana_ZY_BJE4360_combined__sorted.bam"] <- "GH_M2"
my_data[my_data == "XT10_WZ_no_adapt._sorted.bam"] <- "LT_F1"
my_data[my_data == "XT11_WW_trim_no_adapt_scafconcat_sorted.bam"] <- "LT_F2"
my_data[my_data == "XT1_ZY_no_adapt._sorted.bam"] <- "LT_M1"
my_data[my_data == "XT7_WY_no_adapt__sorted.bam"] <- "LT_M2"
my_data[my_data == "F_Nigeria_EUA0331_combined__sorted.bam"] <- "NG_F1"
my_data[my_data == "F_Nigeria_EUA0333_combined__sorted.bam"] <- "NG_F2"
my_data[my_data == "M_Nigeria_EUA0334_combined__sorted.bam"] <- "NG_M1"
my_data[my_data == "M_Nigeria_EUA0335_combined__sorted.bam"] <- "NG_M2"


my_data[my_data == "all_calcaratus_sorted.bam"] <- "Cal"
my_data[my_data == "mello_GermSeq_sorted.bam"] <- "Mello"

#looking for the highest value in both comparisons so I can use the same scale for both

max_y_lim<-max(my_data$shared_with_cal,my_data$shared_with_mello)
upper_y<-(round_any(max_y_lim,1000000))/1000000

# use this section to all sample list
sample_list<-c('SL_F1',
               'SL_F2',
               'SL_M1',
               'SL_M2',
               'LB_F1',
               'IC_F1',
               'GH_F1',
               'GH_M1',
               'GH_M2',
               'LT_F1',
               'LT_F2',
               'LT_M1',
               'LT_M2',
               'NG_F1',
               'NG_F2',
               'NG_M1',
               'NG_M2')


# slect data for pattern 1 and 2 for each of the individuals

for (sample in sample_list) {

current_Sample<-sample

pat1_name<-paste(current_Sample,"_pattern1",sep = "")

pat1<-subset(my_data, p1 == current_Sample)

#melt df
pat1_long <- pat1 %>%
  pivot_longer(
    cols = `shared_with_cal`:`shared_with_mello`,
    names_to = "Comparison",
    values_to = "value"
  )


assign(pat1_name,pat1_long)

pat2_name<-paste(current_Sample,"_pattern2",sep = "")
pat2<-subset(my_data, p2 == current_Sample)

#melt df
pat2_long <- pat2 %>%
  pivot_longer(
    cols = `shared_with_cal`:`shared_with_mello`,
    names_to = "Comparison",
    values_to = "value"
  )

assign(pat2_name,pat2_long)

}

# combine all pattern 1 dataframes together

library(stringr)

pattern_df_1<-mget(ls(pattern = '*pattern1')) %>%
  setNames(str_replace(names(.), ".*_(\\d+)$", "pattern\\1")) %>% 
  bind_rows(.id = "pattern") %>%
  relocate(pattern, .after = last_col())

pattern_df_2<-mget(ls(pattern = '*pattern2')) %>%
  setNames(str_replace(names(.), ".*_(\\d+)$", "pattern\\1")) %>% 
  bind_rows(.id = "pattern") %>%
  relocate(pattern, .after = last_col())

#have to switch p1 and p2 on pattern_df_2 to match it with pattern_df_1


colnames(pattern_df_2)<-c("p2","p1","p3","p4","Comparison","value","pattern")

#combine those two big dataframes together

both_pattern_together<-mget(ls(pattern = 'pattern_df_*')) %>%
  setNames(str_replace(names(.), ".*_(\\d+)$", "pattern_type\\1")) %>% 
  bind_rows(.id = "pattern_type") %>%
  relocate(pattern_type, .after = last_col())

# Combine identifiers to make it easier to compare - pop1***

both_pattern_together$total_comparison <- paste(both_pattern_together$Comparison,both_pattern_together$pattern_type,sep = "_")

#use sample list order as levels

both_pattern_together$p1=factor(both_pattern_together$p1,levels = sample_list)
both_pattern_together$p2=factor(both_pattern_together$p2,levels = sample_list)

#plot

square_plot<-ggplot(both_pattern_together, aes(x=p1,y = value/1000000,fill=total_comparison)) +
  geom_bar(position="dodge",stat='identity')+
  facet_wrap(~ p2,nrow = 19)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  xlab("Individual B") + ylab("No. of bp (millions)")+
  scale_fill_manual(values=c( "red","orange","purple","blue"))+
  #ylim(0,upper_y)+
  theme(plot.title = element_text(hjust = 0.5,size = 10))

ggsave("square_plot.pdf",square_plot,width = 12,height = 11)

#subset samples

selected_sample_list<-c('LB_F1',
               'GH_F1',
               'NG_F1')

subset_to_plot <- both_pattern_together[both_pattern_together$p2 %in% selected_sample_list,]

#plot

subset_square_plot<-ggplot(subset_to_plot, aes(x=p1,y = value/1000000,fill=total_comparison)) +
  geom_bar(position="dodge",stat='identity')+
  facet_wrap(~ p2,nrow = 19)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  xlab("Individual B") + ylab("No. of bp (millions)")+
  scale_fill_manual(values=c( "red","orange","purple","blue"))+
  ylim(0,upper_y)+
  theme(plot.title = element_text(hjust = 0.5,size = 10))

ggsave("subset_square_plot.pdf",subset_square_plot,width = 12,height = 4)


# now I have to change pop order and values responsible for them accordingly

# my_data_rearranged<-my_data[order(match(my_data$p1,sample_list)),]

# #make dummy list
# new_pop_2<-rep(1: nrow(my_data))
# new_pop_1<-rep(1: nrow(my_data))
#
# #not changing pop3 anyway
# new_pop_3<-my_data$p3
#
# #not changing pop4 anyway
# new_pop_4<-my_data$p4
#
#
# #this does not change
# #new_p1_eq_p2<-cal_data$Pop_A_equals_Pop_B
#
# #these do change
# new_shared_with_cal<-rep(1: nrow(my_data))
# new_shared_with_mello<-rep(1: nrow(my_data))
#
#
# # now I have to change pop order and values responsible for them accordingly
#
# for (sample in 1:length(sample_list_reversed)) {
#   print(sample_list_reversed[sample])
#
#   for (location in 1:nrow(my_data)) {
#     if ((my_data$p1[location])==sample_list_reversed[sample]) {
#       new_pop_1[location]<-sample_list_reversed[sample]
#       new_pop_2[location]<-my_data$p2[location]
#       new_shared_with_cal[location]<-my_data$shared_with_cal[location]
#       new_shared_with_mello[location]<-my_data$shared_with_mello[location]
#     }
#      else if ((my_data$p2[location])==sample_list_reversed[sample]) {
#        new_pop_1[location]<-sample_list_reversed[sample]
#        new_pop_2[location]<-my_data$p1[location]
#        new_shared_with_cal[location]<-my_data$shared_with_cal[location]
#        new_shared_with_mello[location]<-my_data$shared_with_mello[location]
#      }
#   }
#
# }

# my_data_rearranged<-my_data_rearranged[order(match(my_data_rearranged$p1,sample_list)),]
# my_data_rearranged<-my_data_rearranged[order(match(my_data_rearranged$p2,sample_list)),]
# 
# #saving a copy to proof read
# 
# write.table(my_data_rearranged,file = './my_data_arranged_df.tsv', sep='\t', col.names = NA)
# 
# #melt df
# my_data_to_plot <- my_data_rearranged %>%
#   pivot_longer(
#     cols = `shared_with_cal`:`shared_with_mello`,
#     names_to = "Comparison",
#     values_to = "value"
#   )
# 
# #use sample list order as levels
# 
# my_data_to_plot$p1=factor(my_data_to_plot$p1,levels = sample_list)
# my_data_to_plot$p2=factor(my_data_to_plot$p2,levels = sample_list)
# 
# my_plot<-ggplot(my_data_to_plot, aes(x=p1,y = value/1000000,fill=Comparison)) +
#   geom_bar(position="dodge",stat='identity')+
#   facet_wrap(~ p2,nrow = 19)+
#   theme_classic()+
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
#   xlab("Individual B") + ylab("No. of bp (millions)")+
#   scale_fill_manual(values=c( "red","purple"))+
#   ylim(0,upper_y)+
#   theme(plot.title = element_text(hjust = 0.5,size = 10))
# 
# ggsave("my_plot.pdf",my_plot,width = 12,height = 11)
# 
# # ********** remove half of the data for triangle plot *************
# 
# my_data_for_triangle<-my_data_to_plot
# 
# 
# i<-1
# while (i < nrow(my_data_for_triangle)) {
#   checking_row_col_1<-my_data_for_triangle$p1[i]
#   checking_row_col_2<-my_data_for_triangle$p2[i]
# 
#   j<-1
#   while (j < nrow(my_data_for_triangle)) {
#     print(j)
#     j<-j+1
#     if (checking_row_col_1==my_data_for_triangle$p2[j] & checking_row_col_2==my_data_for_triangle$p1[j]) {
#       #list_to_remove <- append(list_to_remove, j)
#       my_data_for_triangle <- my_data_for_triangle[-c(j), ]
#       break
#     }
#   }
#   i<-i+1
# }
# 
# # extract rest of the data form my_data dataframe seperately
# 
# install.packages("sqldf")
# require(sqldf)
# 
# my_data_for_triangle_2<-sqldf('SELECT * FROM my_data_to_plot EXCEPT SELECT * FROM my_data_for_triangle')
# 
# #now I have to change p1 and p2 so I can merge these dataframes
# my_data_for_triangle_2_renamed<-`colnames<-`(my_data_for_triangle_2,c('p2','p1','p3','p4','Comparison','value'))
# 
# #and merge them
# all_data_together<-merge(my_data_for_triangle, my_data_for_triangle_2_renamed, by=c("p1","p2","p3","p4","Comparison"))
# 
# #melt df
# my_triangle_data_to_plot <- all_data_together %>%
#   pivot_longer(
#     cols = `value.x`:`value.y`,
#     names_to = "Comparison2",
#     values_to = "value2"
#   )
# 
# # Combine identifiers to make it easier to compare - pop1***
# 
# my_triangle_data_to_plot$total_comparison <- paste(my_triangle_data_to_plot$Comparison,my_triangle_data_to_plot$Comparison2)
# 
# 
# my_plot_triangle<-ggplot(my_data_for_triangle, aes(x=p1,y = value/1000000,fill=Comparison)) +
#   geom_bar(position="dodge",stat='identity')+
#   facet_wrap(~ p2,nrow = 19)+
#   theme_classic()+
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
#   xlab("Individual B") + ylab("No. of bp (millions)")+
#   scale_fill_manual(values=c( "red","purple"))+
#   ylim(0,upper_y)+
#   theme(plot.title = element_text(hjust = 0.5,size = 10))
# 
# ggsave("my_plot_triangle.pdf",my_plot_triangle,width = 12,height = 11)
# 
# my_plot_final<-ggplot(my_triangle_data_to_plot, aes(x=p1,y = value2/1000000,fill=total_comparison)) +
#   geom_bar(stat = "identity", position = "dodge")+
#   facet_wrap(~ p2,nrow = 19)+
#   theme_classic()+
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
#   xlab("Individual B") + ylab("No. of bp (millions)")+
#   scale_fill_manual(values=c( "red","orange",'purple','blue'))+
#   ylim(0,upper_y)+
#   theme(plot.title = element_text(hjust = 0.5,size = 10))
# 
# ggsave("my_plot_final.pdf",my_plot_final,width = 12,height = 11)



```








