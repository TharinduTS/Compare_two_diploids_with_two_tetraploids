# Compare_two_diploids_with_two_tetraploids

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





