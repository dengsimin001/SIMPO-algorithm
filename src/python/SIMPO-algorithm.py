#!/usr/bin/python
# -*- coding: utf-8 -*-
import io
from itertools import islice
from scipy import stats

# cg_value_file comes from the NCBI GEO (Gene Expression Omnibus) database, it contains DNA methylation data from different samples.cg_gene_type_file comes from Illumina HumanMethylation450 BeadChip's platform, the first column is Illumina target ID, the second column is the gene name, and the third column is the gene region of the probe.
cg_value_file = './data/GSE40279_average_beta-5.txt' 
cg_gene_type_file = './data/GPL(3).txt'


# Division of the cg probes. The promoter regions including TSS1500 while the other regions including gene body, 3'UTR, 5'UTR, and 1stExon.
gene_type_class = ['TSS1500']

# We only selected genes with other region-located and promoter region-located probes greater than or equal to five for further calculation.
threshold = 5


# the header titile of cg_value_file
header_title = [];

# Read cg_value_file and save it to a dictionary with the format {" probe ": [v1,v2..]}.
cg_value_dict = {}
with open(cg_value_file,'r') as cg_value_file_data:
    header_line = cg_value_file_data.readline()
    header_title = header_line.strip().split('\t')
    header_title.pop(0)
    for line in islice(cg_value_file_data,0,None):
        txtArray = line.strip().split('\t')
        if txtArray[0] not in cg_value_dict.keys():
            cgkey = txtArray.pop(0)
            cg_value_dict[cgkey] = txtArray
           

# Read cg_gene_type_file and convert the data to a dictionary in the format of {' probe ':{' gene name1' : 'gene region1',' gene name2' : 'gene region2'}...}.
cg_gene_type_dict = {};
with open(cg_gene_type_file,'r') as cg_gene_type_file_data:
    for line in islice(cg_gene_type_file_data,1,None):
        txtArr = line.strip().split('\t')
        if len(txtArr) < 3:
            continue;

        if(txtArr[0] not in cg_value_dict.keys()):
            continue

        gene_type_dict = {}
        if txtArr[0] not in cg_gene_type_dict.keys():
            gene_name_arr = txtArr[1].strip().split(';')
            gene_type_arr = txtArr[2].strip().split(';')
            
            if len(gene_name_arr) == len(gene_type_arr):
                for i in range(len(gene_name_arr)):
                    gene_type_dict[gene_name_arr[i]] = gene_type_arr[i]

                if len(gene_type_dict) != 0:
                    cg_gene_type_dict[txtArr[0]] = gene_type_dict

# Output intermediate file, which contains probe, gene name and gene region.
tmp_cg_gen_type_file = "./result/tmp_cg_gene_type_file.txt";
with open(tmp_cg_gen_type_file,'w') as tmp_file:
    tmp_file.write("probe\tgene_name\tgene_region\n")
    for (key,value) in cg_gene_type_dict.items():
        for (k,v) in value.items():
            txtTemp = key + '\t' + k +'\t' + v +'\n'
            # print(txtTemp)
            tmp_file.write(txtTemp)


# Convert the format of 'probe-gene name-gene region' into 'gene name-promoter region-probe' or 'gene name-other regions-probe' format.
gene_type_cg_dict = {}
for (key,value) in cg_gene_type_dict.items():

    # key: probe
    # k: gene name
    # v: gene region
    for (k,v) in value.items():
        if  k not in gene_type_cg_dict.keys():
            tempDict = {}
            tempDict.setdefault('class1',set())
            tempDict.setdefault('class2',set())
            if v in gene_type_class:
                tempDict['class1'].add(key)
            else:
                tempDict['class2'].add(key)
            gene_type_cg_dict[k] = tempDict
        else:
            if v in gene_type_class:
                gene_type_cg_dict[k]['class1'].add(key)
            else:
                gene_type_cg_dict[k]['class2'].add(key)

# Limit the number probes for each gene according to the threshold above.
gene_type_cg_limit_dict = {};
for (key,value) in gene_type_cg_dict.items():
    if(len(value['class1']) < threshold or len(value['class2']) < threshold):
        continue
    gene_type_cg_limit_dict[key] = value

# Output intermediate file, including the correspondence between genes and promoter regions or the other regions
tmp_gene_type_cg_file = "./result/tmp_gene_type_cg_file.txt"
with open(tmp_gene_type_cg_file,'w') as tmp_file:
    tmp_file.write('gene_name\tclass_type\tcount\tprobes\n')
    for (key,value) in gene_type_cg_limit_dict.items():
        for (k,v) in value.items():
            tmp_file.write(key+'\t'+k +'\t'+str(len(v))+'\t'+str(v)+'\n')

# T-test for probes in the promoter region of the gene and other regions
SIMPO_algorithm_result_dict = {}
for (k,value) in gene_type_cg_limit_dict.items():

    gene_class1_value_list = []
    gene_class2_value_list = []

    for cgitem in value['class1']:
        gene_class1_value_list.append(cg_value_dict[cgitem])

    for cgitem in value['class2']:
        gene_class2_value_list.append(cg_value_dict[cgitem])
    
    class1_value_list = list(zip(*gene_class1_value_list))
    class2_value_list = list(zip(*gene_class2_value_list))

    # T-test on the sample of the dataset
    ttestArr = []
    for i in range(len(class1_value_list)):
        t1 = list(class1_value_list[i])
        t1 = [ float(x) for x in t1 ]
        t2 = list(class2_value_list[i])
        t2 = [ float(x) for x in t2 ]
        ttest = stats.ttest_ind(t1, t2, equal_var = False)
        ttestArr.append(ttest.statistic)
    SIMPO_algorithm_result_dict[k] = ttestArr

# Output the SIMPO score of each gene in all samples
SIMPO_algorithm_result_file = './result/SIMPO_algorithm_result.txt'
with open(SIMPO_algorithm_result_file,'w') as resultfile:    
    headerTxt = 'gene_name'
    for title in header_title:
        headerTxt = headerTxt +'\t' +  title 
    headerTxt = headerTxt + '\n'
    resultfile.write(headerTxt)
    for (key,value) in SIMPO_algorithm_result_dict.items():
        txtTemp = key 
        for i in range(len(value)):
            txtTemp = txtTemp + '\t'  + str(value[i])
        txtTemp += '\n'
        resultfile.write(txtTemp)

    

