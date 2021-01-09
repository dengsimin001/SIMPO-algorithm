# SIMPO-algorithm
**@author**  *`Yuan Quan`, `Fengji Liang`, `Si-Min Deng`, `Yuexing Zhu`, `Ying Chen`, `Jianghui Xiong`*

The algorithm is related to the research *Mining the Selective Remodeling of DNA Methylation in Promoter Regions to Identify Robust Gene-Level Associations with Phenotype*. In this study, we have proposed the statistical difference of DNA Methylation between Promoter and Other Body Region (SIMPO) algorithm to mining gene-level DNA methylation associations with phenotype.</br>


## Content

```
SIMPO-algorithm
└─src
    └─python
        ├─data
        └─result
```

### Quick Start

#### Python

The project code depends on the `scipy` module to run, use the `pip` command to install:

```shell
pip install scipy
```

you can run it by the following command:

```shell
python ./SIMPO-algorithm.py
```

The result is output to the `./result` directory

```
tmp_cg_gene_type_file.txt
tmp_gene_type_cg_file.txt

# The final output file
SIMPO_algorithm_result.txt
```

##### Code implementation instructions:

1. Input file: `src/python/data/test_GPL.txt`, `src/python/data/test_GSE40279_average_beta.txt`

2. Output file: `src/python/result/tmp_cg_gene_type_file.txt.txt`,  `src/python/result/tmp_gene_type_cg_file.txt`, `src/python/result/SIMPO_algorithm_result.txt`

 The final output file `src/python/result/SIMPO_algorithm_result.txt` contains the SIMPO score of each gene in all samples.








