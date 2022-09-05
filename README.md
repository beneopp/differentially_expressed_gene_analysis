# differentially_expressed_gene_analysis

# Introduction

The following code was implemented based on an assignment for a Computational Genomics course at Carnegie Mellon University. This assignment write-up can be found in hw_2.pdf. The code in the code folder demonstrates my skills in normalization (see normalization.py), differential gene expression analysis (see dispersion.py and de_genes.py), and cell type differentiation (see classification.py).

# Code Outline

normalization.py - I normalized gene expression data using the following methods: RPKM, TPM, and size factor normalizaqtion. The code also shows side-by-side boxplots for every normalization method and the raw count data. These figures display the count data for each sample to determine how well the methods normalize the data.

dispersion.py - I generated scatter plots to assist in identifying genes that were significantly differentially expressed.

de_genes.py - I implemented the Benjamini-Hochberg correction on p-values.

classification.py - I implemented Principal Component Analysis (PCA) to analyze the difference in gene expression profile among the different cell types.
