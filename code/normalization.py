# normalization.py
# HW2, Computational Genomics, Spring 2022
# andrewid: 

# WARNING: Do not change the file name; Autograder expects it.

# to implement, execute the following command: python normalization.py ../text_file/ReadCounts.txt ../text_file/GeneLengths.txt

import sys
import numpy as np
import matplotlib.pyplot as plt

PER_MILLION = 1/1000000
PER_KILOBASE = 1/1000


# Do not change this function signature
def rpkm(raw_counts, gene_lengths):
    """Find the normalized counts for raw_counts
    Returns: a matrix of same size as raw_counts
    """
    gene_lengths = np.reshape(gene_lengths, (-1, 1))
    reads_per_kilobase = raw_counts/(gene_lengths * PER_KILOBASE)
    million_reads_mapped = np.sum(raw_counts, axis=0).reshape((1, -1)) * PER_MILLION
    return reads_per_kilobase / million_reads_mapped


# Do not change this function signature
def tpm(raw_counts, gene_lengths):
    """Find the normalized counts for raw_counts
    Returns: a matrix of same size as raw_counts
    """
    gene_lengths = np.reshape(gene_lengths, (-1, 1))
    reads_per_kilobase = raw_counts / (gene_lengths * PER_KILOBASE)
    rpk_sum = np.sum(reads_per_kilobase, axis=0).reshape((1, -1)) * PER_MILLION
    return reads_per_kilobase / rpk_sum

   
# define any helper function here


# Do not change this function signature


def size_factor(raw_counts):
    """Find the normalized counts for raw_counts
    Returns: a matrix of same size as raw_counts
    """
    size_factor = compute_size_factor(raw_counts)
    size_factor = np.reshape(size_factor, (1, -1))
    return raw_counts / size_factor


def compute_size_factor(raw_counts):
    geometric_mean = np.product(raw_counts, axis=1) ** (1/raw_counts.shape[1])
    geometric_mean = np.reshape(geometric_mean, (-1, 1))
    ratios = raw_counts / geometric_mean
    return np.median(ratios, axis=0)
    

if __name__=="__main__":
    raw_counts=np.loadtxt(sys.argv[1])
    gene_lengths=np.loadtxt(sys.argv[2])
    
    rpkm1=rpkm(raw_counts, gene_lengths)
    tpm1=tpm(raw_counts, gene_lengths)
    size_factor1=size_factor(raw_counts)
    np.savetxt("output/size_factor_normalized_counts.txt", size_factor1)

    # TODO: write plotting code here
    data_dict = {"Raw Counts": ["raw_counts.png", raw_counts], "RPKM": ["rkpm.png", rpkm1], "TPM": ["tpm.png", tpm1],
                 "Size Factor": ["size_factor.png", size_factor1]}

    for key, value in data_dict.items():

        f = "output/" + value[0]
        data = value[1]

        log_data = np.log2(data)

        fig, ax = plt.subplots()
        ax.boxplot(log_data)

        if key == "Raw Counts":
            plt.title("Boxplot of " + key)
        else:
            plt.title("Boxplot of " + key + " Normalized Data")

        plt.xlabel("Sample")
        plt.ylabel("Log2 Counts")

        plt.savefig(f)
        plt.clf()
