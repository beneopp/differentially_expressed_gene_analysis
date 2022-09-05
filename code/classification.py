# classification.py
# HW2, Computational Genomics, Spring 2022
# andrewid: 

# WARNING: Do not change the file name; Autograder expects it.

# to implement, execute the following command: python normalization.py ../text_file/ReadCounts.txt ../text_file/GeneLengths.txt

import sys

import numpy as np
from scipy.sparse import csc_matrix, save_npz, load_npz

from matplotlib import pyplot as plt
import pandas as pd
import seaborn as sns

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.svm import SVC

# import torch
# from torch import nn
# from torch.utils.data import Dataset, DataLoader


def get_top_gene_filter(data, n_keep = 2000):
    """Select top n_keep most dispersed genes.

    Args:
        data (n x m matrix): input gene expression data of shape num_cells x num_genes
        n_keep (int): number of genes to be kepted after filtration; default 2000

    Returns:
        filter (array of length n_keep): an array of column indices that can be used as an
            index to keep only certain genes in data. Each element of filter is the column
            index of a highly-dispersed gene in data.
    """

    dispersions = calculate_dispersion(data)
    assert len(dispersions) == data.shape[1]
    genes_to_keep = np.argsort(dispersions)[-n_keep:]
    return genes_to_keep


def calculate_dispersion(counts):
    stds = np.var(counts, axis=0)
    means = np.mean(counts, axis=0)
    dispersions = stds / means
    dispersions[np.isnan(dispersions)] = 0
    return dispersions


def reduce_dimensionality_pca(filtered_train_gene_expression, filtered_test_gene_expression, n_components = 20):
    """Train a PCA model and use it to reduce the training and testing data.
    
    Args:
        filtered_train_gene_expression (n_train x num_top_genes matrix): input filtered training expression data 
        filtered_test_gene_expression (n_test x num_top_genes matrix): input filtered test expression data 
        
    Return:
        (reduced_train_data, reduced_test_data): a tuple of
            1. The filtered training data transformed to the PC space.
            2. The filtered test data transformed to the PC space.
    """
    filtered_test_gene_expression = StandardScaler().fit_transform(filtered_test_gene_expression )
    filtered_train_gene_expression = StandardScaler().fit_transform(filtered_train_gene_expression)
    combined_array = np.append(filtered_train_gene_expression, filtered_test_gene_expression, axis=0)

    pca = PCA(n_components=n_components)
    pca.fit(combined_array)

    reduced_train_data = pca.transform(filtered_train_gene_expression)
    reduced_test_data = pca.transform(filtered_test_gene_expression)

    return reduced_train_data, reduced_test_data


def plot_transformed_cells(reduced_train_data, train_labels):
    """Plot the PCA-reduced training data using just the first 2 principal components.
    
    Args:
        reduced_train_data (n_train x num_components matrix): reduced training expression data
        train_labels (array of length n_train): array of cell type labels for training data
        
    Return:
        None

    """
    cell_types = np.unique(train_labels)
    for cell_type in cell_types:
        data_in_cell_type = reduced_train_data[train_labels == cell_type]
        plt.scatter(data_in_cell_type[:, 0], data_in_cell_type[:, 1], label=cell_type)

    plt.xlabel("PC 1")
    plt.ylabel("PC 2")
    plt.title("Plot Transformed Cells")
    plt.legend()
    plt.savefig("output/plot_trans.png")


def train_and_evaluate_svm_classifier(reduced_train_data, reduced_test_data, train_labels, test_labels):
    """Train and evaluate a simple SVM-based classification pipeline.
    
    Before passing the data to the SVM module, this function scales the data such that the mean
    is 0 and the variance is 1.
    
    Args:
        reduced_train_data (n_train x num_components matrix): reduced training expression data
        train_labels (array of length n_train): array of cell type labels for training data
        
    Return:
        (classifier, score): a tuple consisting of
            1. classifier: the trained classifier
            2. The score (accuracy) of the classifier on the test data.

    """

    clf = make_pipeline(StandardScaler(), SVC())
    clf.fit(reduced_train_data, train_labels)

    pred_vals = clf.predict(reduced_test_data)
    score = np.sum(pred_vals == test_labels) / len(test_labels)

    return clf, score


if __name__ == "__main__":
    train_gene_expression = load_npz(sys.argv[1]).toarray()
    test_gene_expression = load_npz(sys.argv[2]).toarray()
    train_labels = np.load(sys.argv[3])
    test_labels = np.load(sys.argv[4])
    
    top_gene_filter = get_top_gene_filter(train_gene_expression)
    filtered_test_gene_expression = test_gene_expression[:, top_gene_filter]
    filtered_train_gene_expression = train_gene_expression[:, top_gene_filter]

    mode = sys.argv[5]
    if mode == "svm_pipeline":
        reduced_train_data, reduced_test_data = reduce_dimensionality_pca(filtered_train_gene_expression,
                                                                          filtered_test_gene_expression)

        np.savetxt("output/data.txt", reduced_train_data)

        plot_transformed_cells(reduced_train_data, train_labels)

        clf, test_score = train_and_evaluate_svm_classifier(reduced_train_data, reduced_test_data,
                                                            train_labels, test_labels)

        pred_train_vals = clf.predict(reduced_train_data)

        train_score = np.sum(pred_train_vals == train_labels) / len(train_labels)

        print("train accuracy:", train_score)
        print("test accuracy:", test_score)

