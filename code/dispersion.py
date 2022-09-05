import numpy as np
import matplotlib.pyplot as plt


def calculate_dispersion(counts):
    stds = np.std(counts, axis=1)
    means = np.mean(counts, axis=1)
    dispersion = stds / means
    dispersion[np.isnan(dispersion)] = 0
    return means, dispersion


def plot_results(results):
    x_cases = np.log2(results["case"][:, 0])
    y_cases = np.log2(results["case"][:, 1])
    plt.scatter(x_cases, y_cases, color="red", label="Cases")

    x_controls = np.log2(results["control"][:, 0])
    y_controls = np.log2(results["control"][:, 1])
    plt.scatter(x_controls, y_controls, color="blue", label="Controls")

    plt.legend()
    plt.xlabel("Log2 Mean")
    plt.ylabel("Log2 Dispersion")
    plt.title("Dispersion vs Mean")

    plt.savefig("output/mean_vs_disp.png")


def compute_foldchanges(results):
    case_means = results["case"][:, 0]
    control_means = results["control"][:, 0]
    log_foldchanges = np.log2(case_means / control_means)
    return log_foldchanges, case_means


def plot_foldchanges(log_foldchange):
    plt.clf()
    plt.scatter(np.log2(case_means), log_foldchange)
    plt.xlabel("Log2 Case Mean")
    plt.ylabel("Log2 Fold Change")
    plt.title("Log2 Case Means vs. Log2 Fold Changes")
    plt.savefig("output/fold_changes.png")


def find_up_regulated_genes(log_foldchanges, gene_names):
    up_regulated_genes = log_foldchanges[log_foldchanges > 0]
    top_gene_index = np.argsort(up_regulated_genes)[-10:][::-1]
    top_gene_names = gene_names[log_foldchanges > 0][top_gene_index]
    print("top_gene_names")
    print(", ".join(top_gene_names))
    print()


def find_down_regulated_genes(log_foldchanges, gene_names):
    down_regulated_genes = log_foldchanges[log_foldchanges < 0]
    down_regulated_index = np.argsort(down_regulated_genes)[0:10]
    down_gene_names = gene_names[log_foldchanges < 0][down_regulated_index]
    print("down gene names")
    print(", ".join(down_gene_names))


norm_counts = np.loadtxt("output/size_factor_normalized_counts.txt")

gene_names = np.loadtxt("../text_files/GeneNames.txt", dtype=str)

labels = np.loadtxt("../text_files/labels.txt").astype(str)
labels[labels == "1.0"] = "case"
labels[labels == "2.0"] = "control"

results = dict()
num_of_genes = norm_counts.shape[0]

for data_type in ["case", "control"]:

    sample_index = labels == data_type

    data_type_results = np.zeros((num_of_genes, 2))
    data = norm_counts[:, sample_index]

    means, dispersion = calculate_dispersion(data)

    data_type_results[:, 0] = means
    data_type_results[:, 1] = dispersion

    results[data_type] = data_type_results

plot_results(results)
log_foldchanges, case_means = compute_foldchanges(results)
plot_foldchanges(log_foldchanges)
find_up_regulated_genes(log_foldchanges, gene_names)
find_down_regulated_genes(log_foldchanges, gene_names)







