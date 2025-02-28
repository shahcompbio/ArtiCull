"""
generate_labels.py

This module provides functions to generate labels for input data and save the results to a specified output directory.
It includes functionality for validating input arguments, processing data frames, applying assignment functions, and creating plots.

Functions:
    generate_labels(input_file, output_dir, clone1, clone2, alpha):
        Generates labels for the input data and saves the results to the specified output directory.

    shared_subclonal(mut, alpha):
        Determines if a mutation is shared subclonal based on binomial tests.

    shared_clonal(mut, alpha):
        Determines if a mutation is clonal in two different samples based on binomial tests.

    get_assignment(mut, alpha):
        Determines the assignment of a mutation based on copy numbers and clonality.

    process_df(input_file, clone1, clone2):
        Processes the input file to generate a DataFrame with specific columns for two clones.

    create_plot(df, clone1, clone2, output_dir):
        Creates and saves a scatter plot of CCF (Cancer Cell Fraction) values for two clones.

    validate_arguments(input_file, output_dir):
        Validates the input file and output directory for existence and permissions.
"""

from matplotlib import pyplot # type: ignore
import seaborn as sns # type: ignore
import pandas as pd # type: ignore
import os
import scipy.stats # type: ignore

def generate_labels(input_file, output_dir, clone1, clone2, alpha):
    """
    Generates labels for the input data and saves the results to the specified output directory.
    Args:
        input_file (str): Path to the input file containing the data.
        output_dir (str): Directory where the output files will be saved.
        clone1 (str): Name of the first clone to be used in processing.
        clone2 (str): Name of the second clone to be used in processing.
        alpha (float): Parameter to be used in the assignment function.
    Returns:
        None
    This function performs the following steps:
    1. Validates the input arguments.
    2. Processes the input data frame using the specified clones.
    3. Applies the assignment function to each row of the data frame in parallel.
    4. Saves the resulting data frame with assignments to a TSV file in the output directory.
    5. Prints the count of each assignment group.
    6. Creates and saves a plot based on the processed data.
    """


    validate_arguments(input_file, output_dir)
    df = process_df(input_file, clone1, clone2)
    df['assignment'] = df.parallel_apply(lambda x: get_assignment(x, alpha), axis=1)
    df.to_csv("{}/assignments.tsv".format(output_dir), sep='\t', index=False)
    print(df.groupby('assignment').count())
    create_plot(df, clone1, clone2, output_dir)

def shared_subclonal(mut, alpha):
    """
    Determine if a mutation is shared subclonal based on binomial tests.
    Parameters:
    mut (dict): A dictionary containing mutation data with the following keys:
        - 'var_1': Variant allele count in sample 1.
        - 'tot_1': Total allele count in sample 1.
        - 'var_2': Variant allele count in sample 2.
        - 'tot_2': Total allele count in sample 2.
        - 'cn_1': Copy number in sample 1.
        - 'cn_2': Copy number in sample 2.
    alpha (float): Significance level for the binomial tests.
    Returns:
    bool: True if the mutation is shared subclonal, False otherwise.
    """

    v1, t1 = mut['var_1'], mut['tot_1']
    v2, t2 = mut['var_2'], mut['tot_2']
    c1, c2 = mut['cn_1'], mut['cn_2']

    x = 0.8
    x_adj = x/c1 - (2*x*eps)/c1 + eps
    clonal1 = scipy.stats.binomtest(v1, n=t1, p=x_adj, alternative='less').pvalue
    x = 0.8
    x_adj = x/c2 - (2*x*eps)/c2 + eps
    clonal2 = scipy.stats.binomtest(v2, n=t2, p=x_adj, alternative='less').pvalue

    private1 = scipy.stats.binomtest(v1, n=t1, p=eps, alternative='greater').pvalue

    private2 = scipy.stats.binomtest(v2, n=t2, p=eps, alternative='greater').pvalue

    return all([p < alpha for p in [clonal1, clonal2, private1, private2]]) #and mut['ccf_1'] < 0.75 and mut['ccf_2'] < 0.75

def shared_clonal(mut, alpha):
    """
    Determines if a mutation is clonal in two different samples based on binomial tests.
    Args:
        mut (dict): A dictionary containing mutation data with the following keys:
            - 'var_1' (int): Number of variant reads in sample 1.
            - 'tot_1' (int): Total number of reads in sample 1.
            - 'var_2' (int): Number of variant reads in sample 2.
            - 'tot_2' (int): Total number of reads in sample 2.
            - 'cn_1' (int): Copy number in sample 1.
            - 'cn_2' (int): Copy number in sample 2.
        alpha (float): Significance level for the binomial tests.
    Returns:
        bool: True if the mutation is clonal in both samples (p-value < alpha), False otherwise.
    """

    v1, t1 = mut['var_1'], mut['tot_1']
    v2, t2 = mut['var_2'], mut['tot_2']
    c1, c2 = mut['cn_1'], mut['cn_2']
    x = 0.6
    x_adj = x/c1 - (2*x*eps)/c1 + eps
    clonal1 = scipy.stats.binomtest(v1, n=t1, p=x_adj, alternative='greater').pvalue
    x = 0.6
    x_adj = x/c2 - (2*x*eps)/c2 + eps
    clonal2 = scipy.stats.binomtest(v2, n=t2, p=x_adj, alternative='greater').pvalue

    return (clonal1 < alpha and clonal2 < alpha)

def get_assignment(mut, alpha):
    """
    Determines the assignment of a mutation based on copy numbers and clonality.
    Args:
        mut (dict): A dictionary containing mutation information with keys 'cn_1' and 'cn_2' representing copy numbers.
        alpha (float): A threshold value used to determine clonality.
    Returns:
        int: Returns -1 if either copy number is 0 or if the mutation does not meet clonality criteria.
            Returns 0 if the mutation is shared subclonal.
            Returns 1 if the mutation is shared clonal.
    """

    c1, c2 = mut['cn_1'], mut['cn_2']
    if c1 == 0 or c2 == 0: return -1
    if shared_subclonal(mut, alpha): return 0
    elif shared_clonal(mut, alpha): return 1
    else: return -1



def process_df(input_file, clone1, clone2):
    """
    Processes the input file to generate a DataFrame with specific columns for two clones.
    Args:
        input_file (str): Path to the input file in tab-separated format.
        clone1 (int): Identifier for the first clone.
        clone2 (int): Identifier for the second clone.
    Returns:
        pd.DataFrame: A DataFrame containing the processed data with columns for chromosome,
                    position, and specified values for the two clones. Rows with missing
                    values are dropped.
    """

    input = pd.read_table(input_file)
    df = pd.DataFrame()
    df['chrm'] = input['chrm'].astype('str')
    df['pos'] = input['pos']

    for val in ['var', 'tot', 'cn', 'ccf']:
        df['{}_1'.format(val)] = input['{}_{}'.format(val, clone1)]
        df['{}_2'.format(val)] = input['{}_{}'.format(val, clone2)]

    df = df.dropna()
    return df

def create_plot(df, clone1, clone2, output_dir):
    """
    Creates and saves a scatter plot of CCF (Cancer Cell Fraction) values for two clones.
    Parameters:
    df (pandas.DataFrame): DataFrame containing the data to be plotted. It must have columns 'assignment', 'ccf_1', and 'ccf_2'.
    clone1 (str): Name or identifier for the first clone.
    clone2 (str): Name or identifier for the second clone.
    output_dir (str): Directory where the plot will be saved.
    The function sorts the DataFrame by the 'assignment' column and creates a scatter plot with different colors for different assignments.
    The plot is saved as 'labeled_ccf.pdf' in the specified output directory.
    """

    pyplot.gcf().set_size_inches(12,12)
    df = df.sort_values(by = 'assignment')
    c = df['assignment']
    col = ['r', 'g', 'b']
    for d,l in {-1:'?', 0:'Artifact', 1:'Real'}.items():
        plot_df = df[df['assignment'] == d]
        pyplot.scatter(plot_df['ccf_1'],plot_df['ccf_2'], alpha = 0.2, s = 2, label=l)

    pyplot.gcf().set_size_inches(5,5)
    sns.despine()
    pyplot.legend()
    pyplot.xlabel('CCF Clone {}'.format(clone1))
    pyplot.ylabel('CCF Clone {}'.format(clone2))
    pyplot.savefig("{}/labeled_ccf.pdf".format(output_dir), bbox_inches='tight')


def validate_arguments(input_file, output_dir):
    """
    Validates the input file and output directory for existence and permissions.
    Args:
        input_file (str): Path to the input file that needs to be validated.
        output_dir (str): Path to the output directory that needs to be validated.
    Raises:
        AssertionError: If the input file does not exist or cannot be read.
        AssertionError: If the output directory does not exist or cannot be written to.
    """

    assert os.path.isfile(input_file), f"Input file {input_file} does not exist."
    assert os.access(input_file, os.R_OK), f"Input file exists, but cannot be read due to permissions: {input_file}"
    assert os.path.isdir(output_dir), f"Output directory {output_dir} does not exist."
    assert os.access(output_dir, os.W_OK), f"Output directory exists, but cannot be written to due to permissions: {output_dir}"