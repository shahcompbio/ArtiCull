"""
This module provides functions for preprocessing variant and copy number data to 
compute cancer cell fractions (CCFs) for each variant in each clone.
It includes functionality for validating input arguments, processing signal data, 
parsing copy number data, calculating variant counts, computing CCFs, and plotting results.

Functions:
    preprocessing(maf, bam_dirs, signals_dir, output_dir, fullbam, cellclone_file, hscn_file, use_cached_cn):
        Preprocesses variant and copy number data to compute CCFs for each variant in each clone.
"""

import os
from os import path
from collections import Counter
import pandas as pd # type: ignore
pd.options.mode.chained_assignment = None #Suppress SettingWithACopy Warning
import subprocess
from matplotlib import pyplot # type: ignore
import seaborn as sns # type: ignore

from articull._utils.io import get_variants
from articull._utils.bams import match_variants_to_filenames, get_sam, generate_reads

def preprocessing(maf, bam_dirs, signals_dir, output_dir, fullbam, cellclone_file, hscn_file, use_cached_cn, filter_vcf=True):
    """
    Preprocesses variant and copy number data to compute CCFs for each variant in each clone.

    Args:
        maf (str): Path to the MAF (Mutation Annotation Format) file containing variant information.
        bam_dirs (list): List of directories containing BAM files for each sample.
        signals_dir (str): Directory containing signal files for copy number analysis.
        output_dir (str): Directory where the output files will be saved.
        fullbam (bool): Flag indicating whether to use full BAM files.
        cellclone_file (str): Path to the cell clone file.
        hscn_file (str): Path to the HSCN (High-Resolution Copy Number) file.
        use_cached_cn (bool): Flag indicating whether to use cached copy number data.

    Raises:
        RuntimeError: If neither signals_dir nor both cellclone_file and hscn_file are provided.

    Returns:
        None
    """
    
    _validate_arguments(maf, bam_dirs, output_dir, signals_dir, cellclone_file, hscn_file)

    print("1. Reading variants from: {}".format(maf))
    df = get_variants(maf, filter_vcf=filter_vcf)

    print("2. Process copy number information")
    if signals_dir: 
        cellclone_file, hscn_file = _process_signals(df, signals_dir, output_dir, use_cached_cn)
        df, clonemap, clone_ids = _parse_copynumber(df, cellclone_file, hscn_file) 
    elif cellclone_file and hscn_file:
        df, clonemap, clone_ids = _parse_copynumber(df, cellclone_file, hscn_file) 
    else:
        raise RuntimeError('Missing files: Must provide either --signals_dir argument, or --cell_clones AND --hscn')

    print("3. Extracting clone variant counts from: {}".format(bam_dirs))
    df = _get_clone_var_counts(df, bam_dirs, clonemap, clone_ids, fullbam)

    print("4. Computing CCFs")
    df = _compute_ccfs(df, clone_ids)

    result_filename = os.path.join(output_dir, "var_counts.tsv")
    print("5. Outputting result to: {}".format(result_filename))
    df.to_csv(result_filename, sep = '\t', index=False)

    plot_filename = os.path.join(output_dir, "clone_ccfs.png")
    print("6. Creating clone CCF plot at: {}".format(plot_filename))
    _plot_ccfs(df, plot_filename, clone_ids)


def _validate_arguments(maf, bam_dirs, output_dir, signals_dir, cellclone_file, hscn_file):
    """
    Validates the input and output file paths and directories for the preprocessing step.

    Args:
        maf (str): Path to the input MAF file.
        bam_dirs (list of str): List of paths to the input BAM files.
        output_dir (str): Path to the output directory.
        signals_dir (str, optional): Path to the signals directory.
        cellclone_file (str, optional): Path to the cell to clone mapping file.
        hscn_file (str, optional): Path to the signals HSCN file.

    Raises:
        AssertionError: If any of the input files do not exist or cannot be read.
        AssertionError: If the output directory does not exist or cannot be written to.
        AssertionError: If the signals directory, cell to clone mapping file, or signals HSCN file do not exist or cannot be read.
    """
    assert os.path.isfile(maf), f"Input file {maf} does not exist."
    assert os.access(maf, os.R_OK), (
        f"Input file exists, but cannot be read due to permissions: {maf}"
    )
    assert os.path.isdir(output_dir), f"Output directory {output_dir} does not exist."
    assert os.access(output_dir, os.W_OK), (  
        f"Output directory exists, but cannot be written to due to permissions: {output_dir}"
    )
    for bam in bam_dirs:
        assert os.path.isfile(bam), f"Input bam file {bam} does not exist."
        assert os.access(bam, os.R_OK), (
            f"Input bam file exists, but cannot be read due to permissions: {bam}"
        )
    if signals_dir:
        assert os.path.isdir(signals_dir), f"Signals directory {signals_dir} does not exist."
        assert os.access(signals_dir, os.R_OK), (
            f"Signals directory exists, but cannot be read due to permissions: {signals_dir}"
        )
    if cellclone_file:
        assert os.path.isfile(cellclone_file), f"Cell to clone mapping file {cellclone_file} does not exist."
        assert os.access(cellclone_file, os.R_OK), (
            f"Cell to clone mapping file exists, but cannot be read due to permissions: {cellclone_file}"
        )
    if hscn_file:
        assert os.path.isfile(hscn_file), f"Signals HSCN file {hscn_file} does not exist."
        assert os.access(hscn_file, os.R_OK), (
            f"Signals HSCN file exists, but cannot be read due to permissions: {hscn_file}"
        )
    if output_dir:
        assert os.path.isdir(output_dir), f"Output directory {output_dir} does not exist."
        assert os.access(output_dir, os.R_OK), (
            f"Output directory exists, but cannot be read due to permissions: {output_dir}"
        )


def _process_signals(df, signals_dir, output_dir, use_cached):
    """
    Processes signal data and generates a cell-to-clone map.

    Args:
        df (pandas.DataFrame): The dataframe containing the signal data.
        signals_dir (str): The directory where the signal files are located.
        output_dir (str): The directory where the output files will be saved.
        use_cached (bool): If True, use the cached cell-to-clone map if it exists.

    Returns:
        tuple: A tuple containing the path to the cell-to-clone map file and the path to the signals CNS file.
    """
    signals_result = path.join(signals_dir, 'signals.Rdata')
    signals_cns = path.join(signals_dir, 'hscn.csv.gz')
    outfile = os.path.join(output_dir, "signals_clones.tsv")

    if use_cached and os.path.isfile(outfile):
        print("Using existing cell-to-clone map: {}".format(outfile))
        return outfile, signals_cns
    else:
        directory = path.dirname(path.realpath(path.expanduser(__file__)))
        command = "Rscript {}/extract_cell2clone.R {} {}".format(directory, signals_result, outfile)
        subprocess.check_call(command.split(' '))
        return outfile, signals_cns


def _parse_copynumber(df, cellclone_file, hscn_file):
    """
    Parses copy number data and maps it to clones.

    Args:
        df (pd.DataFrame): DataFrame containing the data to be processed.
        cellclone_file (str): Path to the file containing cell to clone mapping.
        hscn_file (str): Path to the file containing copy number data.

    Returns:
        pd.DataFrame: Updated DataFrame with copy number columns added.
        pd.DataFrame: DataFrame containing the cell to clone mapping.
        list: List of unique clone IDs.
    """
    clonemap = pd.read_table(cellclone_file, index_col=0)
    clone_ids = sorted(clonemap['clone_id'].unique())
    ids_map = clonemap['clone_id']

    try:
        cn_df = pd.read_table(hscn_file, usecols=['cell_id', 'chr', 'start', 'end', 'Maj', 'Min'], sep=',')
    except ValueError:
        cn_df = pd.read_table(hscn_file, usecols=['cell_id', 'chr', 'start', 'end', 'A', 'B'], sep=',')
        cn_df['Maj'] = cn_df['A']
        cn_df['Min'] = cn_df['B']

    # Cell to clone map table    
    cn_df['clone']= cn_df['cell_id'].map(ids_map)
    clone_cns = cn_df.groupby(['chr', 'start', 'end', 'clone'])[['Maj', 'Min']].median().reset_index()
    clone_cns['CN'] = clone_cns['Maj'] + clone_cns['Min']

    del cn_df # large file
    clone_cn_dfs = {}

    for chrm in clone_cns['chr'].unique():
        clone_cns1 = clone_cns[clone_cns['chr'] == chrm]
        clone_cns1 = clone_cns1.pivot(index=['start', 'end'], columns = 'clone', values = 'CN')
        clone_cns1 = clone_cns1[sorted(clone_cns1.columns)]
        clone_cns1.index = pd.IntervalIndex.from_tuples(clone_cns1.index, closed='both')
        clone_cn_dfs[chrm] = clone_cns1

    columns = ['cn_{}'.format(c) for c in clone_ids]

    def get_cn(x):
        try:
            result = clone_cn_dfs[x['chrm']].loc[int(x['pos'])].values
            return result

        except KeyError:
            return [float('nan')]*len(columns)

    result = df.parallel_apply(get_cn, axis=1, result_type="expand")
    df[columns] = result

    return df, clonemap, clone_ids

def _get_clone_var_counts(df, data_dirs, clonemap, clone_ids, fullbam=False):
    """
    Calculates variant and total read counts for each variant in each clone in the given DataFrame.

    Args:
        df (pandas.DataFrame): DataFrame containing variant information.
        data_dirs (list of str): List of directories containing SAM/BAM files.
        clonemap (pandas.DataFrame): DataFrame mapping cell barcodes to clone IDs.
        clone_ids (list of str): List of clone IDs to be considered.
        fullbam (bool, optional): If True, use the first directory in data_dirs for all filenames. Defaults to False.

    Returns:
        pandas.DataFrame: DataFrame with variant and total read counts for each clone.
    """
    def get_clone_list(reads):
        cells = [read for read in reads]
        clones= []
        for cell in cells:
            try:
                clones.append(clonemap.loc[cell]['clone_id'])
            except:
                clones.append('0')
        return clones

    def process_variant(x):
        var_sams = [get_sam(data_dir, x['filename']) for data_dir in data_dirs]

        reads = []
        for sam in var_sams: reads+=[r for r in generate_reads(sam, x)]

        read_groups = [r[0].alignment.get_tag('CB:Z:') for r in reads]
        clones = get_clone_list(read_groups)

        var_count = Counter()
        all_count = Counter()

        for c, r in zip(clones, reads):
            if r[1]: var_count[c]+=1
            all_count[c] += 1

        var_data = []
        for c in clone_ids:
            var_data += [var_count[c], all_count[c]]

        return var_data

    if fullbam: df['filename'] = data_dirs[0]
    else:
        df = match_variants_to_filenames(df, data_dirs)

    columns = []
    for c in clone_ids: columns+=['var_{}'.format(c), 'tot_{}'.format(c)]

    df[columns] = df.parallel_apply(process_variant, axis=1, result_type="expand")

    del df['filename']
    return df

def _compute_ccfs(df, clone_ids):
    """
    Computes cancer cell fractions (CCFs) for given clone IDs and adds them to the DataFrame.

    Args:
        df (pandas.DataFrame): The input DataFrame containing the necessary columns for computation.
        clone_ids (list): A list of clone IDs for which CCFs need to be computed.

    Returns:
        pandas.DataFrame: The DataFrame with added CCF columns for each clone ID.

    The function expects the DataFrame to have the following columns for each clone ID:
    - 'var_<clone_id>': Variant allele frequency for the clone.
    - 'cn_<clone_id>': Copy number for the clone.
    - 'tot_<clone_id>': Total copy number for the clone.
    The computed CCF for each clone ID will be added to the DataFrame as 'ccf_<clone_id>'.
    """
    for clone in clone_ids:
        ccf = 'ccf_{}'.format(clone)
        var = 'var_{}'.format(clone)
        cn = 'cn_{}'.format(clone)
        tot = 'tot_{}'.format(clone)

        df[ccf] = df[var] *  df[cn] / df[tot]

    return df

def _plot_ccfs(df, output, clone_ids):
    """
    Plots the cancer cell fractions (CCFs) for given clone IDs using pair plots.

    Args:
        df (pandas.DataFrame): DataFrame containing the CCF data for different clones.
        output (str): The file path where the plot image will be saved.
        clone_ids (list): List of clone IDs to be plotted.

    Returns:
        None
    """
    sns.set_context('paper', font_scale=1.5)
    df_plot = pd.DataFrame()
    for c in clone_ids:
        df_plot[c] = df['ccf_{}'.format(c)]

    g = sns.pairplot(df_plot, kind='scatter', plot_kws={'alpha':0.2, 's':8}, corner=True)
    g.set(xlim=(-0.1,1.1), ylim = (-0.1,1.1))
    pyplot.savefig(output, bbox_inches = 'tight')