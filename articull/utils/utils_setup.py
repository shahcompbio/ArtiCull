"""
utils_setup.py

This module provides utility functions for setting up the computational environment.
It includes functionality for setting the number of cores for OpenBLAS and MKL libraries,
and initializing the pandarallel library for parallel processing with pandas.

Functions:
    setup_ncores(ncores=None):
        Sets up the number of cores to be used by OpenBLAS and MKL libraries.

    setup_pandarallel(mode, ncores=None):
        Initializes the pandarallel library with the specified mode and number of cores.
"""

import os
from pandarallel import pandarallel # type: ignore

def setup_ncores(ncores=None):
    """
    Sets up the number of cores to be used by OpenBLAS and MKL libraries.
    This function sets the environment variables 'OPENBLAS_NUM_THREADS' and 'MKL_NUM_THREADS'
    to control the number of threads used by these libraries. This is primarily done to limit
    the number of cores used by numpy, which otherwise uses as many cores as it wants. 
    Be a good citizen and don't hog all the cores on the cluster.

    Note:
        This function must be called before numpy is imported, or else numpy will use as many
        cores as it wants.

    Args:
        ncores (int, optional): The number of cores to be used. If not provided, the default
                                number of cores will be used.

    Returns:
        None
    """
    import psutil
    threads_per_core = psutil.cpu_count() / psutil.cpu_count(logical=False)

    if ncores:
        os.environ['OPENBLAS_NUM_THREADS'] = str(ncores * threads_per_core)
        os.environ['MKL_NUM_THREADS'] = str(ncores * threads_per_core)

def setup_pandarallel(mode, ncores=None):
    """
    Initializes the pandarallel library with the specified mode and number of cores.

    Args:
        mode (str): The mode in which to run pandarallel. If mode is 'classify', the progress bar will be disabled.
        ncores (int, optional): The number of cores to use for parallel processing. 
                                If not specified, pandarallel will use the default number of cores.

    Note:
        The progress bar is best used when output is running to the command line. 
        When it's being output to a log file, it may result in an excessive number of lines.

    Returns:
        None
    """
    progress_bar = mode != 'classify'
    if ncores:
        pandarallel.initialize(progress_bar=progress_bar, nb_workers=ncores, use_memory_fs=True)
    else:
        pandarallel.initialize(progress_bar=progress_bar, use_memory_fs=False)