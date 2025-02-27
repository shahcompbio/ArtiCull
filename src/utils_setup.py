import os
from pandarallel import pandarallel

def setup_ncores(ncores=None):
    import psutil
    threads_per_core = psutil.cpu_count() / psutil.cpu_count(logical=False)

    if ncores:
        os.environ['OPENBLAS_NUM_THREADS'] = str(ncores * threads_per_core)
        os.environ['MKL_NUM_THREADS'] = str(ncores * threads_per_core)


def setup_pandarallel(mode, ncores=None):
    progress_bar = mode != 'classify'
    if ncores:
        pandarallel.initialize(progress_bar=progress_bar, nb_workers = ncores, use_memory_fs=True)
    else:
        pandarallel.initialize(progress_bar=progress_bar, use_memory_fs=False)
