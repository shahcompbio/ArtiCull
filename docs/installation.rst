Installation
============

Cloning the Repository
----------------------

First, clone the ArtiCull repository from GitHub to your local machine:

.. code-block:: bash

    git clone https://github.com/shahcompbio/ArtiCull.git
    cd ArtiCull

Setting up the Conda Environment
--------------------------------

The recommended way to set up ArtiCull is using conda. Alternatively, you can ensure all packages listed in `requirements.yml` are available.

First, ensure you have `conda <https://docs.conda.io/en/latest/>`_ installed on your machine. Then run:

.. code-block:: bash

    conda env create -f requirements.yml -n articull-env
    conda activate articull-env

Using Docker
------------

A ``Dockerfile`` is provided as an alternative to conda. The image carries only what the
classification path needs, so it is considerably smaller than the full conda environment:
``bedtools`` plus ``numpy``, ``scipy``, ``pandas``, ``pysam``, ``scikit-learn``,
``pandarallel``, ``joblib`` and ``psutil``. Build it from the repository root:

.. code-block:: bash

    docker build -t articull:latest .

The ``models/preprint_model`` weights are baked into the image at
``/opt/articull/models/preprint_model``. The mappability tracks are not — they are ~5GB, so
mount them instead. The image's default ``--resources_dir`` is ``/opt/articull/resources``,
so mounting there means you never have to pass the flag:

.. code-block:: bash

    docker run --rm --shm-size=2g \
        -v "$PWD:/data" \
        -v /path/to/resources:/opt/articull/resources:ro \
        articull:latest classify \
            /data/sample.maf \
            /data/out/sample \
            /opt/articull/models/preprint_model \
            /data/sample.bam \
            --cores 8

``--shm-size`` matters: ``pandarallel`` is initialized with ``use_memory_fs=True`` and Docker
defaults ``/dev/shm`` to 64MB, which is too small for real cohorts.

The image supports ``classify``, ``extract_features``, ``classify_variants`` and
``train_classifier``. The ``train_genlabels`` and ``train_preprocessing`` modes are **not**
included, as they pull in ``matplotlib``/``seaborn`` (and, for ``train_preprocessing``, R via
``scripts/misc/extract_cell2clone.R``). If you need them, add to the image:

.. code-block:: dockerfile

    RUN pip install --no-cache-dir 'matplotlib>=3.7' 'seaborn>=0.12'
    # train_preprocessing with --signals_dir additionally needs R and Rscript on PATH

Note that ``python -m articull`` is run from a source copy on ``PYTHONPATH`` rather than
being ``pip install``\ ed, because ``articull/__init__.py`` currently defines neither a
module docstring nor ``__version__``, both of which the ``flit`` backend in
``pyproject.toml`` declares as dynamic and therefore requires.

Downloading Reference Data Tracks
---------------------------------

Use the provided script to download and process the reference genome mappability track. Note that the download requires ~1GB of space and expands to ~5GB when uncompressed.
By default, files are saved to the `resources` directory unless an alternative output directory is specified.

.. code-block:: bash

    bash scripts/setup_mappability_track.bash [optional: output_directory]

The Docker image ships ``curl`` and ``bigWigToBedGraph`` so the same script can be run inside
a container, writing to a mounted directory:

.. code-block:: bash

    docker run --rm -v /path/to/resources:/resources \
        --entrypoint bash articull:latest \
        /opt/articull/articull/setup_mappability_track.bash /resources

``bigWigToBedGraph`` is only published by UCSC for linux/amd64. On an arm64 host, build the
image with ``--platform linux/amd64`` for this step, or prepare the bedGraphs elsewhere and
mount them.

Currently only `hg19`/`GRCh37` is supported. Support for additional reference genomes coming soon. For other genome versions, please open an issue on the `ArtiCull GitHub repo <https://github.com/shahcompbio/ArtiCull/issues>`_.