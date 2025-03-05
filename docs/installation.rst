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

Downloading Reference Data Tracks
------------------------------

Use the provided script to download and process the reference genome mappability track. Note that the download requires ~1GB of space and expands to ~5GB when uncompressed.
By default, files are saved to the `resources` directory unless an alternative output directory is specified.

.. code-block:: bash

    bash scripts/setup_mappability_track.bash [optional: output_directory]

Currently only `hg19`/`GRCh37` is supported. Support for additional reference genomes coming soon. For other genome versions, please open an issue on the `ArtiCull GitHub repo <https://github.com/shahcompbio/ArtiCull/issues>`_.