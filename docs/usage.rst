Command-line Usage
==================

Classifying variants with ArtiCull involves two steps using a pretrained model (available in the `models` directory).

1. Extracting features for a list of variants
2. Classifying the variants using extracted features with a pre-trained model

Extracting Features
-------------------

The `extract_features` command processes input variants and generates features required for classification.

Usage
^^^^^

.. code-block:: bash

    python articull/articull.py extract_features \
            <input_file> <output> <bams> \
            [--resources_dir <path>] [--cores <ncores>]

Arguments
^^^^^^^^^

- **<input_file>** (required): Candidate variants file (MAF or VCF format)
- **<output>** (required): Output path for extracted features (tab-separated format)
- **<bams>** (required): One or more BAM files containing sequencing data
- **-\-resources_dir** (optional): Path to directory containing folder `mappability` with downloaded mappability tracks [See :doc:`installation`] (default: `resources/hg19_mappability.bedGraph`)
- **-\-cores** (optional): Number of CPU cores for parallel processing (default: all available cores)


Classification
--------------

The `classify` command processes extracted features using a pre-trained model.
We recommend using the model provided in `models/preprint_model/` for general use cases.

Usage
^^^^^

.. code-block:: bash

    python articull/articull.py classify \
            <features> <output_dir> <model_dir> \
            [--chunksize <n>] [--cores <ncores>]

Arguments
^^^^^^^^^


- **<features>** (required): Input file containing extracted features
- **<output_dir>** (required): Directory for classification results
- **<model_dir>** (required): Directory containing pre-trained model (`model.pkl`) and scaler (`scaler.pkl`)
- **-\-chunksize** (optional): Rows per worker for parallel processing (default: 5000)
- **-\-cores** (optional): Number of CPU cores for parallel processing (default: all available cores)
