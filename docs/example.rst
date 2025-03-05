Example Usage
=============

The following example demonstrates how to run ArtiCull.

Example Files
--------------

- `example/example.maf`: Mutation Annotation Format (MAF) file containing candidate variants
- `example/example.bam`: BAM file with sequencing data
- `example/example.bam.bai`: Index files for BAM

Instructions
------------

Follow setup and activate conda environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Follow the directions in :doc:`installation` to create a conda environment and download required genomic tracks.

.. code-block:: bash

    conda env create -f requirements.yml -n articull-env
    conda activate articull-env
    bash scripts/setup_mappability_track.bash [output_directory]

If you have already done this previously, ensure that the `articull-env`` environment is activated.

.. code-block:: bash

    conda activate articull-env

Run classification
^^^^^^^^^^^^^^^^^^^^^

    This will first extract model features and save to `example/example_features.tsv`. Then it will run the pretrained model in `model/preprint_model/` directory to classify and save to `example/example_result.tsv`.

    .. code-block:: bash

        maf=example/example.maf
        bam=example/example.bam
        output_prefix=example/example
        model_dir=models/preprint_model/

        python -m articull classify $maf $output_prefix $model_dir $bam --cores 1

Output
------

The output file, `result.tsv`, is saved in the specified output directory. This is a tab-separated values (TSV) file with the following columns:

.. list-table::
   :header-rows: 1

   * - **Column**
     - **Description**
   * - `chrm`
     - Chromosome where the variant is located.
   * - `pos`
     - Position of the variant on the chromosome.
   * - `ref_allele`
     - Reference allele at the given position.
   * - `alt_allele`
     - Alternate allele identified at the given position.
   * - `result`
     - Classification of the variant (`ARTIFACT`, `PASS`, or `SKIP`).
   * - `prob_artifact`
     - Probability that the variant is an artifact (only provided for classified variants).

Example output:

.. list-table::
   :header-rows: 1

   * - **chrm**
     - **pos**
     - **ref_allele**
     - **alt_allele**
     - **result**
     - **prob_artifact**
   * - 1
     - 201206823
     - G
     - A
     - ARTIFACT
     - 0.9729175263160109
   * - 1
     - 203994039
     - C
     - A
     - PASS
     - 0.01078797299659806
   * - 1
     - 201226655
     - T
     - C
     - SKIP
     - 

Classification criteria:
^^^^^^^^^^^^^^^^^^^^^^^^

- Variants are classified as `ARTIFACT` if `prob_artifact` > 0.5.
- Variants are classified as `PASS` otherwise.
- Variants are classified as `SKIP` if no supporting variant reads are found in the BAM file (e.g., due to realignment during variant calling).