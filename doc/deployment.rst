=====================================
Pipeline configuration and deployment
=====================================

Pipeline structure
------------------

The pipeline currently consists of six independent modules. It also includes a
conda environment file with all the necessary dependencies, as well as the
reference material used in our SFARI project.

:: 

  sv-pipeline
  ├── 00_preprocessing/
  ├── 01_algorithm_integration/
  ├── 02_evidence_assessment/
  ├── 03_variant_filtering/
  ├── 04_variant_resolution/
  ├── 05_annotation/
  ├── README.md
  ├── environment.yaml
  ├── ref/
  └── scripts/
      └── sub_snake.sh

Each module is provided as a Snakemake workflow, and each directory contains
the corresponding Snakefile, configuration files, rules used by the workflow,
and miscellaneous helper scripts.

:: 

  01_algorithm_integration
  ├── README.md
  ├── Snakefile
  ├── cluster.yaml
  ├── config.yaml
  ├── rules
  │   ├── depth_integration.rules
  │   └── pesr_integration.rules
  └── scripts
      ├── make_depth_rdtest_bed.py
      └── make_pesr_rdtest_bed.py


Running the pipeline
--------------------

Automatically deploying these modules in sequence is under development. Until
then, the pipeline requires manual initiation of each module.

.. code-block:: bash

    $ cd $module
    $ vim config.yaml
    $ snakemake

Submitting modules to an LSF cluster
------------------------------------

Modules can optionally be submitted to an LSF cluster through a provided
script that takes advantage of Snakemake's ``--cluster`` functionality.
Support for other scheduling engines through DRMAA is in progress.

.. code-block:: bash

    $ cd $module
    $ vim config.yaml
    $ ../scripts/sub_snake.sh Snakefile

Global configuration
--------------------

Each module is individually configurable through the ``config.yaml`` in the
respective directory.  Several of the configuration variables available here
(e.g. batch key, chromosomes to analyze) should apply globally to all
subworkflows. Implementing this is in progress.

