.. YAML documentation format borrowed from RTD

=====================
01. Cohort clustering
=====================

The algorithm integration module combines variant predictions across each set
of algorithms (i.e. PE/SR and RD), then integrates across the two classes of
evidence.

Module configuration and input
==============================

The configuration file (``config.yaml``) outlines the module's inputs and
parameters.

Batches
-------

The workflow requires a ``batch`` key, describing each sample in the cohort and
the batch to which they belong. An example:

:: 

  sample	group	batch
  11002.fa	11002	Pilot
  11002.mo	11002	Pilot
  11002.p1	11002	Pilot
  11002.s1	11002	Pilot
  11006.fa	11006	Phase1
  11006.mo	11006	Phase1
  11006.p1	11006	Phase1
  11006.s1	11006	Phase1

Note the ``group`` column. Each algorithm was run on some group of samples,
which we are now trying to integrate. The ``group`` column indicates the ID
used for each group in the original algorithm runs, and is used to identify
the corresponding VCF for each sample.

The configuration variable ``batches`` indicates the key to use.

.. code-block:: yaml
    
    batches: ref/batch.key

PE/SR algorithms (VCFs)
-----------------------

Calls from PE/SR algorithms must be provided as standardized VCFs (see Module
0). The filepath to each individual VCF must be labeled with the file's
``source``, i.e., the algorithm which produced it, and the ``group`` of samples
included in the file. As example:

:: 

  delly.11002.vcf.gz
  delly.11006.vcf.gz

The workflow will search for these files in the specified ``input_vcfs``
directory.

.. code-block:: yaml
  
    input_vcfs: ../00_preprocessing/filtered_vcfs/

Depth algorithms (BEDs)
-----------------------

Depth calls should be concatenated into one BED per batch and CNV type, with
each entry corresponding to a call in a single sample (see Module 0 for details).
The filepath to each individual BED must again be formatted accordingly, here
with the ``batch`` of samples and the ``svtype`` included. As example:

:: 

  Phase1.DEL.vcf.gz

The workflow will search for these files in the specified ``input_beds``
directory.

.. code-block:: yaml
  
    input_vcfs: ../00_preprocessing/std_beds/


Workflow wildcards
==================

The Snakemake workflow parallelizes algorithm integration by *batch*, by
*source*, and by *chromosome*. 


==========  ===========
Wildcard    Description
==========  ===========
``group``   ID for original VCFs. 
``source``  Source algorithm. (Permitted options defined in configuration.)
``batch``   Batch. All VCFs belonging to a batch are aggregated together. 
            (Defined in configuration or parsed automatically from batch key.)
``chrom``   Chromosome. (1-22, X, or Y)
``svtype``  SV type. (DEL, DUP, INV, or BND)
==========  ===========

Module inputs
=============


Depth algorithms (BEDs)
-----------------------

.. code-block:: yaml
  
    input_beds: ../00_preprocessing/std_beds/

Module parameters
=================

VCFCluster
----------

The ``vcfcluster`` block sets the corresponding parameters of ``svtools vcfcluster``. See ``svtools vcfcluster -h`` for further details on each parameter.

.. code-block:: yaml
  
    vcfcluster:
        dist: 300    
        frac: 0.1
        blacklist: "ref/b37.lumpy.exclude.4-13.bed.gz"
        svsize: 0
        flags: ""

BEDCluster
----------

The ``bedcluster`` block sets the corresponding parameters of ``svtools bedcluster``. See ``svtools bedcluster -h`` for further details on each parameter.

.. code-block:: yaml

    bedcluster:
        frac: 0.8
        flags: "--merge-coordinates"
