#!/bin/bash
#
# make_workflow_diagram.sh
#
# 
#
# Copyright (C) 2017 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

outdir=doc/figures
snakemake -F --dag | dot -T png > ${outdir}/dag.png
snakemake -F --rulegraph | dot -T png > ${outdir}/rulegraph.png
