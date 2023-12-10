#!/usr/bin/env bash

cat ../tximport/dmel-all-r6.33.gtf | tr " " "\t" | awk 'BEGIN {OFS = '\t'} {if($1 == "X" && $3 == "gene") print $10}' | sed 's/"//g;s/;//g' > dmel_X

cat ../tximport/simV3/Dsim_3.0_genomic.gtf | tr " " "\t" | awk 'BEGIN {OFS = '\t'} {if($1 == "NC_052525.1" && $3 == "gene") print $10}' | sed 's/"//g;s/;//g' > dsim_X

cat ../tximport/dmel-all-r6.33.gtf | tr " " "\t" | awk 'BEGIN {OFS = '\t'} {if($1 == "Y" && $3 == "gene") print $10}' | sed 's/"//g;s/;//g' > dmel_Y
