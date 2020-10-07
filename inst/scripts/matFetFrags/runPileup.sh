#!/bin/bash
samtools mpileup --positions informativeSites.bed -d 8000 --output-QNAME -f /pathToMyDir/Refs/GRCh38Primary.fna -o informativeReads.pileup /pathToMyDir/cffDNA/prelimData/keep/FES-0034-4.keep.bam 
