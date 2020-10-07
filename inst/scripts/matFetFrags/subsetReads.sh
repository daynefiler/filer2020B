#!/bin/bash
picard FilterSamReads FILTER=includeReadList I=/pathToMyDir/cffDNA/prelimData/keep/FES-0034-4.keep.bam O=/pathToMyDir/cffDNA/prelimData/matFetFrags/maternalReads.bam READ_LIST_FILE=maternalReads.txt

picard FilterSamReads FILTER=includeReadList I=/pathToMyDir/cffDNA/prelimData/keep/FES-0034-4.keep.bam O=/pathToMyDir/cffDNA/prelimData/matFetFrags/fetalReads.bam READ_LIST_FILE=fetalReads.txt

