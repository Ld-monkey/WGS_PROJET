#!/bin/bash

echo "---- WGS ----"

# Decompresser l'échantionage G
ECHG_R1_FILE=fastq/EchG_R1.fastq
ECHG_R2_FILE=fastq/EchG_R2.fastq

# Verification que les fichiers sont bien decompresser.
if [ -f "$ECHG_R1_FILE" ] && [ -f "$ECHG_R2_FILE" ]
then
    echo "Les fichier EchG_R1 et EchG_R2 sont décompressés."
else
    gunzip fastq/EchG*.fastq.gz
fi

# 1 - Identifier et quantifier les bactéries présentes dans votre échantillon.
# Aligner les reads avec all_genome.fasta + end to end avec mode rapide.
# retourner un fichier sam.

./soft/bowtie2 --help

if [ -d "result" ] && [] && ; then
    echo "Le dossier results a été créé."
    echo "Le dossier rebuilt_genome a été créé."
    echo "le dossier sam_output a été crée."
else
    mkdir results
    mkdir results/rebuilt_genome
    mkdir results/sam_output
fi

#./soft/bowtie2 -1 -2 -S results/sam_output/sam_out.sam


