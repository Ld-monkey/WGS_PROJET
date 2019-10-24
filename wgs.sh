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


