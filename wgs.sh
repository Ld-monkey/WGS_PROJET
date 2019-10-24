#!/bin/bash

echo "---- WGS ----"

# Decompresser l'échantionage G
ECHG_FILE=fastq/EchG_R1.fastq
if test -f "$ECHG_FILE";then
    echo "Les fichier EchG_R1 et EchG_R2 sont décompressés."
else
    gunzip fastq/EchG*.fastq.gz
fi

