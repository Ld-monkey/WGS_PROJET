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

# Verifie que les bons dossiers sont créés.
if [ -d "results" ] && [ -d "results/rebuilt_genome" ] && [ -d "results/sam_output" ]; then
    echo "Le dossier results a été créé."
    echo "Le dossier rebuilt_genome a été créé."
    echo "le dossier sam_output a été crée."
else
    mkdir results
    mkdir results/rebuilt_genome
    mkdir results/sam_output
fi

# 1-a : Indexer la banque à l'aide de bowtie2-build
# Vérifier que le fichier a bien été créé.
if [ -f "results/rebuilt_genome/built_genome.1.bt2" ]; then
    echo "Le génome a été built."
else
    echo "Le génome va etre crée."
    ./soft/bowtie2-build databases/all_genome.fasta results/rebuilt_genome/built_genome
fi

# Alignement
if [ -f "results/sam_output/sam_echG.sam" ]; then
    echo "sam crée."
else
    echo "sam non crée."
    ./soft/bowtie2 -p 6 --end-to-end -x results/rebuilt_genome/built_genome -1 fastq/EchG_R1.fastq -2 fastq/EchG_R2.fastq -S results/sam_output/sam_echG.sam
fi
#
echo "Alignement est fini"
# 2 - Calculer l'abondance de chaque bactérie en analysant le fichier sam.
# 2 - a : Convertir votre fichier sam en bam (format binaire) (samtools view) + compression rapide

# Make samtools si ce n'est pas fait.
make soft/samtools-1.6/

# Vérification du dossier sam_to_bam.
if [ -d "results/sam_to_bam" ]; then
    echo "Le dossier sam_to_bam a déjà été créé."
else
    mkdir results/sam_to_bam
    echo "Le dossier sam_to_bam est créé."
fi

# Création du bam.
if [ -f "results/sam_to_bam/output.bam" ]; then
    echo "bam déjà crée."
else
    echo "bam non crée."
    ./soft/samtools-1.6/samtools view -b --thread 6 -1 results/sam_output/sam_echG.sam -o results/sam_to_bam/output.bam
    echo "Samtool view a créé le binaire."
fi

# 2 - b : Trier votre bam (samtools sort).
# Vérification du dossier sort_bam.
if [ -d "results/sort_bam" ]; then
    echo "Le dossier sort_bam a déjà été créé."
else
    mkdir results/sort_bam
    echo "Le dossier sort_bam est créé."
fi

# Création du bam trié.
if [ -f "results/sort_bam/sort_output.bam" ]; then
    echo "bam trié déjà crée."
else
    echo "bam trié non crée."
    ./soft/samtools-1.6/samtools sort results/sam_to_bam/output.bam -o results/sort_bam/sort_output.bam
    echo "Le bam est trié."
fi

# 2 - c : Indexer le fichier bam.
if [ -f "results/sort_bam/sort_output.bam.bai" ]; then
    echo "bam indexé déjà crée."
else
    echo "bam indexé non crée."
    ./soft/samtools-1.6/samtools index results/sort_bam/sort_output.bam
    echo "Le fichier a été indexé."
fi

# 2 - d : Extraction du comptage (samtools idxstats)

# Verification de l'existence du dossier idxstats.
if [ -d "results/idxstats" ]; then
    echo "Le dossier idxstats a déjà été créé."
else
    mkdir results/idxstats
    echo "Le dossier idxstats est créé."
fi

# Création des statistique avec idxstats.
if [ -f "results/idxstats/stats.tsv" ]; then
    echo "Les statistiques sont déjà crées."
else
    echo "Les statistiques sont en cours."
    ./soft/samtools-1.6/samtools idxstats results/sort_bam/sort_output.bam > results/idxstats/stats.tsv
    echo "Statistique depuis le fichier trié est réalisé."
fi

# 2 - f : L’association gi ->.
mkdir results/association
echo "Création du dossier association."
grep ">" databases/all_genome.fasta|cut -f 2 -d ">" >results/association/association.tsv
echo "Fin de création de association.tsv"

# 3 : Assembler le génome avec megahit avec un k-mer max 21.
# megahit en mode de consommation mémoire minimum.
mkdir results/assembler_genome
./soft/megahit -1 fastq/EchG_R1.fastq  -2 fastq/EchG_R2.fastq --mem-flag 0 --k-max 21 --k-min 21 -o results/assembler_genome/
