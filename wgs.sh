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

# Verification de l'existence du dossier assembler_genome .
if [ -d "results/assembler_genome" ]; then
    echo "Le génome des bactéries sont déjà assemblés avec megahit."
else
    ./soft/megahit -1 fastq/EchG_R1.fastq -2 fastq/EchG_R2.fastq --mem-flag 0 --k-max 21 --k-min 21 -o results/assembler_genome
    echo "Le dossier assembler_genome est créé."
    echo "Fin de l'assemblage avec megahit"
fi

# 4 : Prédire les gènes présents sur vos contigs avec prodigal.
# sortie en format fasta.
mkdir results/prediction_gene

if [ -d "results/prediction_gene" ]; then
    echo "Le dossier prediction_gene est déjà créé."
else
    mkdir results/prediction_gene
fi

if [ -f "results/prediction_gene/pred_file.fasta" ]; then
    echo "La prédiction des gènes existe déjç."
else
    ./soft/prodigal -i results/assembler_genome/final.contigs.fa -d results/prediction_gene/pred_file.fasta
    echo "Fin de la prédiction"
fi

# 5 : Sélectionner les gènes “complets”.
mkdir results/genes_complet
sed "s:>:*\n>:g" results/prediction_gene/pred_file.fasta | sed -n "/partial=00/,/*/p"|grep -v "*" > results/genes_complet/genes_full.fna
echo "Fin de la sélection complète des gènes."

# 6 : Annoter les gènes “complets” contre la banque resfinder (database/resfinder.fna) à l’aide de blastn. Vous sélectionnerez à l’aide de blast, les gènes avec une identité de >=80% et une couverture >= 80% pour une evalue supérieure à 1E-3.
mkdir results/blastn
./soft/blastn -query results/genes_complet/genes_full.fna -db databases/resfinder.fna -evalue 0.001 -perc_identity 80 -qcov_hsp_perc 80 -outfmt "6 qseqid sseqid pident qcovs evalue" -best_hit_score_edge 0.00001 -out results/blastn/EchG_blastn.txt

echo "Fin du programme."
