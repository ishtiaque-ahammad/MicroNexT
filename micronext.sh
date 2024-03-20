#!/bin/bash

SECONDS=0
mkdir 3.annotated_genome 4.metabolic_pathways 1.major_output_files 5.antibiotic_resistant_genes 6.biosynthetic_gene_clusters 
rm -rf species_check

echo "Welcome to Microbial Next generation sequence analysis Tool (MicroNexT)!"
sleep 5s
read -e -p "Please provide your forward read file name: " R1
read -e -p "Please provide your reverse read file name: " R2
read -e -p "Please provide your genus name: " genus
read -e -p "Please provide your species name: " species
read -e -p "Please provide number of CPU cores: " core
echo "Thanks for providing the neccessary information! Let's start!"
sleep 5s

echo "Starting QC of ${R1} and ${R2}!"
sleep 5s

#Fastp
#Activate the fastp environment
eval "$(conda shell.bash hook)"
conda activate fastp

fastp -i ${R1}  -I ${R2} -o trimmed_paired_${R1}.fastq -O trimmed_paired_${R2}.fastq
mv fastp.html quality_control_summary.html
mv quality_control_summary.html ./1.major_output_files/
rm fastp.json

echo "Starting Genome Assembly!"                                                   
sleep 5s 
#Unicycler
#Activate the unicycler environment
eval "$(conda shell.bash hook)"
conda activate unicycler

unicycler -1 trimmed_paired_${R1}.fastq -2 trimmed_paired_${R2}.fastq -o 2.assembled_genome -t ${core}
rm ./2.assembled_genome/*.log
mv ./2.assembled_genome/assembly.fasta ./2.assembled_genome/assembled_genome.fasta

echo "Genome Assembly Finished!"
sleep 5s

echo "Starting Genome Annotation!"
sleep 5s

#Bbmap
#Activate the bbmap environment
eval "$(conda shell.bash hook)"
conda activate bbmap

reformat.sh in=./2.assembled_genome/assembled_genome.fasta out=./2.assembled_genome/assembled_genome_final.fasta minlength=1000

cp ./2.assembled_genome/assembled_genome_final.fasta ./1.major_output_files/assembled_genome_final.fasta

#Prokka
#Activate the prokka environment
eval "$(conda shell.bash hook)"
conda activate prokka

prokka \
    --outdir ./3.annotated_genome/ \
    --force \
    --addgenes \
    --addmrna \
    --compliant \
    --mincontiglen 200 \
    --rfam \
    --rnammer \
    --cpus ${core} \
    --prefix annotation \
    --genus ${genus} \
    --species ${species} \
    --centre "NIB" \
    ./1.major_output_files/assembled_genome_final.fasta

echo "Genome Annotation Finished!"

rm -f trimmed_paired_${R1} trimmed_unpaired_${R1} trimmed_paired_${R2} trimmed_unpaired_${R2}
cp ./3.annotated_genome/*.faa ./1.major_output_files/proteins.faa
cp ./3.annotated_genome/*.ffn ./1.major_output_files/genes.ffn
cp ./3.annotated_genome/*.fna ./1.major_output_files/annotated_genome.fna
cp ./3.annotated_genome/*.gbk ./1.major_output_files/annotated_genome.gbk
cp ./3.annotated_genome/*.gff ./1.major_output_files/annotated_genome.gff
cp ./3.annotated_genome/*.fsa ./1.major_output_files/annotated_genome.fsa
cp ./3.annotated_genome/*.tbl ./1.major_output_files/annotation_table.tbl
cp ./3.annotated_genome/*.tsv ./1.major_output_files/annotation_table_main.tsv
cp ./3.annotated_genome/*.txt ./1.major_output_files/annotation_summary.txt
rm ./3.annotated_genome/*.log

echo "Starting Pathway Analysis!"
sleep 5s 

#Prodigal
#Activate the prodigal environment
eval "$(conda shell.bash hook)"
conda activate prodigal

prodigal -i ./1.major_output_files/assembled_genome_final.fasta -o ./2.assembled_genome/prodigal_proteins.gff3 -d ./2.assembled_genome/prodigal_proteins.fa -f gff -q


#microbeannotator
#Activate the microbeannotator environment
eval "$(conda shell.bash hook)"
conda activate microbeannotator

microbeannotator -i ./3.annotated_genome/annotation.faa -d /media/bioinfo/Data_3/genome_assembly/software/microbe_annotator/MicrobeAnnotator_DB -o ./4.metabolic_pathways/ -m diamond -p 1 -t ${core} --light
rm ./4.metabolic_pathways/metabolic_summary__barplot.pdf
cp ./4.metabolic_pathways/*.pdf ./1.major_output_files/
cp ./4.metabolic_pathways/*.tab ./1.major_output_files/


echo "Metabolic Pathway Analysis Finished!"
sleep 5s 

echo "Starting Antibiotic Resistant Gene Analysis!"
sleep 5s 

#Abricate
#Activate the abricate environment
eval "$(conda shell.bash hook)" 
conda activate abricate

abricate ./1.major_output_files/assembled_genome_final.fasta > ./5.antibiotic_resistant_genes/antibiotic_resistant_genes.tsv
cp ./5.antibiotic_resistant_genes/antibiotic_resistant_genes.tsv ./1.major_output_files/


echo "Antibiotic Resistant Gene Analysis Finished!"
sleep 5s 

echo "Starting Secondary Metabolite Identification!"
sleep 5s 

#Antismash
#Activate the antismash environment
eval "$(conda shell.bash hook)"
conda activate antismash

antismash --cpus ${core} ./1.major_output_files/assembled_genome_final.fasta --genefinding-gff3 ./2.assembled_genome/prodigal_proteins.gff3 --output-dir ./6.biosynthetic_gene_clusters --cb-knownclusters --tigrfam --pfam2go --html-start-compact --output-basename biosynthetic_gene_clusters -v

mv ./6.biosynthetic_gene_clusters/index.html ./6.biosynthetic_gene_clusters/biosynthetic_gene_clusters.html
cp ./6.biosynthetic_gene_clusters/biosynthetic_gene_clusters.html ./1.major_output_files/

#echo "Secondary Metabolite Identification Finished!"
#sleep 5s 

echo "Starting CRISPR-Cas site Identification!"
sleep 5s 

#Cctyper
#Activate the cctyper environment
eval "$(conda shell.bash hook)"
conda activate cctyper

cctyper ./1.major_output_files/assembled_genome_final.fasta ./7.crispr_sites -t ${core}
mv ./7.crispr_sites/plot.png ./7.crispr_sites/crispr_cas_plot.png
mv ./7.crispr_sites/spacers/*.fa ./7.crispr_sites/spacers/crispr_spacers.fa
cp ./7.crispr_sites/crisprs_all.tab ./1.major_output_files/
cp ./7.crispr_sites/crispr_cas_plot.png ./1.major_output_files/
cp ./7.crispr_sites/spacers/crispr_spacers.fa ./1.major_output_files/

echo "CRISPR-Cas site Identification Finished!"
sleep 5s 

duration=$SECONDS

echo "Thanks for using MicroNexT! It took $(($duration / 60)) minutes and $(($duration % 60)) seconds to run!"
