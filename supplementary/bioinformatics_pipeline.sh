#!/usr/bin/env bash
# ==============================================================================
# PIPELINE: Metagenomic Profiling of Cattle Gut Microbiome & Resistome
# Project: Humanisation of Cattle Gut Microbiome (Urban Stray vs. Rural Domestic)
# Author: Angsuman Das (PhD Research)
# ==============================================================================

set -euo pipefail

# Configuration & Resource Allocation
THREADS=16
RAW_DIR="00_raw_reads"
QC_DIR="01_clean_reads"
HOST_FILTER_DIR="02_decontaminated_reads"
TAXONOMY_DIR="03_taxonomy"
FUNCTIONAL_DIR="04_functional"
RESISTOME_DIR="05_resistome"
ASSEMBLY_DIR="06_assembly"

# Reference Database Paths (Set appropriate local database paths)
HOST_INDEX="references/bos_indicus_bowtie2/bos_indicus"
KRAKEN2_DB="/databases/kraken2_gtdb_r214"
CARD_DB="/databases/card_data"
HUMANN_NUCLEOTIDE_DB="/databases/chocophlan"
HUMANN_PROTEIN_DB="/databases/uniref90"

mkdir -p "${QC_DIR}" "${HOST_FILTER_DIR}" "${TAXONOMY_DIR}" "${FUNCTIONAL_DIR}" "${RESISTOME_DIR}" "${ASSEMBLY_DIR}"

sample="$1"
R1="${RAW_DIR}/${sample}_R1.fastq.gz"
R2="${RAW_DIR}/${sample}_R2.fastq.gz"

echo "=== [Step 1/6] Quality Trimming with fastp for sample: ${sample} ==="
fastp \
  --in1 "${R1}" --in2 "${R2}" \
  --out1 "${QC_DIR}/${sample}_clean_R1.fastq.gz" \
  --out2 "${QC_DIR}/${sample}_clean_R2.fastq.gz" \
  --detect_adapter_for_pe \
  --qualified_quality_phred 20 \
  --unqualified_percent_limit 20 \
  --length_required 75 \
  --thread "${THREADS}" \
  --html "${QC_DIR}/${sample}_fastp.html" \
  --json "${QC_DIR}/${sample}_fastp.json"

echo "=== [Step 2/6] Host DNA Depletion (Bowtie2 vs. Host Genome) ==="
bowtie2 \
  -p "${THREADS}" \
  -x "${HOST_INDEX}" \
  -1 "${QC_DIR}/${sample}_clean_R1.fastq.gz" \
  -2 "${QC_DIR}/${sample}_clean_R2.fastq.gz" \
  --un-conc-gz "${HOST_FILTER_DIR}/${sample}_microbial_R%.fastq.gz" \
  -S /dev/null

CLEAN_R1="${HOST_FILTER_DIR}/${sample}_microbial_R1.fastq.gz"
CLEAN_R2="${HOST_FILTER_DIR}/${sample}_microbial_R2.fastq.gz"

echo "=== [Step 3/6] Taxonomic Profiling (Kraken2 + Bracken) ==="
kraken2 \
  --db "${KRAKEN2_DB}" \
  --threads "${THREADS}" \
  --paired "${CLEAN_R1}" "${CLEAN_R2}" \
  --output "${TAXONOMY_DIR}/${sample}.kraken" \
  --report "${TAXONOMY_DIR}/${sample}.kreport"

bracken \
  -d "${KRAKEN2_DB}" \
  -i "${TAXONOMY_DIR}/${sample}.kreport" \
  -o "${TAXONOMY_DIR}/${sample}_bracken_species.tsv" \
  -w "${TAXONOMY_DIR}/${sample}_bracken_species.kreport" \
  -r 150 -l S -t 10

echo "=== [Step 4/6] Resistome & Virulome Profiling (RGI / CARD) ==="
# Convert paired reads to interleaved or analyze contigs post-assembly
rgi bwt \
  --read_one "${CLEAN_R1}" \
  --read_two "${CLEAN_R2}" \
  --output_file "${RESISTOME_DIR}/${sample}_rgi_reads" \
  --threads "${THREADS}" \
  --clean

echo "=== [Step 5/6] Functional Pathway Profiling (HUMAnN 3.0) ==="
# Concatenate clean reads for HUMAnN input
cat "${CLEAN_R1}" "${CLEAN_R2}" > "${HOST_FILTER_DIR}/${sample}_merged.fastq.gz"

humann \
  --input "${HOST_FILTER_DIR}/${sample}_merged.fastq.gz" \
  --output "${FUNCTIONAL_DIR}/${sample}_humann" \
  --threads "${THREADS}" \
  --nucleotide-database "${HUMANN_NUCLEOTIDE_DB}" \
  --protein-database "${HUMANN_PROTEIN_DB}" \
  --remove-temp-output

humann_renorm_table \
  --input "${FUNCTIONAL_DIR}/${sample}_humann/${sample}_merged_pathabundance.tsv" \
  --output "${FUNCTIONAL_DIR}/${sample}_humann/${sample}_pathabundance_relab.tsv" \
  --units relab

echo "=== [Step 6/6] Metagenomic Assembly & Contig Binning (MEGAHIT) ==="
megahit \
  -1 "${CLEAN_R1}" \
  -2 "${CLEAN_R2}" \
  -t "${THREADS}" \
  --min-contig-len 1000 \
  -o "${ASSEMBLY_DIR}/${sample}_megahit"

echo "=== Pipeline Completed Successfully for ${sample} ==="
