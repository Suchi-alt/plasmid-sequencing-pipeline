#!/bin/bash

# Exit on error
set -e
set -x  # Enable debugging for detailed output

# Define file names
reference_file="pl2070.fa"
merged_fastq="merged.fastq"
filtered_fastq="filtered_output.fastq"
aligned_sam="aligned_reads.sam"
aligned_bam="aligned_reads.bam"
sorted_bam="sorted_reads.bam"
unmapped_bam="unmapped.bam"
unmapped_fasta="unmapped.fasta"
unmapped_aligned_sam="unmapped_aligned.sam"
unmapped_sorted_bam="unmapped_sorted.bam"
mapped_vcf="mapped_variants.vcf.gz"
unmapped_vcf="unmapped_variants.vcf.gz"
renamed_unmapped_vcf="unmapped_variants_renamed.vcf.gz"
final_vcf="final_variants.vcf.gz"
filtered_variants_vcf="filtered_variants.vcf"
filtered_variants_vcf_gz="filtered_variants.vcf.gz"
consensus_file="consensus.fa"

# Step 1: Merge FASTQ files
echo "Merging FASTQ files..."
if ls *.fastq 1> /dev/null 2>&1; then
    cat *.fastq > "$merged_fastq"
else
    echo "No FASTQ files found to merge!"
    exit 1
fi

# Step 2: Quality filtering (with Porechop)
echo "Running Porechop for barcode removal..."
if command -v porechop > /dev/null; then
    porechop -i "$merged_fastq" -o "$filtered_fastq" --format fastq -v 2 | tee porechop.log
else
    echo "Porechop not found! Please install it or add it to your PATH."
    exit 1
fi

# Step 3: Generate quality statistics using NanoPlot
echo "Generating quality statistics with NanoPlot..."
if command -v NanoPlot > /dev/null; then
    NanoPlot --fastq "$filtered_fastq" --outdir nanoplot_output --verbose | tee nanoplot.log
else
    echo "NanoPlot not found! Please install it or add it to your PATH."
    exit 1
fi

# Step 4: Align reads to reference using Minimap2
echo "Aligning reads to reference using Minimap2..."
if command -v minimap2 > /dev/null; then
    minimap2 -ax map-ont "$reference_file" "$filtered_fastq" > "$aligned_sam"
else
    echo "Minimap2 not found! Please install it or add it to your PATH."
    exit 1
fi

# Step 5: Convert SAM to BAM, sort, and index
echo "Processing BAM file..."
if command -v samtools > /dev/null; then
    samtools view -Sb "$aligned_sam" > "$aligned_bam"
    samtools sort "$aligned_bam" -o "$sorted_bam"
    samtools index "$sorted_bam"
else
    echo "Samtools not found! Please install it or add it to your PATH."
    exit 1
fi

# Step 6: Validate alignment
echo "Validating alignment with samtools flagstat..."
samtools flagstat "$sorted_bam" > alignment_stats.txt

# Step 7: Extract unmapped reads
echo "Extracting unmapped reads..."
samtools view -b -f 4 "$sorted_bam" > "$unmapped_bam"
samtools fasta "$unmapped_bam" > "$unmapped_fasta"

# Step 8: Align unmapped reads to reference
echo "Aligning unmapped reads to reference using Minimap2..."
minimap2 -ax map-ont "$reference_file" "$unmapped_fasta" > "$unmapped_aligned_sam"

# Step 9: Convert unmapped SAM to sorted BAM
echo "Converting unmapped SAM to sorted BAM..."
samtools view -bS "$unmapped_aligned_sam" | samtools sort -o "$unmapped_sorted_bam"
samtools index "$unmapped_sorted_bam"

# Step 10: Variant calling for mapped reads
echo "Calling variants for mapped reads..."
bcftools mpileup -d 300000000 --no-BAQ --min-BQ 0 -Ou -f "$reference_file" "$sorted_bam" | bcftools call -c -M --ploidy 1 -Oz -o "$mapped_vcf"
bcftools index "$mapped_vcf"

# Step 11: Variant calling for unmapped reads
echo "Calling variants for unmapped reads..."
bcftools mpileup -d 300000000 --no-BAQ --min-BQ 0 -Ou -f "$reference_file" "$unmapped_sorted_bam" | bcftools call -c -M --ploidy 1 -Oz -o "$unmapped_vcf"
bcftools index "$unmapped_vcf"

# Step 12: Rename the sample in the unmapped VCF
echo "Renaming sample in unmapped VCF..."
bcftools reheader --samples <(echo unmapped) "$unmapped_vcf" -o "$renamed_unmapped_vcf"
bcftools index "$renamed_unmapped_vcf"

# Step 13: Merging variants from both calls
echo "Merging variants from mapped and renamed unmapped VCFs..."
bcftools merge -o "$final_vcf" "$mapped_vcf" "$renamed_unmapped_vcf"
bcftools index "$final_vcf"

# Step 14: Filter the VCF to include only variant records (SNPs, indels)
echo "Filtering the VCF to include only variant records..."
bcftools view -v snps,indels "$final_vcf" > "$filtered_variants_vcf"

# Step 15: Compress the filtered VCF with bgzip
echo "Compressing the filtered VCF..."
bgzip -c "$filtered_variants_vcf" > "$filtered_variants_vcf_gz"
bcftools index "$filtered_variants_vcf_gz"

# Step 16: Generate consensus sequence applying only the variants
echo "Generating consensus sequence as nucleotide sequence using filtered variants..."
bcftools consensus -f "$reference_file" "$filtered_variants_vcf_gz" > "$consensus_file"

# Step 17: Show the 3 variants (optional, for verification)
echo "Showing the detected variants:"
bcftools view "$filtered_variants_vcf_gz"

echo "Pipeline completed successfully!"
