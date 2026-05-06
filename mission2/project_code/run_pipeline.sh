#!/usr/bin/env bash
# Mission 2 – shell pipeline
# Run from project root: bash project_code/run_pipeline.sh

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DATA="$ROOT/data"
WORK="$ROOT/work"

SAMTOOLS=/usr/local/bin/samtools
BEDTOOLS=/Users/angela/opt/anaconda3/bin/bedtools
BIOAWK=/Users/angela/opt/anaconda3/bin/bioawk
GREP=/usr/local/bin/ggrep   # GNU grep (supports -P for Perl regex / \t)
SED=/usr/local/bin/gsed      # GNU sed  (macOS sed doesn't support \t in [])

mkdir -p "$WORK"

echo "======================================================"
echo "STEP 1 – Start codons (+ strand, TSL 1)"
echo "======================================================"
$GREP -P '\tstart_codon\t' "$DATA/gencode.gtf" \
  | $GREP -P '\t\+\t' \
  | $GREP 'transcript_support_level "1"' \
  | $SED -e 's/\t[^\t]*transcript_id "\([^"]*\)".*$/\t\1/g' \
  > "$WORK/gencode-start.gtf"
echo "  Lines: $(wc -l < "$WORK/gencode-start.gtf")"

echo ""
echo "======================================================"
echo "STEP 2 – Plus-strand exons"
echo "======================================================"
$GREP -P '\texon\t' "$DATA/gencode.gtf" \
  | $GREP -P '\t\+\t' \
  | $SED -e 's/\t[^\t]*transcript_id "\([^"]*\)".*$/\t\1/g' \
  > "$WORK/gencode-plusexon.gtf"
echo "  Lines: $(wc -l < "$WORK/gencode-plusexon.gtf")"

echo ""
echo "======================================================"
echo "STEP 3 – Exons containing start codons"
echo "======================================================"
$BEDTOOLS intersect \
  -a "$WORK/gencode-start.gtf" \
  -b "$WORK/gencode-plusexon.gtf" \
  -wa -wb \
  | awk -F'\t' -v OFS='\t' '$9 == $18 { print $10, $13-1, $14, $18, $4-1, $16 }' \
  | sort -k1,1 -k2,3n -k4,4 \
  > "$WORK/gencode-exons-containing-startcodon.bed"
echo "  Lines: $(wc -l < "$WORK/gencode-exons-containing-startcodon.bed")"

echo ""
echo "======================================================"
echo "STEP 4 – Filter RPF BAM (+ strand, length >= 25 nt)"
echo "======================================================"
( $SAMTOOLS view -H "$DATA/RPF-siLuc.bam";
  $SAMTOOLS view -F20 "$DATA/RPF-siLuc.bam" \
    | $BIOAWK -c sam '{ if (length($seq) >= 25) print $0 }' ) \
  | $SAMTOOLS view -b -o "$WORK/filtered-RPF-siLuc.bam"
ls -lh "$WORK/filtered-RPF-siLuc.bam"

echo ""
echo "======================================================"
echo "STEP 5 – 5'-end coverage (bedgraph)"
echo "======================================================"
$BEDTOOLS genomecov \
  -ibam "$WORK/filtered-RPF-siLuc.bam" \
  -bg -5 \
  > "$WORK/fivepcounts-RPF-siLuc.bed"
echo "  Lines: $(wc -l < "$WORK/fivepcounts-RPF-siLuc.bed")"

echo ""
echo "======================================================"
echo "STEP 6 – Intersect 5'-ends with start-codon exons"
echo "======================================================"
$BEDTOOLS intersect \
  -a "$WORK/fivepcounts-RPF-siLuc.bed" \
  -b "$WORK/gencode-exons-containing-startcodon.bed" \
  -wa -wb -nonamecheck \
  > "$WORK/fivepcounts-filtered-RPF-siLuc.txt"
echo "  Lines: $(wc -l < "$WORK/fivepcounts-filtered-RPF-siLuc.txt")"

echo ""
echo "Shell pipeline done. Run: python3 project_code/plot.py"
