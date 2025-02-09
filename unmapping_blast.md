
## Extraction of Unmapped Reads with minimap2 + BLASTn Classification

Commands for host DNA removal -> optaining non-host reads (unmapped reads) and for BLASTn classification 

### 1. Convert BAM to FASTQ

```bash
samtools fastq input.bam > output.fastq.gz
```

### 2. Map Reads and Extract Unmapped Reads

```bash
for i in *.gz; do \
  ref=$(echo "${i}" | cut -f1 -d"_"); \
  out=$(echo "${i}" | cut -f1 -d"."); \
  ~/Software/minimap2/minimap2 -ax map-hifi -t 15 \
    ~/Scratch/mites_proj/raw_data/genomes/"${ref}".fasta "${i}" | \
  samtools view -b -f 4 - | \
  samtools sort -o ~/Scratch/mites_proj/localBlast/raw_reads/"${out}".unmapped.bam; \
done
```

here for a single file (Ass_reads.fastq.gz) 

```bash
~/Software/minimap2/minimap2 -ax map-hifi -t 15 \
  ~/Scratch/mites_proj/raw_data/genomes/Ass.fasta \
  Ass_reads.fastq.gz | samtools view -b -f 4 - | \
  samtools sort -o ~/Scratch/mites_proj/localBlast/raw_reads/Ass_reads.unmapped.bam
```

### 3. Convert BAM to FASTA
```bash
samtools fasta input.bam > output.fasta
```

### 4. Create BLASTn database 

```bash
makeblastdb -in SILVA_138.2_SSURef_NR99_tax_silva.fasta -dbtype nucl -out silvaDB
```

### 5. Run BLASTn on unmapped reads 

here: e-value cutoff is to 1e-10; and top 5 matching sequences per query 
```bash
for i in /home/tkoch/Scratch/mites_proj/localBlast/raw_reads/*.fasta; do
    b=$(basename "$i" .unmapped.fasta);
    blastn -query "$i" \
           -db /home/tkoch/Scratch/mites_proj/localBlast/db/silvaDB \
           -evalue 1e-10 \
           -max_target_seqs 5 \
           -outfmt "6 qseqid sseqid stitle evalue bitscore pident length mismatch gapopen qstart qend sstart send" \
           -out /home/tkoch/Scratch/mites_proj/localBlast/results_blast/unmapped_blast/"$b".blast.results \
           -num_threads 10;
done
```

### 6. Formatting BLASTn output 

to filter BLASTn output files by high percent identity >90; alignment lengths >650 (bp) 
```bash
for i in /home/tete/Documents/Project/Blast/Unmapped/*.tsv; do 
    b=$(basename "$i" _reads.blast.results.tsv)
    
    cat "$i" | \
    cut -f1,3- | \
    tr " " ";" | \
    awk '{if ($5 > 90 && $6 > 650) print $0}' | \
    sort -k1,1 -k5,5nr | \
    sort -k1,1 -u -s > /home/tete/Documents/Project/Blast/Unmapped/Formatted/"$b"_formatted_unmapped.txt

done
```

