## For Length and Quality Check

### 1. Seqkit for read length information
for raw read length (PacBio HiFi sequences); consists of both host DNA and unmapped reads 

```bash
for i in /home/tkoch/Scratch/mites_proj/raw_data/reads/*.gz; do
	seqkit stats "$i" >> /home/tkoch/Scratch/mites_proj/raw_data/reads/stats_reads/all_stats_raw_reads.txt 
done
```
for unmapped reads 
```bash
for i in /home/tkoch/Scratch/mites_proj/localBlast/fasta_unmapped_reads/*.fasta; do
	seqkit stats "$i" >> /home/tkoch/Scratch/mites_proj/localBlast/fasta_unmapped_reads/stats_reads/all_stats_unmapped_reads.txt 
done
```
