# RNA-Seq-Pipeline
Built under Bash version 4.3.48 and R version 3.5.1

### Dependencies 
Unix Packages: SRA Toolkit version 2.9.0, BBTools v38, STAR v2.6.9c, RSEM v1.3.1,  GNU parallel 2017 (optional)<br/>
R Packages: tximport, EdgeR, limma, ggplot, ggpubr, reshape2

### System Requirement
Minimum: 16GB RAM <br/>
Recommended: 32GB+ RAM

### Prerequisite Files
1. FASTA: Genome Assembly of your choice as a single FASTA file. <br/>
2. GTF: Gene Annotation File in GTF format
3. TX2GENE: A Transcript ID to Gene ID Conversion table for RSEM. IMPORTANT: Has to match IDs in your GTF file.

### Usage
1. Create a Tab-delimited file with SRA Accession IDs in first column <br/>
2. Specify your inputs in Config.sh File <br/>
3. Navigate to the RNA-Seq-Pipeline directory 
4. Run the following command: <br/>
```bash
./Run.sh
```


