# Config / Input File ------------------------------------------------

# Tab Delimited Files Specifying Sequence Read IDs in first column; Unique Sample Names in column 2; Phenotype in column 3.
inputSRAFile="/media/Linux/Workspace/BvPBm/pheno.txt"

# This is where your files will be saved
workDir="/media/Linux/Workspace/BvPBm"

#Your Project Name
project="my_project"

#Required Packages
	#Specify path to Packages: sra-toolkit, BBTools, STAR and RSEM
	sra_toolkit="/media/Linux/Packages/sratoolkit.2.9.0"
	bbtools="/media/Linux/Packages/BBTools-v38"
	star="/media/Linux/Packages/STAR-2.6.0c"
	rsem="/media/Linux/Packages/RSEM-1.3.1"

# Reference Genome File (FASTA) - Required only if you lack the RSEM reference / Star Index
refGenomeFile="/media/Documents/1.MetaAnalysis/Genome/Homo_sapien/hg38.v27/Genome/hg38_genome.fa"  

#Reference Gene Annotation File (GTF) - Required only if you lack the RSEM reference / Star Index
geneAnnotationFile="/media/Documents/1.MetaAnalysis/Genome/Homo_sapien/hg38.v27/Genome/gencode.v27.annotation.gtf"

# Transcript to Gene Conversion Table - Transcript ID in Col 1, Gene ID in Col2
tx2gene="/media/Documents/1.MetaAnalysis/Genome/Homo_sapien/hg38.v27/Genome/tx2gene"

# Gene Summarisation Parameters
	#Phenotypes to Compare separated by comma i.e "Diseased-Control,Treated-Control"
	contrast="Case-Control, Case2-Control"

#Computational Power
	#RAM in Gigabytes to be used for Trimming
	trimmingMemory=1

	#Number of Processor Cores to use for Alignment etc.
	numthreads=1

	#Your Machine's total RAM minus 2-4 GB
	totalRAM=14

# Advanced Settings -----------------------------------------------------------------------------

  # TestMode will download first  N reads, this is useful for testing. 
  testMode="false"
  maxReads=200
  
  #Optional Package Parallel - Requires GNU Parallel Package
	parallel_download="false"
  
  # You change this file for specific library prep used e.g. truseq, nextera, etc
  adapterFile="$bbtools/resources/adapters_no_transposase.fa.gz"
  
	# STAR 2 pass alignment improves detection of splice variants. But increases processing times.
	secondpass="false"  

  # Location of Star Index 
  starIndexDir="$working_Dir/starIndex"

  # Location of RSEM Reference
  rsemRefDir="$working_Dir/rsemReference"
  	#Prefix (if any) for RSEM reference files
  	rprefix="rsemRef"  

	#Specify following if Star Index is missing
		#Maximum read length-1 (typically 99 for illumina sequencers)
		overhang=99
			# Optional: Indexing SNPs takes a very long time for larger genomes 
			indexSNPs="false"
			# If you set indexSNP = "true",  Specify VCF file containing consensus SNPs
			consensusSNPFile="path/to/file.vcf"


	#Filter out Specific Genes
	filterGenes="false"
		#Genes to Filter out - 1 gene per line
		genefilter="/path/to/genelist"

#Edit these if you have limited disk space --------
	# Raw File Save Location
	rawDir="$workDir/1.Raw"

	# Trimmed File Save Location
	trimDir="$workDir/2.Trimmed"

	# Aligned File Save Location
	bamDir="$workDir/3.Aligned"

	# Read Count Save Location
	rsemDir="$workDir/4.ReadCounts"

	# Differential Expression Save Location
	diffExpDir="$workDir/5.Differential_Expression"
	
	# A folder to save large temporary files. 
	sra_tempdir="$workDir/tmp"
