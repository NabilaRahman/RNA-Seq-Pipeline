source ./Config.sh

# VALIDATION ROUTINES -------------------------------------------------
set -o nounset
command || { echo $(date) "ERROR: Command failed. Check your script for errors"; exit; }

# Working Directory
if [[ ! -d $workDir || -z $workDir ]]
	then
	echo $(date) "Your specified working directory does not exist. Please input a valid path."
	exit
fi

#Read Sequence Accession ID
if [[ ! -f $inputSRAFile || -z $inputSRAFile ]]; then echo $(date) "Accession file does not exist. Create a tab delimited file with SRA accession in first column."; fi

# Check if sra-toolkit is installed
if [[ ! -f "$sra_toolkit/bin/fastq-dump" || -z $sra_toolkit/bin/fastq-dump ]]
	then
	echo $(date) "ERROR: Downloading SRA files from NCBI requires 'sra-toolkit' package. Please install 'sra-toolkit 2.9.0' or higher and specify the location in the config file"
	exit
fi

# PARALLEL
paraexist=`dpkg -s parallel | grep 'Status: install ok installed'`
if [[ $parallel_download == "true" && -z $paraexist || $paraexist == "" ]]
	then
	echo $(date) "ERROR: UNIX package 'GNU parallel' not found. Set 'parallel_download' to \"false\" to download sequentially. Otherwise install the GNU parallel 2017 or higher."
	exit
fi

# BBTools
if [[ ! -f "$bbtools/bbduk.sh" || -z $bbtools/bbduk.sh ]] 
then
	echo $(date) "ERROR: BBTools not found. Please install 'BBTools v38' or higher and specify the package location in the config file"
	exit
fi

#STAR
if [[ ! -f $star/source/STAR || -z $star/source/STAR ]]; then
	echo $(date) "ERROR: STAR package not found.  Please install 'STAR 2.6.0c' or higher and specify the package location in the config file"
	exit
fi

#RSEM 
if [[ ! -f $rsem/rsem-calculate-expression || -z $rsem/rsem-calculate-expression ]]; then
	echo $(date) "ERROR: RSEM package not found. Please install 'RSEM 1.3.1' or higher and specify the package location in the config file"
	exit
fi

#STAR Index
if [[ ! -f $starIndexDir/SAindex || -z $starIndexDir/SAindex ]]
	then 
	if [[ $indexSNPs == "true" ]] && [[ ! -f $consensusSNPFile ||  -z $consensusSNPFile ]]
		then
		echo "ERROR: STAR index missing. To create STAR index with SNP index specify genome Consensus File (.VCF), or set 'indexSNPs' to 'false'"
		exit
	fi
	
	if [ $totalRAM -gt 32  ]
		then sparse=1
	elif	[ $totalRAM -eq 32  ]
		then sparse=3
	elif [ $totalRAM -ge 16 ]
		then sparse=4
	elif [ $totalRAM -ge 8 ]
		then sparse=5
	else 
		echo "ERROR: Need more memory to run STAR indexer" 
		exit 
	fi
	
	for i in refGenomeFile geneAnnotationFile 
		do
		if [[ ! -f ${!i} || -z ${!i} ]]
			then
			echo "ERROR: STAR Index does not exist at specified directory 'starIndexDir', please specify location of your STAR Index or alternatively specify '$i', to create a new STAR Index"
			exit
		fi
	done
fi


#Tx2Gene
if [[ ! -f $tx2gene || -z $tx2gene ]]
	then
	echo "ERROR: 'tx2gene' does not exist at specified directory. Please create a Transcript ID to Gene ID conversion table that matches the IDs in your Gene Annotation File. This is required for Read Counting by the RSEM package."
	exit
	fi
fi

#RSEM REFERENCE
if [[ ! -f $rsemRefDir/$rprefix.n2g.idx.fa || -z $rsemRefDir/$rprefix.n2g.idx.fa ]]
	then 
	for i in refGenomeFile geneAnnotationFile
		do
		if [[ ! -f ${!i} || -z ${!i} ]]
			then
			echo "ERROR: RSEM Reference does not exist at specified directory 'rsemDir', please specify location of RSEM Reference or alternatively specify '$i', to create a new RSEM Reference"
			exit
		fi
	done
fi

#BOOLEANS
#Boolean
for i in testMode parallel_download secondpass indexSNPs filterGenes
do
	if [[ ${!i} != "true" ]]
		then
			if [[ ${!i} != "false" ]]
				then 
				echo "ERROR: '$i'  must be \"true\" or \"false\""
				exit
			fi	
	fi
done

#NUMERIC
integer='^[0-9]+$'
for i in trimmingMemory numthreads totalRAM overhang maxReads
	do
	if ! [[ ${!i} =~ $integer ]] ; then
	   echo "'ERROR: '$i' is not an integer. Please set an integer value for '$i'"
	   exit
	fi
done

# Process inputs --------------------------------------------------
# Set SRA save directory
if [[ "" != "$sra_tempdir" &&  -z $sra_tempdir && -f $HOME/.ncbi/user-settings.mkfg ]]
	then 
	mkdir -p $sra_tempdir
	echo "/repository/user/main/public/root = '$sra_tempdir'" > $HOME/.ncbi/user-settings.mkfg
fi  

validated="true"
mkdir -p $diffExpDir
export validated
export sparse
