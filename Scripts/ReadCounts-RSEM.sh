source ./Config.sh
rsemRef=$rsemRefDir/$rprefix

# Check For RSEM Index
if [[ ! -f $rsemRefDir/$rprefix.n2g.idx.fa ]] 
	then
mkdir -p $rsemRef
	$rsem/rsem-prepare-reference \
		--gtf $geneAnnotationFile \
		--star $refGenomeFile \
		$rsemRef/$rprefix \
		--num-threads $numthreads \
		--star-sjdboverhang 99
else
 	echo $(date) "RSEM Reference Found..."
fi


# READ COUNTING - RSEM
mkdir -p $rsemDir

for i in `ls $bamDir | grep -E ".Aligned.toTranscriptome.out.bam\$" | sed 's/\..*//' | uniq`
	do
	if [[ ! -f $rsemDir/$i/$i.isoforms.results && -f $bamDir/$i.single.Aligned.toTranscriptome.out.bam ]]
		then
		echo $(date) "Quantifying $i as Single Read"
		mkdir -p $rsemDir/$i
		$rsem/rsem-calculate-expression \
			--bam $bamDir/$i.single.Aligned.toTranscriptome.out.bam \
			$rsemRef \
			-p $numthreads \
			--no-bam-output $rsemDir/$i/$i

	elif [[ ! -f $rsemDir/$i/$i.isoforms.results && -f $bamDir/$i.paired.Aligned.toTranscriptome.out.bam  ]]; then
		echo $(date) "Quantifying $i as Paired End"
		mkdir -p $rsemDir/$i
		$rsem/rsem-calculate-expression \
			--bam \
			--paired-end $bamDir/$i.paired.Aligned.toTranscriptome.out.bam \
			$rsemRef \
			-p $numthreads \
			--no-bam-output $rsemDir/$i/$i
	else
		echo $(date) "Skipping $i, already quantified"
	fi
done
