source ./Config.sh
STAR="$star/bin/Linux_x86_64/STAR" 

# STAR - ALIGN TO REFERENCE GENOME ----------------------------------------------------------------------

if [[ ! -f $starIndexDir/SAindex || -z $starIndexDir/SAindex ]]
	then
	echo "Creating Star Index"
	#Index SNPs
	if [[ $indexSNPs == "true" ]]
		then
		genomecon="--genomeConsensusFile $consensusSNPFile"
	else
		genomecon=""
	fi

	#Create Star Index if required
	mkdir -p $starIndexDir
	$STAR --runThreadN $numthreads /
		--runMode genomeGenerate /
		--genomeDir $starIndexDir /
		--genomeFastaFiles $refGenomeFile /
		--genomeSAsparseD $sparse $genomecon /
		--sjdbGTFfile $geneAnnotationFile /
		--sjdbOverhang $overhang
fi

#First Pass
mkdir -p $bamDir

saveBam='None'
singlePrefix='temp'
pairedPrefix='temp'
genomeKeep="--genomeLoad LoadAndKeep"
rmTrimSingle=''
rmTrimPaired=''
genomeloaded="false"
if [[ $secondpass != "true" ]]
	then
	saveBam='BAM SortedByCoordinate --quantMode TranscriptomeSAM  --outSAMunmapped Within'
	singlePrefix="single"
	pairedPrefix="paired"
	genomeKeep=""
	rmTrimSingle='&&  rm $trimDir/$i\_AQtrim.fastq'
	rmTrimPaired='&& rm $trimDir/$i\_AQtrim_1.fastq && rm $trimDir/$i\_AQtrim_2.fastq'
fi


echo $(date) "Aligning to Reference Genome - STAR First Pass"
for i in `awk -F"\t" '{print $1}' $inputSRAFile`; 
do 
	if [[ -f $trimDir/$i\_AQtrim.fastq && -f $trimDir/stats/$i.qtrim.txt && ! -f $bamDir/$i.$singlePrefix.SJ.out.tab ]] 
		then
		echo $(date) "$i: First Pass - Start"
		$STAR --readFilesIn $trimDir/$i\_AQtrim.fastq \
			--runThreadN $numthreads \
			--runMode alignReads \
			--genomeDir $starIndexDir \
			--outSAMtype $saveBam \
			--outFileNamePrefix $bamDir/$i.$singlePrefix. \
			$genomeKeep $rmTrimSingle
			if [[ "$singlePrefix" == "temp" ]]
				then
				genomeloaded="true"
			fi
			echo $(date) "$i: First Pass - End"
			echo
	elif [[ -f $bamDir/$i.$singlePrefix.SJ.out.tab ]]
		then
		 echo $(date) "$i already aligned. Skipping First Pass Alignment for  $i"

	elif [[ -f $trimDir/$i\_AQtrim_1.fastq && -f $trimDir/stats/$i.qtrim.txt && ! -f $bamDir/$i.$pairedPrefix.SJ.out.tab ]];     
		then
		echo $(date) "$i: First Pass - Start"
		$STAR --readFilesIn $trimDir/$i\_AQtrim_1.fastq $trimDir/$i\_AQtrim_2.fastq \
		--runThreadN $numthreads \
		--runMode alignReads \
		--genomeDir $starIndexDir \
		--outSAMtype $saveBam \
		--outFileNamePrefix $bamDir/$i.$pairedPrefix. \
		$genomeKeep
		if [[ "$pairedPrefix" == "temp" ]]
				then
				genomeloaded="true"
		fi
		echo $(date) "$i: First Pass - End"
		echo
	elif [[ -f $bamDir/$i.$pairedPrefix.SJ.out.tab ]];
		then
		echo $(date) "$i already aligned. Skipping First Pass Alignment for  $i"
	else
	  	echo $(date) "Cannot Align $i, No fastq file found"
	fi
done

if [[ $secondpass == "true" && $genomeloaded == "true" ]]
	then
	#Remove Genome from memory  
	$STAR --genomeDir $starIndexDir --genomeLoad Remove;
	cat $bamDir/*.temp.*SJ.out.tab > $bamDir/all.SJ.out.tab
fi
	
	#Second Pass
if [[ $secondpass == "true" ]]
	then
	for i in `awk -F"\t" '{print $1}' $inputSRAFile` 
	do 
		if [[ -f $trimDir/$i\_AQtrim.fastq && -f $trimDir/stats/$i.qtrim.txt && ! -f $bamDir/$i.single.SJ.out.tab ]] 
			then
			echo $(date) "$i: Second Pass - Start"
			$STAR --readFilesIn $trimDir/$i\_AQtrim.fastq \
			--runThreadN $numthreads \
			--runMode alignReads \
			--sjdbFileChrStartEnd $bamDir/all.SJ.out.tab \
			--genomeDir $starIndexDir \
			--outSAMtype BAM SortedByCoordinate \
			--quantMode TranscriptomeSAM \
			--outSAMunmapped Within  \
			--outFileNamePrefix $bamDir/$i.single. \
			&&  rm $trimDir/$i\_AQtrim.fastq
			echo $(date) "$i: Second Pass - End"
			echo		
		elif [[ -f $bamDir/$i.single.SJ.out.tab ]]
			then
			echo $(date) "$i already aligned. Skipping Second Pass Alignment for  $i"

		elif [[ -f $trimDir/$i\_AQtrim_1.fastq && -f $trimDir/stats/$i.qtrim.txt && ! -f $bamDir/$i.paired.SJ.out.tab ]]
			then
			echo $(date) "$i: Second Pass - Start"
			$STAR --readFilesIn $trimDir/$i\_AQtrim_1.fastq $trimDir/$i\_AQtrim_2.fastq  \
			--runThreadN $numthreads \
			--runMode alignReads \
			--quantMode TranscriptomeSAM \
			--sjdbFileChrStartEnd $bamDir/all.SJ.out.tab \
			--genomeDir $starIndexDir \
			--outSAMtype BAM SortedByCoordinate \
			--outFileNamePrefix $bamDir/$i.paired.  \
			--outSAMunmapped Within \
			&& rm $trimDir/$i\_AQtrim_1.fastq && rm $trimDir/$i\_AQtrim_2.fastq
			echo $(date) "$i: Second Pass - End"
			echo
		elif [[ -f $bamDir/$i.paired.SJ.out.tab ]]
			then
			echo $(date) "$i already aligned. Skipping Second Pass Alignment for  $i"
		else
			echo $(date) "Cannot Align $i, No fastq file found"
		fi
	done
fi
