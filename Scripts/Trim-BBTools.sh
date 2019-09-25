source ./Config.sh
export PATH="$bbtools/:$PATH"

# TRIM ADAPTORS / CONTAMINANTS -------------------------------------------------
mkdir -p $trimDir
for i in `awk -F"\t" '{print $1}' $inputSRAFile`; 
	do 
	#Untrimmed Single Read Files
	if [[ -f $rawDir/$i.fastq.gz && ! -f $trimDir/stats/$i.qtrim.txt ]]
	then
		echo $(date) "Trimming Started $i"
		$bbtools/bbduk.sh -Xmx$trimmingMemory\g \
			in=$rawDir/$i.fastq.gz \
			out=$trimDir/$i\_atrim.fastq \
			ref=$adapterFile \
			ktrim=r \
			mink=11 \
			k=23 \
			hdist=1 \
			stats=$trimDir/stats/$i.adaptrim.txt \
			overwrite=t tpe tbo && \
		rm $rawDir/$i.fastq.gz    
		$bbtools/bbduk.sh -Xmx$trimmingMemory\g \
			in=$trimDir/$i\_atrim.fastq \
			out=$trimDir/$i\_AQtrim.fastq \
			qtrim=rl \
			trimq=10 \
			stats=$trimDir/stats/$i.qtrim.txt \
			overwrite=t && \
		rm $trimDir/$i\_atrim.fastq;
		echo $(date) "Trimming Ended $i"
	elif [[ ! -f $rawDir/$i.fastq.gz && -f $trimDir/stats/$i.qtrim.txt ]];  
		then
		echo $(date) "$i already trimmed. Skipping Trimming for  $i"

	#Untrimmed Pair-End Files
	elif [[ -f $rawDir/$i\_1.fastq.gz && ! -f $trimDir/stats/$i.qtrim.txt ]];  
		then
		echo $(date) "Trimming Started $i" 
		$bbtools/bbduk.sh -Xmx$trimmingMemory\g \
			in1=$rawDir/$i\_1.fastq.gz \
			in2=$rawDir/$i\_2.fastq.gz \
			out1=$trimDir/$i\_atrim_1.fastq \
			out2=$trimDir/$i\_atrim_2.fastq \
			ref=$adapterFile \
			ktrim=r \
			mink=11 \
			k=23 \
			hdist=1 \
			stats=$trimDir/stats/$i.adaptrim.txt \
			overwrite=t \
			tpe tbo  && \
		rm $rawDir/$i\_1.fastq.gz && rm $rawDir/$i\_2.fastq.gz
		$bbtools/bbduk.sh -Xmx$trimmingMemory\g \
			in1=$trimDir/$i\_atrim_1.fastq \
			in2=$trimDir/$i\_atrim_2.fastq \
			out1=$trimDir/$i\_AQtrim_1.fastq \
			out2=$trimDir/$i\_AQtrim_2.fastq \
			qtrim=rl \
			trimq=10 \
			overwrite=t \
			stats=$trimDir/stats/$i.qtrim.txt && \
		rm $trimDir/$i\_atrim_1.fastq && rm $trimDir/$i\_atrim_2.fastq
		echo $(date) "Trimming Ended $i"

	elif [[ ! -f $rawDir/$i\_1.fastq.gz && -f $trimDir/stats/$i.qtrim.txt ]];  
	then
		echo $(date) "$i already trimmed. Skipping Trimming for  $i"
	else
		echo $(date) "Cannot Trim $i. No fastq file found"
	fi
done
