source ./Config.sh
export PATH="$sra_toolkit/:$PATH"

# DOWNLOAD & CONVERT RAW RNA-SEQ FILES ---------------------------
mkdir -p $rawDir
echo $(date) "Downloading SRA files. This may take a long time."

testParam=""
sra='sra'
if [[ $testMode == "true" ]]
	then 
	testParam="-X $maxReads"
	sra='sra.cache'
fi


if [[ $parallel_download == "true" ]]
	then
	export rawDir
	export sra_toolkit
	export testParam
	export sra_tempdir
	export sra
	echo $(date) "Downloading in Parallel"
	awk -F"\t" '{print $1}' $inputSRAFile | parallel --no-notice -a - '/
	if [[ ! -f $rawDir/{}.log ]]
    	then
    	$sra_toolkit/bin/fastq-dump /
			--split-3 $testParam -O $rawDir --gzip {} && /
		echo $(date) "Successfully downloaded: {}.fastq,gz" > $rawDir/{}.log && /
		rm $sra_tempdir/sra/{}.$sra 2> /dev/null;
	else
    	echo $(date) "Skipping {}. Already Downloaded";
    fi /
	'

elif [[ $parallel_download != "true" ]]
	then
	echo $(date) "Downloading one by one"
	for i in `awk -F"\t" '{print $1}' $inputSRAFile`
		do
		if [[ ! -f $rawDir/$i.log ]]; then
			echo $(date) "$i: Downloading... "
			$sra_toolkit/bin/fastq-dump --split-3 $testParam -O $rawDir --gzip $i
			echo $(date) "$i: Download Complete" > $rawDir/$i.log
			rm $sra_tempdir/sra/$i.$sra 2> /dev/null
		else
  			echo $(date) "Skipping {}. Already Downloaded"
		fi
done
fi 
