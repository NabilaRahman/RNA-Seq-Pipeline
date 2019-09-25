source ./Scripts/ValidateInput.sh
functions="./Scripts/Functions.R"

export inputFile  species project tx2gene gene2sym filterGenes genefilter functions rsemDir diffExpDir contrast project

if [[ $validated == "true" ]]
	then
	./Scripts/DownloadFastQ.sh && \
	./Scripts/Trim-BBMap.sh && \
	./Scripts/Align-STAR.sh && \
	./Scripts/ReadCounts-RSEM.sh && \
	Rscript ./Scripts/SummariseGene.R
fi
