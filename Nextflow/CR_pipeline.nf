#!/usr/bin/env nextflow

inputDir = params.inputDir
dataDir = params.dataDir
FragSize = params.FragSize
resultsDir = params.resultsDir
HomerSize = params.HomerSize
pipelineMode = params.pipelineMode
big_task_cpus = params.big_task_cpus

process File_sheet_match {
    label 'small_task'
    input:
    tuple val(sampleId), path(reads)

    output:
    tuple env(dirname), path(reads)
    script:
    """
    dirname=$sampleId #dirname taken from input
    dirname=\${dirname%%_*} #Shortened to appropriate length
    """
}

process MoveToDataDir {
    publishDir "$params.dataDir" , mode: 'copy', overwrite: true
    label 'small_task'

    input:
    tuple val(sampleId), file(reads), val(targetID)

    output:
    path "$targetID", emit: workingDir
    stdout emit: workingDirName

    script:
    """
    mkdir -p $sampleId #Make the data directory
    mv ${reads} $sampleId 
    mv $sampleId $targetID
    cd $targetID
    rawname=${sampleId}*.fastq.gz
    for file in \$rawname ; do
    mv \$file ${targetID}_\$file
    done;
    mkdir -p demultiplexed_reads
    mv *.fastq.gz demultiplexed_reads
    cd ..
    touch $targetID
    printf $targetID
    """
}

process TrimFiles {
  publishDir "${params.dataDir}" , mode: 'copy', overwrite: true
  label 'big_task'

  input:
  path workingDir
  val workingDirName
  
  output:
  path "$workingDir/trimmed_reads/*_R1_*.fq.gz", emit: trimmedfile1
  path "$workingDir/trimmed_reads/*_R2_*.fq.gz", emit: trimmedfile2
  path "$workingDir/trimmed_reads/*R1*_trimming_report.txt", emit: trimreport1
  path "$workingDir/trimmed_reads/*R2*_trimming_report.txt", emit: trimreport2
  path workingDir, includeInputs: true, emit: workingDir
  val workingDirName, emit: workingDirName

  script:
  """
  module purge
  module load Trim_Galore
  cd $workingDir
  file1="demultiplexed_reads/*R1*.fastq.gz"
  file2="demultiplexed_reads/*R2*.fastq.gz"
  cmd="trim_galore --paired \$file1 \$file2 --fastqc --gzip -o ./trimmed_reads"
  eval \${cmd}
  """
}

process MapFiles {
  publishDir "${params.dataDir}/${workingDirName}" , mode: 'copy', overwrite: true, pattern: "*.{bam,bai,log,out}", enabled: "$publish19"
  publishDir "${params.dataDir}/${workingDirName}/${refgenome}" , mode: 'copy', overwrite: true, pattern: "*.{bam,bai,log,out}", enabled: "$publish38"
  label 'big_task'

  input:
  path trimmedfile1
  path trimmedfile2
  path trimreport1
  path trimreport2
  path workingDir
  val workingDirName
  tuple val(refgenome),val(publish19),val(publish38)
  
  output:
  path "*_Hs_mappingstder_${refgenome}.out"
  path "*.samtoolstats.Hs.${refgenome}.log"
  path "*.flagstat.${refgenome}.log"
  path workingDir, includeInputs: true
  tuple val (workingDirName), path ("*_Hs_sorted.MAPQ30.${refgenome}.bam"), path ("*_Hs_sorted_${refgenome}.bam"),path ("*.bai"), emit: deduptuple
  tuple val(refgenome),val(publish19),val(publish38)

  script:
  HGIndex=params.HGIndex[task.process]
  """
  module purge
  module load Bowtie2/2.4.5-GCC-11.3.0 
  module load SAMtools/1.16.1-GCC-11.3.0

  basename="$workingDirName"

  mkdir -p tmpfiles

  cmd="(bowtie2 -p $big_task_cpus --local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -X 700 -x $HGIndex -1 $trimmedfile1 -2 $trimmedfile2 | samtools sort -m 5G -T tmpfiles -O bam -@ $big_task_cpus -o \${basename}_Hs_sorted_${refgenome}.bam) 3>&1 1>&2 2>&3 | tee \${basename}_Hs_mappingstder_${refgenome}.out"
  eval \${cmd}
  
  samtools index -@ $big_task_cpus \${basename}_Hs_sorted_${refgenome}.bam

  samtools stats \${basename}_Hs_sorted_${refgenome}.bam > \${basename}.samtoolstats.Hs.${refgenome}.log
  samtools flagstat \${basename}_Hs_sorted_${refgenome}.bam > \${basename}.flagstat.${refgenome}.log

  samtools view -b -q 30 \${basename}_Hs_sorted_${refgenome}.bam > \${basename}_Hs_sorted.MAPQ30.${refgenome}.bam
  samtools index -@ $big_task_cpus \${basename}_Hs_sorted.MAPQ30.${refgenome}.bam
  """
}

process MultiQC {
  publishDir "${params.dataDir}/${workingDirName}" , mode: 'copy', overwrite: true, pattern: "*multiQC*", enabled: "$publish19"
  publishDir "${params.dataDir}/${workingDirName}/${refgenome}" , mode: 'copy', overwrite: true, pattern: "*multiQC*", enabled: "$publish38"
  label 'small_task'

  input:
  path mappingstder
  path samtoolslog
  path flagstatlog
  path workingDir
  tuple val (workingDirName), path (MAPQ30bam), path (sortedbam),path (bai)
  tuple val(refgenome),val(publish19),val(publish38)
  
  output:
  path workingDir, includeInputs: true, emit: workingDir
  path "${workingDirName}_${refgenome}_multiQC.html"
  path "${workingDirName}_${refgenome}_multiQC_data"
  val workingDirName, emit: workingDirName
  env duppercent, emit: Dedup_percent
  script:
  """
  module purge
  module load MultiQC

  basename=$workingDirName
  
  multiqc -f . -n \${basename}_${refgenome}_multiQC

  duppercent=\$(cut -f14 \${basename}_${refgenome}_multiQC_data/multiqc_general_stats.txt | sort -r | head -2 | tail -1)
  """
}

process Dedup {
  publishDir "${params.dataDir}/${workingDirName}" , mode: 'copy', overwrite: true, pattern: "*dedup*", enabled: "$publish19"
  publishDir "${params.dataDir}/${workingDirName}/${refgenome}" , mode: 'copy', overwrite: true, pattern: "*dedup*", enabled: "$publish38"
  label 'big_task'

  input:
  tuple val (workingDirName), val(Dedup_percent), path (MAPQ30bam), path (sortedbam),path (bai)
  tuple val(refgenome),val(publish19),val(publish38)

  output:
  path "*dedup.MAPQ30.${refgenome}.bam", emit: postDedupMAPQ30BAM
  path "*dedup_${refgenome}.bam"
  path "*dedup*.log", emit: deduppedLogs
  path "*dedup*.bai", emit: dedupBAI
  val workingDirName, emit: workingDirName
  path "metrics_${refgenome}.txt", optional: true

  script:
  
  if (Dedup_percent >= 10)
  """
  basename=$workingDirName
  metrics=\${basename}_dedup_metrics_${refgenome}.txt
  
  module purge
  module load picard

  java -jar \$EBROOTPICARD/picard.jar MarkDuplicates I="./\${basename}_Hs_sorted_${refgenome}.bam" O="./\${basename}_Hs_sorted_dedup_${refgenome}.bam" M="./\${metrics}" VALIDATION_STRINGENCY=SILENT OPTICAL_DUPLICATE_PIXEL_DISTANCE=12000 REMOVE_DUPLICATES=true

  module purge
  module load SAMtools

  #indexeren
  samtools index -@ $big_task_cpus \${basename}_Hs_sorted_dedup_${refgenome}.bam

  samtools stats \${basename}_Hs_sorted_dedup_${refgenome}.bam > \${basename}.samtoolstats.dedup.Hs.${refgenome}.log
  samtools flagstat \${basename}_Hs_sorted_dedup_${refgenome}.bam > \${basename}.dedup.flagstat.${refgenome}.log

  samtools view -b -q 30 \${basename}_Hs_sorted_dedup_${refgenome}.bam > \${basename}_Hs_sorted_dedup.MAPQ30.${refgenome}.bam
  samtools index -@ $big_task_cpus \${basename}_Hs_sorted_dedup.MAPQ30.${refgenome}.bam

  samtools stats \${basename}_Hs_sorted_dedup.MAPQ30.${refgenome}.bam > \${basename}.samtoolstats.dedup.MAPQ30.Hs.${refgenome}.log
  samtools flagstat \${basename}_Hs_sorted_dedup.MAPQ30.${refgenome}.bam > \${basename}.dedup.MAPQ30.flagstat.${refgenome}.log
  """
  
  else
  """
  echo "No deduplication required."
  basename=$workingDirName
  mv \${basename}_Hs_sorted_${refgenome}.bam \${basename}_Hs_sorted_nodedup_${refgenome}.bam
  mv \${basename}_Hs_sorted_${refgenome}.bam.bai \${basename}_Hs_sorted_nodedup_${refgenome}.bam.bai
  echo "Duplicate percentage ($Dedup_percent%) was lower than 10% so no deduplication was performed" > \${basename}.samtoolstats.nodedup.Hs.${refgenome}.log
  mv \${basename}_Hs_sorted.MAPQ30.${refgenome}.bam \${basename}_Hs_sorted_nodedup.MAPQ30.${refgenome}.bam
  mv \${basename}_Hs_sorted.MAPQ30.${refgenome}.bam.bai \${basename}_Hs_sorted_nodedup.MAPQ30.${refgenome}.bam.bai
  echo "Duplicate percentage ($Dedup_percent%) was lower than 10% so no deduplication was performed" > \${basename}.samtoolstats.nodedup.MAPQ30.Hs.${refgenome}.log
  """
}

process SplitFragments {
  publishDir "${params.dataDir}/${workingDirName}" , mode: 'copy', overwrite: true, pattern: "*{low,high}*.{bam,bam.bai}", enabled: "$publish19"
  publishDir "${params.dataDir}/${workingDirName}/${refgenome}" , mode: 'copy', overwrite: true, pattern: "*{low,high}*.{bam,bam.bai}", enabled: "$publish38"
  publishDir "${params.resultsDir}/Fragmentsizes", mode: 'move', overwrite: true, pattern: "*.png"

  label 'big_task'

  input:
  path postDedupMAPQ30BAM
  val workingDirName
  path dedupBAI
  tuple val(refgenome),val(publish19),val(publish38)

  output:
  path "*.low${FragSize}.${refgenome}.bam", emit: SmallFragBAM
  path "*.high${FragSize}.${refgenome}.bam", emit: LargeFragBAM
  path "*${FragSize}.${refgenome}.bam.bai", emit: splitFragBAI
  path postDedupMAPQ30BAM, includeInputs: true, emit: postDedupMAPQ30BAM
  val workingDirName, emit: workingDirName
  path "*.png"

  script:
  """
  filename=$postDedupMAPQ30BAM
  basename=\${filename%%.${refgenome}.bam}
  
  module purge
  module load SAMtools

  samtools view -h \$filename | \
    awk -v LEN=$FragSize '{if (\$9 <= LEN && \$9 >= -(LEN) && \$9 != 0 || \$1 ~ /^@/) print \$0}' | \
    samtools view -bh - > \${basename}.low${FragSize}.${refgenome}.bam
  samtools index \${basename}.low${FragSize}.${refgenome}.bam;

  samtools view -h \$filename | \
    awk -v LEN=$FragSize '{if (\$9 >= LEN && \$9 <= 1000 && \$9 != 0 || \$1 ~ /^@/) print \$0}' | \
    samtools view -bh - > \${basename}.high${FragSize}.${refgenome}.bam
  samtools index \${basename}.high${FragSize}.${refgenome}.bam;

  module purge
  module load deepTools

  bamPEFragmentSize -b \$filename -hist Fragmentsize_\${basename}.${refgenome}.png -T "Fragment size of PE seq data"

  """
}

process BlackFiltering {
  publishDir "${params.dataDir}/${workingDirName}" , mode: 'copy', overwrite: true, pattern: "*blackfiltered*", enabled: "$publish19"
  publishDir "${params.dataDir}/${workingDirName}/${refgenome}" , mode: 'copy', overwrite: true, pattern: "*blackfiltered*", enabled: "$publish38"
  label 'big_task'

  input:
  path SmallFragBAM
  path LargeFragBAM
  val workingDirName
  path postDedupMAPQ30BAM
  path splitFragBAI
  tuple val(refgenome),val(publish19),val(publish38)

  output:
  path "*blackfiltered_${refgenome}.bam", emit: BlackfilteredBAMS
  path "*.MAPQ30_blackfiltered_${refgenome}.bam", emit: PeakCallingBAMS
  path "*blackfiltered_${refgenome}.bam.bai", emit: BlackfilteredBAIS
  stdout emit: BlackFilterBaseName

  script:
  BlackFilterPath=params.BlackFilterPath[task.process]
  """
  filename=$postDedupMAPQ30BAM
  basename=\${filename%%.${refgenome}.bam}
  echo \$basename
  module purge
  module load BEDTools
  bedtools intersect -v -a $postDedupMAPQ30BAM -b $BlackFilterPath > \${basename}_blackfiltered_${refgenome}.bam
  bedtools intersect -v -a $SmallFragBAM -b $BlackFilterPath > \${basename}.low120_blackfiltered_${refgenome}.bam
  bedtools intersect -v -a $LargeFragBAM -b $BlackFilterPath > \${basename}.high120_blackfiltered_${refgenome}.bam

  module purge
  module load SAMtools
  # Index last filtered file
  samtools index \${basename}_blackfiltered_${refgenome}.bam
  samtools index \${basename}.low120_blackfiltered_${refgenome}.bam
  samtools index \${basename}.high120_blackfiltered_${refgenome}.bam

  """
}

process RPKM_normalizing {
  publishDir "${params.resultsDir}/bwfiles_normRPKM" , mode: 'move', overwrite: true
  label 'big_task'

  input:
  path BlackfilteredBAMS
  val BlackFilterBaseName
  path BlackfilteredBAIS
  tuple val(refgenome),val(publish19),val(publish38)

  output:
  path "*.bw"

  script:
  """
  basename=$BlackFilterBaseName
  
  module purge
  module load deepTools
  bamCoverage -p $big_task_cpus -b \${basename}_blackfiltered_${refgenome}.bam --binSize 10 --normalizeUsing RPKM -of bigwig -o \${basename}_blackfiltered_${refgenome}.bw
  bamCoverage -p $big_task_cpus -b \${basename}.low120_blackfiltered_${refgenome}.bam --binSize 10 --normalizeUsing RPKM -of bigwig -o \${basename}.low120_blackfiltered_${refgenome}.bw
  bamCoverage -p $big_task_cpus -b \${basename}.high120_blackfiltered_${refgenome}.bam --binSize 10 --normalizeUsing RPKM -of bigwig -o \${basename}.high120_blackfiltered_${refgenome}.bw

  """
}

process IGV {
  publishDir "${params.resultsDir}/IGVfiles" , mode: 'move', overwrite: true
  label 'small_task'

  input:
  path BlackfilteredBAMS
  val BlackFilterBaseName
  path BlackfilteredBAIS
  tuple val(refgenome),val(publish19),val(publish38)

  output:
  path "*.tdf"
  
  script:
  IGVGenome=params.IGVGenome[task.process]
  """
  basename=$BlackFilterBaseName
  
  module purge
  module load IGV

  igvtools count \${basename}_blackfiltered_${refgenome}.bam \${basename}_blackfiltered_${refgenome}.tdf $IGVGenome

  igvtools count \${basename}.low120_blackfiltered_${refgenome}.bam \${basename}.low120_blackfiltered_${refgenome}.tdf $IGVGenome

  igvtools count \${basename}.high120_blackfiltered_${refgenome}.bam \${basename}.high120_blackfiltered_${refgenome}.tdf $IGVGenome
  """
}

process PeakCalling {
  publishDir "${params.resultsDir}/Macs2_peakcalling/${Target}_vs_${Control}.MAPQ30" , mode: 'copy', overwrite: true
  label 'big_task'

  input:
  path BlackfilteredBAMS
  tuple val(Control),val(Target),val(peakType)
  path BlackfilteredBAIS
  tuple val(refgenome),val(publish19),val(publish38)

  output:
  path "diff_peaks*"
  tuple env(peakDirName), path ("*.MAPQ30_keepdupauto_${refgenome}_chr_peaks.narrowPeak"), path ("*.MAPQ30.low120_keepdupauto_${refgenome}_chr_peaks.narrowPeak"), emit: PeakTuple, optional: true

  script:
  if (peakType == "narrowPeak")
  
  """
  module purge
  module load MACS2

  bashcontrol="${Control}*.MAPQ30_blackfiltered_${refgenome}.bam"
  smallfragcontrol="${Control}*.MAPQ30.low120_blackfiltered_${refgenome}.bam"

  bashtarget="${Target}*.MAPQ30_blackfiltered_${refgenome}.bam"
  smallfragtarget="${Target}*.MAPQ30.low120_blackfiltered_${refgenome}.bam"
  NumIDprecursor=$Target
  NumIDstep2=\${NumIDprecursor%%_*}
  NumID=\${NumIDstep2##CR}
  
  #MAPQ30
  macs2 callpeak -t \${bashtarget} --format BAMPE -c \${bashcontrol} --outdir . \
   -n diff_peaks_${Target}_vs_${Control}.MAPQ30_keepdupauto_${refgenome} -q 0.05 -g hs --bdg --keep-dup auto
  #Low120
  macs2 callpeak -t \${smallfragtarget} --format BAMPE -c \${smallfragcontrol} --outdir . \
   -n diff_peaks_${Target}_vs_${Control}.MAPQ30.low120_keepdupauto_${refgenome} -q 0.05 -g hs --bdg --keep-dup auto

  #filter MT peaks MAPQ30

  awk '\$NumID !~ /MT/' diff_peaks_${Target}_vs_${Control}.MAPQ30_keepdupauto_${refgenome}_peaks.narrowPeak  > diff_peaks_${Target}_vs_${Control}.MAPQ30_keepdupauto_${refgenome}_noMT_peaks.narrowPeak
  awk '\$NumID !~ /MT/' diff_peaks_${Target}_vs_${Control}.MAPQ30_keepdupauto_${refgenome}_summits.bed > diff_peaks_${Target}_vs_${Control}.MAPQ30_keepdupauto_${refgenome}_noMT_summits.bed

  sed 's/^/chr/' diff_peaks_${Target}_vs_${Control}.MAPQ30_keepdupauto_${refgenome}_peaks.narrowPeak >  diff_peaks_${Target}_vs_${Control}.MAPQ30_keepdupauto_${refgenome}_chr_peaks.narrowPeak

  #filter MT peaks low120
  
  awk '\$NumID !~ /MT/' diff_peaks_${Target}_vs_${Control}.MAPQ30.low120_keepdupauto_${refgenome}_peaks.narrowPeak  > diff_peaks_${Target}_vs_${Control}.MAPQ30.low120_keepdupauto_${refgenome}_noMT_peaks.narrowPeak
  awk '\$NumID !~ /MT/' diff_peaks_${Target}_vs_${Control}.MAPQ30.low120_keepdupauto_${refgenome}_summits.bed > diff_peaks_${Target}_vs_${Control}.MAPQ30.low120_keepdupauto_${refgenome}_noMT_summits.bed

  sed 's/^/chr/' diff_peaks_${Target}_vs_${Control}.MAPQ30.low120_keepdupauto_${refgenome}_peaks.narrowPeak >  diff_peaks_${Target}_vs_${Control}.MAPQ30.low120_keepdupauto_${refgenome}_chr_peaks.narrowPeak

  
  peakDirName="${Target}_vs_${Control}.MAPQ30"

  """
  else if (peakType == "broadPeak")
  """
  module purge
  module load MACS2

  bashcontrol="${Control}*.MAPQ30_blackfiltered_${refgenome}.bam"
  largefragcontrol="${Control}*.MAPQ30.high120_blackfiltered_${refgenome}.bam"

  bashtarget="${Target}*.MAPQ30_blackfiltered_${refgenome}.bam"
  largefragtarget="${Target}*.MAPQ30.high120_blackfiltered_${refgenome}.bam"
  NumIDprecursor=$Target
  NumIDstep2=\${NumIDprecursor%%_*}
  NumID=\${NumIDstep2##CR}

  #MAPQ30
  macs2 callpeak -t \${bashtarget} --format BAMPE -c \${bashcontrol} --outdir . \
   -n diff_peaks_${Target}_vs_${Control}.MAPQ30_keepdupauto_${refgenome} -g hs --keep-dup auto --bdg --broad --broad-cutoff 0.05

  #High120
  macs2 callpeak -t \${largefragtarget} --format BAMPE -c \${largefragcontrol} --outdir . \
   -n diff_peaks_${Target}_vs_${Control}.MAPQ30.high120_keepdupauto_${refgenome} -g hs --keep-dup auto --bdg --broad --broad-cutoff 0.05

  
  #filter MT peaks 

  awk '\$NumID !~ /MT/' diff_peaks_${Target}_vs_${Control}.MAPQ30_keepdupauto_${refgenome}_peaks.broadPeak  > diff_peaks_${Target}_vs_${Control}.MAPQ30_keepdupauto_${refgenome}_noMT_peaks.broadPeak

  awk '\$NumID !~ /MT/' diff_peaks_${Target}_vs_${Control}.MAPQ30.high120_keepdupauto_${refgenome}_peaks.broadPeak  > diff_peaks_${Target}_vs_${Control}.MAPQ30.high120_keepdupauto_${refgenome}_noMT_peaks.broadPeak

  peakDirName="${Target}_vs_${Control}.MAPQ30"
  """
  else
  """
  echo Something is wrong with one of the following: Input:$Control Target:$Target Peakcalling:$peakType please confirm that none of these fields are emtpy.
  """
}

process Homer_findMotif {
  publishDir "${params.resultsDir}/HomerMotifs/" , mode: 'move', overwrite: true
  pod env: 'PATH', value: '$VSC_DATA_VO/PPOL/resources/repos/homer/bin:$PATH'
  errorStrategy 'finish'
  label 'big_task'

  input:
  tuple val (peakDirName), path (NormalPeakFile), path (LowPeakFile)
  tuple val(refgenome),val(publish19),val(publish38)
  
  output:
  path "${peakDirName}*"

  script:
  if (NormalPeakFile.size() > 0 && LowPeakFile.size() > 0)
  """
  
  module purge
  module load Python
  
  #MAPQ30
  mkdir ${peakDirName}_keepdupauto_size${HomerSize}_unmasked_motifoutput_${refgenome}_MAPQ30
  findMotifsGenome.pl $NormalPeakFile $refgenome ./${peakDirName}_keepdupauto_size${HomerSize}_unmasked_motifoutput_${refgenome}_MAPQ30 -size $HomerSize -p $big_task_cpus 2> ./${peakDirName}_keepdupauto_size${HomerSize}_unmasked_motifoutput_${refgenome}_MAPQ30/HomerErrorOut.txt
  
  #MAPQ30.low$FragSize
  mkdir ${peakDirName}_keepdupauto_size${HomerSize}_unmasked_motifoutput_${refgenome}_low$FragSize
  findMotifsGenome.pl $LowPeakFile $refgenome ./${peakDirName}_keepdupauto_size${HomerSize}_unmasked_motifoutput_${refgenome}_low$FragSize -size $HomerSize -p $big_task_cpus 2> ./${peakDirName}_keepdupauto_size${HomerSize}_unmasked_motifoutput_${refgenome}_low$FragSize/HomerErrorOut.txt
  
  rm -f */*.tmp
  """
  else if (NormalPeakFile.size() > 0)
  """
  
  module purge
  module load Python
  
  #MAPQ30
  mkdir ${peakDirName}_keepdupauto_size${HomerSize}_unmasked_motifoutput_${refgenome}_MAPQ30
  findMotifsGenome.pl $NormalPeakFile $refgenome ./${peakDirName}_keepdupauto_size${HomerSize}_unmasked_motifoutput_${refgenome}_MAPQ30 -size $HomerSize -p $big_task_cpus 2> ./${peakDirName}_keepdupauto_size${HomerSize}_unmasked_motifoutput_${refgenome}_MAPQ30/HomerErrorOut.txt

  #MAPQ30.low$FragSize
  mkdir ${peakDirName}_keepdupauto_size${HomerSize}_unmasked_motifoutput_${refgenome}_low$FragSize
  echo "$LowPeakFile contained 0kb of data meaning no peaks were found during peakcalling" > ${peakDirName}_keepdupauto_size${HomerSize}_unmasked_motifoutput_${refgenome}_low$FragSize/pipeline.err
  
  rm -f */*.tmp
  """

  else if (LowPeakFile.size() > 0)
  """
  
  module purge
  module load Python
  
  #MAPQ30.low$FragSize
  mkdir ${peakDirName}_keepdupauto_size${HomerSize}_unmasked_motifoutput_${refgenome}_low$FragSize
  findMotifsGenome.pl $LowPeakFile $refgenome ./${peakDirName}_keepdupauto_size${HomerSize}_unmasked_motifoutput_${refgenome}_low$FragSize -size $HomerSize -p $big_task_cpus 2> ./${peakDirName}_keepdupauto_size${HomerSize}_unmasked_motifoutput_${refgenome}_low$FragSize/HomerErrorOut.txt

  #MAPQ30
  mkdir mkdir ${peakDirName}_keepdupauto_size${HomerSize}_unmasked_motifoutput_${refgenome}_MAPQ30
  echo "$NormalPeakFile contained 0kb of data meaning no peaks were found during peakcalling" > mkdir ${peakDirName}_keepdupauto_size${HomerSize}_unmasked_motifoutput_${refgenome}_MAPQ30/pipeline.err
  
  rm -f */*.tmp
  """

  else 
  """
  mkdir mkdir ${peakDirName}_keepdupauto_size${HomerSize}_unmasked_motifoutput_${refgenome}_MAPQ30
  echo "$NormalPeakFile contained 0kb of data meaning no peaks were found during peakcalling" > mkdir ${peakDirName}_keepdupauto_size${HomerSize}_unmasked_motifoutput_${refgenome}_MAPQ30/pipeline.err
  
  mkdir ${peakDirName}_keepdupauto_size${HomerSize}_unmasked_motifoutput_${refgenome}_low$FragSize
  echo "$LowPeakFile contained 0kb of data meaning no peaks were found during peakcalling" > ${peakDirName}_keepdupauto_size${HomerSize}_unmasked_motifoutput_${refgenome}_low$FragSize/pipeline.err
  
  """
}

process Homer_annotatePeaks {
  publishDir "${params.resultsDir}/Macs2_peakcalling/${peakDirName}" , mode: 'move', overwrite: true
  pod env: 'PATH', value: '$VSC_DATA_VO/PPOL/resources/repos/homer/bin:$PATH'
  errorStrategy 'finish'
  label 'big_task'

  input:
  tuple val (peakDirName), path (NormalPeakFile), path (LowPeakFile)
  tuple val(refgenome),val(publish19),val(publish38)
  
  output:
  path "*_genomeontology", optional: true
  path "*_{homerenhancer,homerpromoter,annotation,err}.txt"

  script:
  if (NormalPeakFile.size() > 0 && LowPeakFile.size() > 0)
  """
  inputname=$NormalPeakFile
  outfilename=\${inputname%%_peaks.narrowPeak}
  
  PATH=\$PATH:/user/data/gent/gvo000/gvo00027/PPOL/resources/repos/homer/.//bin/

  module purge 
  module load Perl
  
  #MAPQ30
  annotatePeaks.pl $NormalPeakFile $refgenome -genomeOntology \${outfilename}_genomeontology/ > \${outfilename}_annotation.txt 2> \${outfilename}_err.txt

  grep promoter-TSS \${outfilename}_annotation.txt > \${outfilename}_homerpromoter.txt || true

  grep -v promoter-TSS \${outfilename}_annotation.txt > \${outfilename}_homerenhancer.txt || true
  
  sed '1d' \${outfilename}_homerenhancer.txt > \${outfilename}_homerenhancer.txt

  #MAPQ30_low$FragSize
  inputname=$LowPeakFile
  outfilename=\${inputname%%_peaks.narrowPeak}

  annotatePeaks.pl $LowPeakFile $refgenome -genomeOntology \${outfilename}_genomeontology/ > \${outfilename}_annotation.txt 2> \${outfilename}_err.txt

  grep promoter-TSS \${outfilename}_annotation.txt > \${outfilename}_homerpromoter.txt || true

  grep -v promoter-TSS \${outfilename}_annotation.txt > \${outfilename}_homerenhancer.txt || true
  
  sed '1d' \${outfilename}_homerenhancer.txt > \${outfilename}_homerenhancer.txt
  """
  else if (NormalPeakFile.size() > 0)
  """
  inputname=$NormalPeakFile
  outfilename=\${inputname%%_peaks.narrowPeak}
  
  PATH=\$PATH:/user/data/gent/gvo000/gvo00027/PPOL/resources/repos/homer/.//bin/

  module purge 
  module load Perl
  
  #MAPQ30
  annotatePeaks.pl $NormalPeakFile $refgenome -genomeOntology \${outfilename}_genomeontology/ > \${outfilename}_annotation.txt 2> \${outfilename}_err.txt

  grep promoter-TSS \${outfilename}_annotation.txt > \${outfilename}_homerpromoter.txt || true

  grep -v promoter-TSS \${outfilename}_annotation.txt > \${outfilename}_homerenhancer.txt || true
  
  sed '1d' \${outfilename}_homerenhancer.txt > \${outfilename}_homerenhancer.txt

  #MAPQ30_low$FragSize
  
  inputname=$LowPeakFile
  outfilename=\${inputname%%_peaks.narrowPeak}

  echo "The size of \$inputname was 0 bytes indicating that no peaks were detected during peakcalling. As a result no Homer-analysis was performed." > \${outfilename}_err.txt
  
  """
  else if (LowPeakFile.size() > 0)
  """
  inputname=$LowPeakFile
  outfilename=\${inputname%%_peaks.narrowPeak}
  
  PATH=\$PATH:/user/data/gent/gvo000/gvo00027/PPOL/resources/repos/homer/.//bin/

  module purge 
  module load Perl
  
  #MAPQ30_low$FragSize
  annotatePeaks.pl $LowPeakFile $refgenome -genomeOntology \${outfilename}_genomeontology/ > \${outfilename}_annotation.txt 2> \${outfilename}_err.txt

  grep promoter-TSS \${outfilename}_annotation.txt > \${outfilename}_homerpromoter.txt || true

  grep -v promoter-TSS \${outfilename}_annotation.txt > \${outfilename}_homerenhancer.txt || true
  
  sed '1d' \${outfilename}_homerenhancer.txt > \${outfilename}_homerenhancer.txt

  #MAPQ30
  
  inputname=$NormalPeakFile
  outfilename=\${inputname%%_peaks.narrowPeak}

  echo "The size of \$inputname was 0 bytes indicating that no peaks were detected during peakcalling. As a result no Homer-analysis was performed." > \${outfilename}_err.txt
  
  """
  else 
  """
  #MAPQ30
  
  inputname=$NormalPeakFile
  outfilename=\${inputname%%_peaks.narrowPeak}

  echo "The size of \$inputname was 0 bytes indicating that no peaks were detected during peakcalling. As a result no Homer-analysis was performed." > \${outfilename}_err.txt
  
  #MAPQ30_low$FragSize
  
  inputname=$LowPeakFile
  outfilename=\${inputname%%_peaks.narrowPeak}

  echo "The size of \$inputname was 0 bytes indicating that no peaks were detected during peakcalling. As a result no Homer-analysis was performed." > \${outfilename}_err.txt
  
  """

}

process RunName{
  label 'small_task'
  input:
  val (runname)
  tuple val(refgenome),val(publish19),val(publish38)
  output:
  env runname

  script:
  """
  nameprecursor=$runname
  intermediate=\${nameprecursor##*/}
  runname=\${intermediate}_$refgenome
  """
}

process summary_MultiQC {
  publishDir "${params.resultsDir}/run_MultiQC/$RunName" , mode: 'move', overwrite: true
  label 'big_task'

  input:
  path BamCollection
  path BaiCollection
  val RunName
  tuple val(refgenome),val(publish19),val(publish38)

  output:
  path "*_multiQC.html"
  path "*_multiQC_data"
  path "*.log"
  path "*fastqc.html"

  script:
  """
    module purge
    module load SAMtools

    rm *{low,high}$FragSize*.bam

    for bamfile in *.bam; do
    basename=\${bamfile%%.bam}
    samtools stats \$bamfile > \${basename}.samtoolstats.Hs.log
    samtools flagstat \$bamfile > \${basename}.flagstat.log
    done;

    module purge
    module load FastQC

    fastqc -t $big_task_cpus $BamCollection

    module purge
    module load MultiQC

    runname=$RunName
    
    multiqc -f . -n \${runname}_multiQC
  """
}

workflow trim {


  main:
  def filepairs = Channel.fromFilePairs("$params.inputDir/*{R1,R2}*.fastq.gz", checkIfExists:true) 
  def SampleTuple = Channel.fromPath(params.SampleSheet, checkIfExists:true) | splitCsv(header: true, strip: true ,sep:';') | map { row-> tuple(row.DNA.replaceAll(/ /,""), row.CR.replaceAll(/ /,"")) }

    File_sheet_match(filepairs)
  | join(SampleTuple)
  | MoveToDataDir 
  | TrimFiles

  emit:
  TrimmedFile1 = TrimFiles.out.trimmedfile1
  TrimmedFile2 = TrimFiles.out.trimmedfile2
  TrimReport1 = TrimFiles.out.trimreport1
  TrimReport2 = TrimFiles.out.trimreport2
  workingDir = TrimFiles.out.workingDir
  workingDirName = TrimFiles.out.workingDirName
}

workflow hg19 {
  take:
  trimmedfile1
  trimmedfile2
  trimreport1
  trimreport2
  workingDir
  workingDirName

  main:
    def PeakCallingSheet = Channel.fromPath(params.PeakcallingSheet, checkIfExists:true) | splitCsv(header:true, strip: true,sep:';') | map { row-> tuple(row.Input.replaceAll(/ /,""), row.Target.replaceAll(/ /,""),row.Peakcalling.replaceAll(/ /,"")) }
    def runname = Channel.fromPath(inputDir)
    
    MapFiles(trimmedfile1,trimmedfile2,trimreport1,trimreport2,workingDir,workingDirName,params.hg19Tuple)
  | MultiQC
    
    MultiQC.out.Dedup_percent
  | toFloat
  | merge (MultiQC.out.workingDirName) {a,b -> tuple(b,a)}
  | join (MapFiles.out.deduptuple)
  | set {Dedup_Input}
  
    Dedup(Dedup_Input,params.hg19Tuple)
    
    SplitFragments(Dedup.out.postDedupMAPQ30BAM,Dedup.out.workingDirName,Dedup.out.dedupBAI,params.hg19Tuple)
    BlackFiltering(SplitFragments.out.SmallFragBAM,SplitFragments.out.LargeFragBAM,SplitFragments.out.workingDirName,SplitFragments.out.postDedupMAPQ30BAM,SplitFragments.out.splitFragBAI,params.hg19Tuple)
  
    RPKM_normalizing(BlackFiltering.out.BlackfilteredBAMS,BlackFiltering.out.BlackFilterBaseName,BlackFiltering.out.BlackfilteredBAIS,params.hg19Tuple)
    IGV(BlackFiltering.out.BlackfilteredBAMS,BlackFiltering.out.BlackFilterBaseName,BlackFiltering.out.BlackfilteredBAIS,params.hg19Tuple)
    
    BlackFiltering.out.BlackfilteredBAMS
  |  collect
  |  set{BamCollection}
   BlackFiltering.out.BlackfilteredBAIS
  |  collect
  |  set{BaiCollection}
    PeakCalling(BamCollection,PeakCallingSheet,BaiCollection,params.hg19Tuple)

    Homer_findMotif(PeakCalling.out.PeakTuple,params.hg19Tuple)
    Homer_annotatePeaks(PeakCalling.out.PeakTuple,params.hg19Tuple)

    RunName(runname,params.hg19Tuple)

    summary_MultiQC(BamCollection,BaiCollection,RunName.out,params.hg19Tuple)

}

workflow hg38 {
  take:
  trimmedfile1
  trimmedfile2
  trimreport1
  trimreport2
  workingDir
  workingDirName

  main:
  def PeakCallingSheet = Channel.fromPath(params.PeakcallingSheet, checkIfExists:true) | splitCsv(header:true, strip: true,sep:';') | map { row-> tuple(row.Input.replaceAll(/ /,""), row.Target.replaceAll(/ /,""),row.Peakcalling.replaceAll(/ /,"")) }
    def runname = Channel.fromPath(inputDir)
    
    MapFiles(trimmedfile1,trimmedfile2,trimreport1,trimreport2,workingDir,workingDirName,params.hg38Tuple)
  | MultiQC
    
    MultiQC.out.Dedup_percent
  | toFloat
  | merge (MultiQC.out.workingDirName) {a,b -> tuple(b,a)}
  | join (MapFiles.out.deduptuple)
  | set {Dedup_Input}
  
    Dedup(Dedup_Input,params.hg38Tuple)
    
    SplitFragments(Dedup.out.postDedupMAPQ30BAM,Dedup.out.workingDirName,Dedup.out.dedupBAI,params.hg38Tuple)
    BlackFiltering(SplitFragments.out.SmallFragBAM,SplitFragments.out.LargeFragBAM,SplitFragments.out.workingDirName,SplitFragments.out.postDedupMAPQ30BAM,SplitFragments.out.splitFragBAI,params.hg38Tuple)
  
    RPKM_normalizing(BlackFiltering.out.BlackfilteredBAMS,BlackFiltering.out.BlackFilterBaseName,BlackFiltering.out.BlackfilteredBAIS,params.hg38Tuple)
    IGV(BlackFiltering.out.BlackfilteredBAMS,BlackFiltering.out.BlackFilterBaseName,BlackFiltering.out.BlackfilteredBAIS,params.hg38Tuple)
    
    BlackFiltering.out.BlackfilteredBAMS
  |  collect
  |  set{BamCollection}
   BlackFiltering.out.BlackfilteredBAIS
  |  collect
  |  set{BaiCollection}
    PeakCalling(BamCollection,PeakCallingSheet,BaiCollection,params.hg38Tuple)

    Homer_findMotif(PeakCalling.out.PeakTuple,params.hg38Tuple)
    Homer_annotatePeaks(PeakCalling.out.PeakTuple,params.hg38Tuple)

    RunName(runname,params.hg38Tuple)

    summary_MultiQC(BamCollection,BaiCollection,RunName.out,params.hg38Tuple)
}

workflow {

  trim()
  if (pipelineMode == "both"){
  hg19(trim.out)
  hg38(trim.out)
  }
  else if (pipelineMode == "hg19"){
  hg19(trim.out)
  }
  else if (pipelineMode == "hg38"){
  hg38(trim.out)
  }
}
