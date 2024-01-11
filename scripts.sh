cat downloadlist.txt | while read line; do name=`echo $line | awk '{print $1}'| awk -F "?" '{print $1}' | awk -F "-" '{print $NF}'`; mwget -c 3 -n 10 "$line" -d ./ -f $name ; done

echo ZPY02$'\t'brain$'\n'ZPY03$'\t'cerebelllum > sample_list.txt

cat sample_list.txt | while read line; do file=`echo $line| awk '{print $1}'`; sample=`echo $line| awk '{print $2}'`; ln -s 20230815-YueGangAo-${file}/01.rawFq/00.mergeRawFq/5hmc/5hmc_raw_1.fq.gz ./${sample}_5hmc_R1.fastq.gz; ln -s 20230815-YueGangAo-${file}/01.rawFq/00.mergeRawFq/5hmc/5hmc_raw_2.fq.gz ./${sample}_5hmc_R2.fastq.gz; done

ln -s ../../scripts/sprite-pipeline-master/configuration_file.txt

source activate h5mC_env
for i in *_R1.fastq.gz; do echo "java -jar /picb/rnomics4/rotation/fuzhican/download/Trimmomatic-0.38/trimmomatic-0.38.jar PE  -threads 15 -phred33 ${i} ${i%_R1.fastq.gz}_R2.fastq.gz -baseout ./${i%_R1.fastq.gz}.fq.gz ILLUMINACLIP:/picb/rnomics4/rotation/fuzhican/download/Trimmomatic-0.38/adapters/TruSeq3-PE-2.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50 > ${i%.fastq.gz}.trim.log 2>&1" ; done > 01_script.sh

ln -s /data/rnomics10/peihong/2022_5hmC/20221102_testData/0_fastq/Split-Pool-jx_3barcode_rev.fa ./

for i in *_R1.fastq.gz; do echo "java -jar ../../scripts/sprite-pipeline-master/scripts/java/BarcodeIdentification_v1.2.0.jar --input1 ${i%_R1.fastq.gz}_1P.fq.gz --input2 ${i%_R1.fastq.gz}_2P.fq.gz --output1 ${i%_R1.fastq.gz}_1P.barcode.fq.gz --output2 ${i%_R1.fastq.gz}_2P.barcode.fq.gz --config configuration_file.txt > ${i%_R1.fastq.gz}_BarcodeIdentif.log 2>&1; ../../scripts/full_barcodes.sh ${i%_R1.fastq.gz}_1P.barcode.fq.gz" ; done > 02_BarcodeIdentification.sh

for i in *_R1.fastq.gz ; do echo ''' cutadapt -j 15 -m 50 -a "file:Split-Pool-jx_3barcode_rev.fa" ''' " -o ${i%*_R1.fastq.gz}_1P.barcode_full_clean.fq.gz ${i%*_R1.fastq.gz}_1P.barcode_full.fq.gz > ${i%*_R1.fastq.gz}_1P.barcode_full_clean.log 2>&1" ; done > 04_remove_3pribarcode.sh


echo brain_5hmc$'\t'../../20230712_4_mousebrain/0_fastq/brain_5hmc_1P.barcode_full_clean.fq.gz,../../20230712_3_mousebrain/0_fastq/brain_5hmc_1P.barcode_full_clean.fq.gz,../../20230712_2_mousebrain/0_fastq/brain_5hmc_1P.barcode_full_clean.fq.gz,../../20230712_mousebrain/0_fastq/brain_5hmc_1P.barcode_full_clean.fq.gz$'\n'cerebelllum_5hmc$'\t'../../20230712_5_mousebrain/0_fastq/cerebelllum_5hmc_1P.barcode_full_clean.fq.gz,../../20230712_4_mousebrain/0_fastq/cerebelllum_5hmc_1P.barcode_full_clean.fq.gz,../../20230712_2_mousebrain/0_fastq/cerebelllum_5hmc_1P.barcode_full_clean.fq.gz,../../20230712_3_mousebrain/0_fastq/cerebelllum_5hmc_1P.barcode_full_clean.fq.gz,../../20230712_mousebrain/0_fastq/cerebelllum-2_5hmc_1P.barcode_full_clean.fq.gz,../../20230712_mousebrain/0_fastq/cerebelllum_5hmc_1P.barcode_full_clean.fq.gz > sample_list2  ## 需要个性化配置

cat sample_list2 | while read line; do name=`echo $line| awk '{print $1}'`; file=`echo $line| awk '{print $2}'`; echo "bowtie2  --trim5 11 -q -p 15  -x /picb/rnomics4/peihong/mm10_index/genome/mm10_all -U $file -S ../1_bowtie2/${name}.sam  > ../1_bowtie2/${name}.log 2>&1; samtools view -@ 10 -bS ../1_bowtie2/${name}.sam | samtools sort -@ 10 > ../1_bowtie2/${name}.bam;   samtools view -@ 10 -F4 -h ../1_bowtie2/${name}.bam | " ''' grep -v "XS:i:"| awk -v OFS="\t" '\''BEGIN{n=1}{if($1 ~ /^@/) print $0; else {split($1, a, "::"); $1=a[1]; $NF=$NF"\tCB:Z:"a[2]"\tUB:Z:"n; n=n+1; print $0}}'\'' ''' "| samtools view -@ 10 -bS > ../1_bowtie2/${name}.tag.bam;  samtools index ../1_bowtie2/${name}.tag.bam;  ../../scripts/removeDup.py ../1_bowtie2/${name}.tag.bam;  samtools view ../1_bowtie2/${name}.tag_rmDup.bam | " ''' awk '\''{print $(NF-1)}'\'' ''' "| sort | uniq -c | sort -k1nr > ../1_bowtie2/${name}.tag_rmDup.cells "; done > 05_bowtie2.sh

