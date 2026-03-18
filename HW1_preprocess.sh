# Failus parsisiunčiau su wget

# pasižiūriu referentinių chromosomų skaičių 

grep ">" GRCh38.primary_assembly.genome.fa | wc -l

# paleidžiam fastqc ant raw failų

fastqc -t 8 -o fastqc_raw *.fastq.gz 

# multiqc su raw failais
cd ~/HW1/raw_data/fastqc_raw
multiqc .
multiqc /home/genetics/HW1/raw_data/fastqc_raw/ -o /home/genetics/HW1/multiqc_report_raw/

# multiqc nekrovė plot'ų, tai su chatgpt pagalba bandžiau sutvarkyti
conda create -n mqc -c conda-forge -c bioconda multiqc
conda activate mqc

#Trimming

trim_galore -o fastq_trimmed -j 8 --quality 20 --length 20 --paired SRR11647692_1.fastq.gz SRR11647692_2.fastq.gz
trim_galore -o fastq_trimmed -j 8 --quality 20 --length 20 --paired SRR11647693_1.fastq.gz SRR11647693_2.fastq.gz
trim_galore -o fastq_trimmed -j 8 --quality 20 --length 20 --paired SRR11647694_1.fastq.gz SRR11647694_2.fastq.gz
trim_galore -o fastq_trimmed -j 8 --quality 20 --length 20 --paired SRR11647702_1.fastq.gz SRR11647702_2.fastq.gz
trim_galore -o fastq_trimmed -j 8 --quality 20 --length 20 --paired SRR11647703_1.fastq.gz SRR11647703_2.fastq.gz
trim_galore -o fastq_trimmed -j 8 --quality 20 --length 20 --paired SRR11647704_1.fastq.gz SRR11647704_2.fastq.gz
# -j 8 pagreitina procesą (dienos metu geriau naudoti 6, o vakare galima ir 8)

ls
cd HW1/raw_data/fastq_trimmed
mkdir fastqc_trimmed

# paleidžiam fastqc ir multiqc ant trimmed failų

fastqc -t 8 -o fastqc_trimmed *.fq.gz 
cd fastqc_trimmed
multiqc .
multiqc /home/genetics/HW1/raw_data/fastq_trimmed/fastqc_trimmed -o /home/genetics/HW1/multiqc_report_trimmed/


# pereinam į references vietoj raw data ir INDEKSUOJAM referentinį genomą
cd ~/HW1/references  
hisat2-build GRCh38.primary_assembly.genome.fa.gz indexed_ref # neveikia, nes zip failas buvo
gunzip GRCh38.primary_assembly.genome.fa.gz # unzipinam
hisat2-build -p 6 GRCh38.primary_assembly.genome.fa indexed_ref

# Mappinam

hisat2 -p 8 -q -x ~/HW1/references/indexed_ref -1 SRR11647692_1_val_1.fq.gz -2 SRR11647692_2_val_2.fq.gz -S mapping/SRR11647692.sam      
hisat2 -p 8 -q -x ~/HW1/references/indexed_ref -1 SRR11647693_1_val_1.fq.gz -2 SRR11647693_2_val_2.fq.gz -S mapping/SRR11647693.sam
hisat2 -p 8 -q -x ~/HW1/references/indexed_ref -1 SRR11647694_1_val_1.fq.gz -2 SRR11647694_2_val_2.fq.gz -S mapping/SRR11647694.sam
hisat2 -p 8 -q -x ~/HW1/references/indexed_ref -1 SRR11647702_1_val_1.fq.gz -2 SRR11647702_2_val_2.fq.gz -S mapping/SRR11647702.sam
hisat2 -p 8 -q -x ~/HW1/references/indexed_ref -1 SRR11647703_1_val_1.fq.gz -2 SRR11647703_2_val_2.fq.gz -S mapping/SRR11647703.sam
hisat2 -p 8 -q -x ~/HW1/references/indexed_ref -1 SRR11647704_1_val_1.fq.gz -2 SRR11647704_2_val_2.fq.gz -S mapping/SRR11647704.sam


# konvertuojam sam į bam

cd ~/HW1/raw_data/fastq_trimmed/mapping
samtools view -b SRR11647692.sam > SRR11647692.bam
samtools view -b SRR11647693.sam > SRR11647693.bam
samtools view -b SRR11647694.sam > SRR11647694.bam
samtools view -b SRR11647702.sam > SRR11647702.bam
samtools view -b SRR11647703.sam > SRR11647703.bam
samtools view -b SRR11647704.sam > SRR11647704.bam

# patikrinam, ar bam susikūrė

ls *.bam

# patikrinu, ar bam susikūrė teisingai. Jei matosi alignment eilutes → BAM sukurtas teisingai.

samtools view SRR11647692.bam | head
samtools view SRR11647693.bam | head
samtools view SRR11647694.bam | head
samtools view SRR11647702.bam | head
samtools view SRR11647703.bam | head
samtools view SRR11647704.bam | head

# ARBA patikrinam visus vienu metu ir jei nieko neišveda - visi BAM geri.

samtools quickcheck *.bam

# jei su bam viskas gerai, ištrinam sam

rm SRR11647692.sam
rm SRR11647693.sam
rm SRR11647694.sam
rm SRR11647702.sam
rm SRR11647703.sam
rm SRR11647704.sam

# išsortinam bam failus

samtools sort SRR11647692.bam -o SRR11647692_sorted.bam
samtools sort SRR11647693.bam -o SRR11647693_sorted.bam
samtools sort SRR11647694.bam -o SRR11647694_sorted.bam
samtools sort SRR11647702.bam -o SRR11647702_sorted.bam
samtools sort SRR11647703.bam -o SRR11647703_sorted.bam
samtools sort SRR11647704.bam -o SRR11647704_sorted.bam

# indeksuojam bam failus

samtools index SRR11647692_sorted.bam
samtools index SRR11647693_sorted.bam
samtools index SRR11647694_sorted.bam
samtools index SRR11647702_sorted.bam
samtools index SRR11647703_sorted.bam
samtools index SRR11647704_sorted.bam

# iš HW, ieškau mapping rate kiekvienam mėginiui ir paskui naudoju grafiko braižymui
cd ~/HW1/raw_data/fastq_trimmed/mapping
samtools flagstat SRR11647692_sorted.bam
samtools flagstat SRR11647693_sorted.bam
samtools flagstat SRR11647694_sorted.bam
samtools flagstat SRR11647702_sorted.bam
samtools flagstat SRR11647703_sorted.bam
samtools flagstat SRR11647704_sorted.bam


# Read Counting (Feature Quantification)

featureCounts -T 4 -p -t exon -g gene_id -s 0 -a ~/HW1/references/gencode.v49.primary_assembly.basic.annotation.gtf.gz -o counts.txt *_sorted.bam


# ALIGNMENT QUALITY ASSESSMENT

# DUPLICATES

cd ~/HW1
mkdir duplicates # iš HW1

# leidžiam iš mapping
cd ~/HW1/raw_data/fastq_trimmed/mapping

samtools collate -@ 4 -O -u SRR11647692_sorted.bam | samtools fixmate -@ 4 -m -u - - | samtools sort -@ 4 -u - | samtools markdup -@ 4 - ~/HW1/duplicates/SRR11647692_markdup.bam
samtools view -c -f 1024 ~/HW1/duplicates/SRR11647692_markdup.bam > ~/HW1/duplicates/SRR11647692_duplicates.txt

samtools collate -@ 4 -O -u SRR11647693_sorted.bam | samtools fixmate -@ 4 -m -u - - | samtools sort -@ 4 -u - | samtools markdup -@ 4 - ~/HW1/duplicates/SRR11647693_markdup.bam
samtools view -c -f 1024 ~/HW1/duplicates/SRR11647693_markdup.bam > ~/HW1/duplicates/SRR11647693_duplicates.txt

samtools collate -@ 4 -O -u SRR11647694_sorted.bam | samtools fixmate -@ 4 -m -u - - | samtools sort -@ 4 -u - | samtools markdup -@ 4 - ~/HW1/duplicates/SRR11647694_markdup.bam
samtools view -c -f 1024 ~/HW1/duplicates/SRR11647694_markdup.bam > ~/HW1/duplicates/SRR11647694_duplicates.txt

samtools collate -@ 4 -O -u SRR11647702_sorted.bam | samtools fixmate -@ 4 -m -u - - | samtools sort -@ 4 -u - | samtools markdup -@ 4 - ~/HW1/duplicates/SRR11647702_markdup.bam
samtools view -c -f 1024 ~/HW1/duplicates/SRR11647702_markdup.bam > ~/HW1/duplicates/SRR11647702_duplicates.txt

samtools collate -@ 4 -O -u SRR11647703_sorted.bam | samtools fixmate -@ 4 -m -u - - | samtools sort -@ 4 -u - | samtools markdup -@ 4 - ~/HW1/duplicates/SRR11647703_markdup.bam
samtools view -c -f 1024 ~/HW1/duplicates/SRR11647703_markdup.bam > ~/HW1/duplicates/SRR11647703_duplicates.txt

samtools collate -@ 4 -O -u SRR11647704_sorted.bam | samtools fixmate -@ 4 -m -u - - | samtools sort -@ 4 -u - | samtools markdup -@ 4 - ~/HW1/duplicates/SRR11647704_markdup.bam
samtools view -c -f 1024 ~/HW1/duplicates/SRR11647704_markdup.bam > ~/HW1/duplicates/SRR11647704_duplicates.txt


# COVERAGE

mkdir ~/HW1/coverage

# leidžiam iš mapping
cd ~/HW1/raw_data/fastq_trimmed/mapping

samtools coverage SRR11647692_sorted.bam > ~/HW1/coverage/SRR11647692_coverage.txt
samtools coverage SRR11647693_sorted.bam > ~/HW1/coverage/SRR11647693_coverage.txt
samtools coverage SRR11647694_sorted.bam > ~/HW1/coverage/SRR11647694_coverage.txt
samtools coverage SRR11647702_sorted.bam > ~/HW1/coverage/SRR11647702_coverage.txt
samtools coverage SRR11647703_sorted.bam > ~/HW1/coverage/SRR11647703_coverage.txt
samtools coverage SRR11647704_sorted.bam > ~/HW1/coverage/SRR11647704_coverage.txt


# GENE BODY COVERAGE

# atsisiunčiam bed failą
cd ~/HW1/reference
wget https://sourceforge.net/projects/rseqc/files/BED/hg38.HouseKeepingGenes.bed.gz
gunzip hg38.HouseKeepingGenes.bed.gz

# ((((((((((((((()))))))))))))))
# reikėjo įsidiegti RSeQC
conda create -n rseqc_env -c conda-forge -c bioconda python=3.10 rseqc
conda activate rseqc_env
geneBody_coverage.py --help
# ((((((((((((((()))))))))))))))

# leidžiamt iš mapping
cd ~/HW1/raw_data/fastq_trimmed/mapping
geneBody_coverage.py -r ~/HW1/references/hg38.HouseKeepingGenes.bed -i *_sorted.bam -o ~/HW1/gene_body_cov/gene_body_cov


# CORRELATION

mkdir ~/HW1/correlation
cd ~/HW1/raw_data/fastq_trimmed/mapping

#koreliacijos matrica

# leidžiu iš mapping
multiBamSummary bins --bamfiles *_sorted.bam -p 8 -o ~/HW1/correlation/out.npz

plotCorrelation -in ~/HW1/correlation/out.npz -c spearman -p heatmap -o ~/HW1/correlation/correlation_heatmap_deeptools.png --plotNumbers

# PCA

mkdir ~/HW1/PCA

plotPCA -in ~/HW1/correlation/out.npz -o ~/HW1/PCA/PCA_plot.png 




# inner distance, clipping profile ir annotated junctions plots
# prieš tai reikia aktyvuoti RSeQC

mkdir ~/HW1/inner_distance
mkdir ~/HW1/clipping_profile
mkdir ~/HW1/annotated_junctions

# inner distance
inner_distance.py -r ~/HW1/references/hg38.HouseKeepingGenes.bed \
-i ~/HW1/raw_data/fastq_trimmed/mapping/*_sorted.bam \
-o ~/HW1/inner_distance/inner_distance

# clipping profile

bam_stat.py -i ~/HW1/raw_data/fastq_trimmed/mapping/*_sorted.bam \
> ~/HW1/clipping_profile/bam_statistics.txt

for bam in ~/HW1/raw_data/fastq_trimmed/mapping/*_sorted.bam
do
    name=$(basename "$bam" _sorted.bam)
    clipping_profile.py -i "$bam" -o ~/HW1/clipping_profile/"$name" -s PE
done

cd ~/HW1/clipping_profile
# sujungia į vieną pdf
pdfunite *R1.pdf clipping_profile_all_samples.pdf
# 6 grafikus sujungia į vieną grafiką
montage *R1.pdf -tile 3x2 -geometry +10+10 clipping_profile_combined.png

# annotated junctions

junction_annotation.py \
-i ~/HW1/raw_data/fastq_trimmed/mapping/*_sorted.bam \
-r ~/HW1/references/hg38.HouseKeepingGenes.bed \
-o ~/HW1/annotated_junctions/junction_annotation

# nustatinėju, ar biblioteka yra stranded ar unstranded

infer_experiment.py \
-i ~/HW1/raw_data/fastq_trimmed/mapping/SRR11647692_sorted.bam \
-r ~/HW1/references/hg38.HouseKeepingGenes.bed





# sudedu visus gautus grafikus į vieną direktoriją

mkdir grafikai_gauti_per_R
mv *.png grafikai_gauti_per_R/