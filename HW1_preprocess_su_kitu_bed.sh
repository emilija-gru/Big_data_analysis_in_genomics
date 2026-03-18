# GENE BODY COVERAGE

# atsisiunčiam bed failą
cd ~/HW1/references
wget https://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens/hg38_GENCODE_V47.bed.gz
gunzip hg38_GENCODE_V47.bed.gz

# ((((((((((((((()))))))))))))))
conda activate rseqc_env
# ((((((((((((((()))))))))))))))

# leidžiamt iš mapping
cd ~/HW1/raw_data/fastq_trimmed/mapping
geneBody_coverage.py -r ~/HW1/references/hg38_GENCODE_V47.bed -i *_sorted.bam -o ~/HW1/gene_body_cov/gene_body_cov_2


# INNER DISTANCE - gavosi tas pats

inner_distance.py -r ~/HW1/references/hg38_GENCODE_V47.bed \
-i ~/HW1/raw_data/fastq_trimmed/mapping/*_sorted.bam \
-o ~/HW1/inner_distance/inner_distance_2

# ANNOTATED JUNCTIONS - gavosi kitaip

junction_annotation.py \
-i ~/HW1/raw_data/fastq_trimmed/mapping/*_sorted.bam \
-r ~/HW1/references/hg38_GENCODE_V47.bed \
-o ~/HW1/annotated_junctions/junction_annotation_2

# nustatinėju, ar biblioteka yra stranded ar unstranded - STRANDED

infer_experiment.py \
-i ~/HW1/raw_data/fastq_trimmed/mapping/SRR11647692_sorted.bam \
-r ~/HW1/references/hg38_GENCODE_V47.bed