# Bulk-RNAseq---miRNA

# I. Download data and tools 
# 1. Download raw data
```bash
mkdir raw
cd /path/to/raw
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR129/006/SRR12952906/SRR12952906.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR129/004/SRR12952904/SRR12952904.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR129/000/SRR12952900/SRR12952900.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR129/003/SRR12952903/SRR12952903.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR129/002/SRR12952902/SRR12952902.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR129/098/SRR12952898/SRR12952898.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR129/099/SRR12952899/SRR12952899.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR129/001/SRR12952901/SRR12952901.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR129/005/SRR12952905/SRR12952905.fastq.gz
```
# 2. Download reference
```bash
## Reference genome 
mkdir reference
p_ref='/path/to/reference'
cd $p_ref
curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCA_012934285.1/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT&filename=GCA_012934285.1.zip" -H "Accept: application/zip"
unzip GCA_012934285.1.zip 
sed -i '/^>/ s/ .*//' /path/to/ncbi_datasets/ncbi_dataset/data/GCA_012934285.1/GCA_012934285.1_ASM1293428v1_genomic.fna  # remove space in fasta sequence name

## Mature miRNA from miRbase 
wget https://www.mirbase.org/ftp/CURRENT/mature.fa.gz -O $p_ref/mature.fa.gz
gunzip /path/to/reference/mature.fa.gz
sed -i '/^>/ s/ .*//' /path/to/reference/mature.fa.gz


## Hairpin
wget https://www.mirbase.org/ftp/CURRENT/hairpin.fa.gz -O $p_ref/hairpin.fa.gz
gunzip $p_ref/hairpin.fa.gz
sed -i '/^>/ s/ .*//' /path/to/reference/hairpin.fa.gz
sed -e '/^[^>]/s/[^acgtunACGTUN]/N/g' $p_ref/hairpin.fa.gz > $p_ref/hairpin_cleanup.fa

```
# 3. Install tools for practice
## Fastqc
```bash
mkdir tools
cd /path/to/tools
## Download FASTQC
### Install & unzip
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip
unzip fastqc_v0.12.1.zip
### Write at the end of the .bashrc file. This command will let you export the path to FastQC, and
execute it everywhere
nano ~/.bashrc
export PATH='path/to/FastQC/':$PATH
### For example:
export PATH='/path/to/tools/FastQC':$PATH
source ~/.bashrc
### Try it by running
fastqc
```
## Install multiqc
```bash
pip install multiqc 
### or using conda to install 
conda install multiqc  # https://multiqc.info/
```

## Install Trimmomatic
```bash
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
unzip Trimmomatic-0.39.zip
```

```bash
wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 -O samtools.tar.bz2
tar -xjvf samtools.tar.bz2
cd samtools-{version}
make
sudo make prefix=/usr/local/bin install
```
## Install bowtie2
```bash
# Using conda 
conda install -c bioconda bowtie2
```

## Install miRDeep2
```bash
git clone https://github.com/rajewsky-lab/mirdeep2.git
mv mirdeep2 0.1.3; cd 0.1.3
source ~/.bashrc
perl install.pl 
```

# II. miRNAseq: raw data processing & aligment
## 1. Set up directory
```bash
```bash
p_raw='/path/to/raw'
p_ref='/path/to/reference/'
p_trim='/path/to/trimming'
p_map='/path/to/mapping'
p_genome='/path/to/ncbi_datasets/ncbi_dataset/data/GCA_012934285.1/'
p_adap='path/to/Trimmomatic-0.39/adapters'
bwa_sam_converter='/path/to/miRDeep2/0.1.3/src/bwa_sam_converter.pl'
hairpin='path/to/hairpin_cleanup.fa'

```
  
## 2. Running fastqc for each sample and multifastqc
```bash
fastqc -o ./fastqc *.fastq.gz # create a folder 'fastqc' if it does not exist to contan output of fastqc
cd fastqc
multiqc .
```
## 3. Raw data processing: Trimming & Filtering 
### Trimming
```bash
for file in "$p_raw"/*.fastq.gz; do
  output="${file%.fastq.gz}_trimmed.fastq.gz"
  output_file="$p_trim/$(basename "$output")"
  trimmomatic SE "$file" "$output_file" -phred33 -threads 4 ILLUMINACLIP:$adap/TruSeq3-SE.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:18
done

# Check fastqc again
mkdir $p_trim/FastQC_trimmed
cd $p_trim
fastqc -o ./FastQC_trimmed *.fastq.gz

# Multiqc again
cd $p_trim/FastQC_trimmed
multiqc . 

```


### Filtering
```bash
## Index with bowtie2 
sed -i '/^>/ s/ .*//' $p_genome/GCA_012934285.1_ASM1293428v1_genomic.fna
mkdir $p_ref/index
cd $p_ref/index
bowtie2-build $p_genome/GCA_012934285.1_ASM1293428v1_genomic.fna bowtie

## Alignment
for file in $p_trim/*_trimmed.fastq.gz; do
  output="${file%_trimmed.fastq.gz}.sam"
  output_file="$p_map/$(basename "$output")"
  bowtie2 -x "$p_ref/index/bowtie" -U "$file" -S "$output_file"
done

``` 

### Process alined sam to convert it to reads_collapsed format
```bash
# fasta - miRDeep2 using "reads_collapsed" format - we need to process alined sam to convert it to reads_collapsed format. As description from original publication, they extract mapped read only 
cd $p_map
for file in  $(ls $p_map); do
  $bwa_sam_converter -i "$file" -c -o "${file%.sam}.collapsed.fa" -a "${file%.sam}.collapsed_vs_genome.arf"
done

mkdir collapse 
mv *.fa collapse
mkdir arf
mv *.arf arf
```

### Running miRDeep2 with a sample
```bash
# Since this step requires a significant amount of time, we can run it with a single sample to observe the results.
miRDeep2='path/to/miRDeep2/0.1.3/src/miRDeep2.pl'
collapsefa='path/to/collapse/SRR12952898.collapsed.fa'
genome='path/to/ncbi_datasets/ncbi_dataset/data/GCA_012934285.1/GCA_012934285.1_ASM1293428v1_genomic.fna'
reads_collapsed_vs_genome='path/to/arf/SRR12952898.collapsed_vs_genome.arf'
hairpin='path/to/hairpin_cleanup.fa'
out='path/to/output'; mkdir -p $out; cd $out
time $miRDeep2 $collapsefa $genome $reads_collapsed_vs_genome $miRNAfa none $hairpin -d -c -v 2 > report.log 

```

### Running miRDeep2 with multi samples
```bash
mkdir -p $out
cd $out
sample_prefixes=("SRR12952898" "SRR12952899" "SRR12952900" "SRR12952901" "SRR12952902" "SRR12952903" "SRR12952904" "SRR12952905" "SRR12952906")
for prefix in "${sample_prefixes[@]}"; do
  collapsefa="/path/to/collapse/$prefix.reads_collapsed.fa"
  reads_collapsed_vs_genome="/path/to/arf/$prefix.reads_collapsed_vs_genome.arf"

  time "$miRDeep2" "$collapsefa" "$genome" "$reads_collapsed_vs_genome" none  "$hairpin" -d -c -v 2>report.log
done
```
