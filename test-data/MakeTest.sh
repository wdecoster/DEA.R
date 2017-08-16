# !/bin/sh
# Script for generating test data for DEA.R script
# wdecoster

if [ -z "$1" ]
  then
      echo "Please specify number of processes to be used by MakeTest.sh"
      echo "Example for 8 cores: 'bash MakeTest.sh 8'"
      exit 1
fi
echo "Using $1 cores for MakeTest.sh"

mkdir MakeTestDir || { echo >&2 "Could not make working directory, do you have write permissions to the currect directory?" ; exit 1 ; }
cd MakeTestDir

echo "Testing if all required software is available and in the path..."
hash STAR 2>/dev/null || { echo >&2 "Dependency Error: MakeTest requires STAR to be installed and in the path."; exit 1; }
hash wget 2>/dev/null || { echo >&2 "Dependency Error: MakeTest requires wget to be installed and in the path."; exit 1; }
echo -e "It seems all required software is present.\n"

echo "Downloading reference genome and annotation from Ensembl..."
wget -nv -O annotation.gtf.gz ftp://ftp.ensembl.org/pub/grch37/release-88/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz \
    && gunzip -f annotation.gtf.gz || { echo "Downloading of annotation failed." ; exit 1; }  &
wget -nv -O genome.fa.gz ftp://ftp.ensembl.org/pub/grch37/release-88/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz \
    && gunzip -f genome.fa.gz || { echo "Downloading of reference genome failed." ; exit 1; }
GTF=$(readlink -f annotation.gtf)
REF=$(readlink -f genome.fa)
echo -e "Reference genome and annotation downloaded.\n"


echo "Downloading fastq files for test data from ENA. This will take a while..."
wget -nv ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/008/SRR1039508/SRR1039508_?.fastq.gz
wget -nv ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/009/SRR1039509/SRR1039509_?.fastq.gz
wget -nv ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/002/SRR1039512/SRR1039512_?.fastq.gz
wget -nv ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/003/SRR1039513/SRR1039513_?.fastq.gz
wget -nv ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/006/SRR1039516/SRR1039516_?.fastq.gz
wget -nv ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/007/SRR1039517/SRR1039517_?.fastq.gz
wget -nv ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/000/SRR1039520/SRR1039520_?.fastq.gz
wget -nv ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/001/SRR1039521/SRR1039521_?.fastq.gz
echo -e "Fastq files downloaded.\n"

echo "Building STAR genome index..."
GDIR=genome_STAR
mkdir $GDIR
STAR --runThreadN $1 \
--runMode genomeGenerate \
--genomeDir $GDIR \
--genomeFastaFiles $REF \
--sjdbGTFfile $GTF
echo -e "STAR genome index complete.\n"

echo "Starting alignment with STAR..."
for sample in $(ls *_*.fastq.gz | cut -d'_' -f1 | sort -u)
do
STAR --runThreadN $1 \
--genomeDir $GDIR \
--readFilesCommand zcat \
--readFilesIn ${sample}_1.fastq.gz ${sample}_2.fastq.gz \
--outFileNamePrefix ${sample}_ \
--outSAMtype BAM SortedByCoordinate ;
done ;
echo -e "Alignment completed.\n"

echo "Creating sample info file..."
echo -e "sample\tcondition\tcell\tsequencing\tstrandedness" > test-sample-info.interm
echo -e "SRR1039508\tCON\tA\tPE\tunstranded" >> test-sample-info.interm
echo -e "SRR1039509\tTRT\tA\tPE\tunstranded" >> test-sample-info.interm
echo -e "SRR1039512\tCON\tB\tPE\tunstranded" >> test-sample-info.interm
echo -e "SRR1039513\tTRT\tB\tPE\tunstranded" >> test-sample-info.interm
echo -e "SRR1039516\tCON\tC\tPE\tunstranded" >> test-sample-info.interm
echo -e "SRR1039517\tTRT\tC\tPE\tunstranded" >> test-sample-info.interm
echo -e "SRR1039520\tCON\tD\tPE\tunstranded" >> test-sample-info.interm
echo -e "SRR1039521\tTRT\tD\tPE\tunstranded" >> test-sample-info.interm
echo "file" > bams
for f in *.bam ; do readlink -f $f >> bams ; done
paste bams test-sample-info.interm > test-sample-info.txt
rm test-sample-info.interm bams
SINFO=$(readlink -f test-sample-info.txt)
echo "Sample info file ready."

echo "Cleaning up intermediate files..."
rm *.fastq.gz *.out
cd ..
echo "MakeTest.sh is finished, ready to start testing DEA.R using:"
echo "DEA.R $SINFO $GTF"
