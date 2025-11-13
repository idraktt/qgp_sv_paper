#!/usr/bin/bash
module load bcftools/1.9
module load vcftools/0.1.15
module load svtyper/v0.7.0
module load java/1.8.121
module load vcflib/master
module load SAMtools/1.3
module load ANNOTSV/2.0
module load htslib/1.9

echo "Starting SIDRA SV Discovery";

genotype_samples() {
	while IFS='' read -r line || [[ -n "$line" ]]; do
        fname=$(basename "$line");
        name=$(echo $fname | cut -f 1 -d '.');
        IFS='- ' read -r -a array <<< "-t /";
        family="${array[0]}";
        #java -jar -Xmx64G $software_path/htsSidra-1.0-jar-with-dependencies.jar adjustALT $reports_path/sample_final/$name.hq.vcf;
	echo "cut -d$'\t' -f1-9 $reports_path/sample_final/$name.id.vcf > $reports_path/sample_final/$name.cnvs.vcf;svtyper -i $reports_path/sample_final/$name.cnvs.vcf -B $line > $reports_path/sample_final/$name.gt.vcf;DUPHOLD_SAMPLE_NAME=$name;$software_path/duphold -v $reports_path/sample_final/$name.gt.vcf -b $line -f $ref_data -t 4 -o  $reports_path/sample_final/$name.gt.duphold.vcf;bcftools view -i '(SVTYPE = \"DEL\" & (FMT/DHFFC[0] < 0.7 & FMT/GT[0] != \"0/0\")) | (SVTYPE = \"DUP\" | (FMT/DHBFC[0] > 1.2 & FMT/GT[0] != \"0/0\")) | (SVTYPE = \"INV\" & FMT/GT[0] != \"0/0\")' $reports_path/sample_final/$name.gt.duphold.vcf > $reports_path/sample_final/$name.hq.vcf;" | bsub -n 2 -e $reports_path/logs/$name-err.txt -o $reports_path/logs/$name-out.txt -q ngstest -P ${projectname};
	done < "$1"
}


DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

echo $DIR;

uuid=$(cat /dev/urandom | tr -cd 'a-f0-9' | head -c 16);

echo $uuid;

bam_list=$1
work_dir=$2
projectname=$3
vcf_snp_path=$4

ref_data="/home/ealiyev_qgp/QF-QBB-RES-ACC-0032/ref_data/Homo_sapiens/hs37d5/Sequences/WholeGenomeSequence/hs37d5.fa"


software_path=$DIR;
cnvnator_chromosomes_path="$software_path/cnvnator_chromosomes";

breakdancer_path="$work_dir/breakdancer"
breakseq_path="$work_dir/breakseq"
cnvnator_path="$work_dir/cnvnator"
delly_path="$work_dir/delly"
erds_path="$work_dir/erds"
genomestrip_path="$work_dir/genomestrip"
lowmq_path="$work_dir/lowmq"
manta_path="$work_dir/manta"
svaba_path="$work_dir/svaba"
speedseq_path="$work_dir/speedseq"
wham_path="$work_dir/wham"

reports_path="$work_dir/reports"

module load svtyper/v0.7.0

genotype_samples $bam_list