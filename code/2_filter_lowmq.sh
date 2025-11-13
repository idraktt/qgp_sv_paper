#!/usr/bin/bash
module load SURVIVOR/v1.0.5
module load bcftools/1.9
module load vcftools/0.1.15
module load svtyper/v0.7.0
module load java/1.8.121
module load vcflib/master
module load speedseq/0.1.2
module load CNVnator/v0.3.3
module load delly/0.7.8
module load svaba/v0.2.1
module load wham/master
module load SAMtools/1.3
module load breakdancer/1.4.5
module load python/2.7
module load ANNOTSV/2.0

filter_calls_lowmq() {
        while IFS='\n' read -r line || [[ -n "$line" ]]; do
        fname=$(basename "$line");
        name=$(echo $fname | cut -f 1 -d '.');
        IFS='- ' read -r -a array <<< "$name";
        family="${array[0]}";
        # echo -e "SURVIVOR filter $breakdancer_path/raw/$name.vcf $software_path/0.bed 50 10000000 -1 -1 $breakdancer_path/lowmq/$name.lowmq.vcf" | bsub -n 1 -e $reports_path/logs/$name-err.txt -o $reports_path/logs/$name-out.txt  -q qgp -P ${projectname};
        # echo -e "SURVIVOR filter $breakseq_path/raw/$name.vcf $software_path/0.bed 50 10000000 -1 -1 $breakseq_path/lowmq/$name.lowmq.vcf" | bsub -n 1 -e $reports_path/logs/$name-err.txt -o $reports_path/logs/$name-out.txt  -q qgp -P ${projectname};
        # echo -e "SURVIVOR filter $cnvnator_path/raw/$name.vcf $software_path/0.bed 50 10000000 -1 -1 $cnvnator_path/lowmq/$name.lowmq.vcf" | bsub -n 1 -e $reports_path/logs/$name-err.txt -o $reports_path/logs/$name-out.txt  -q qgp -P ${projectname};
        # echo -e "SURVIVOR filter $delly_path/raw/$name.vcf $software_path/0.bed 50 10000000 -1 -1 $delly_path/lowmq/$name.lowmq.vcf" | bsub -n 1 -e $reports_path/logs/$name-err.txt -o $reports_path/logs/$name-out.txt  -q qgp -P ${projectname};
        # echo -e "SURVIVOR filter $erds_path/raw/$name.vcf $software_path/0.bed 50 10000000 -1 -1 $erds_path/lowmq/$name.lowmq.vcf" | bsub -n 1 -e $reports_path/logs/$name-err.txt -o $reports_path/logs/$name-out.txt  -q qgp -P ${projectname};
        # echo -e "SURVIVOR filter $genomestrip_path/raw/$name.vcf $software_path/0.bed 50 10000000 -1 -1 $genomestrip_path/lowmq/$name.lowmq.vcf" | bsub -n 1 -e $reports_path/logs/$name-err.txt -o $reports_path/logs/$name-out.txt  -q qgp -P ${projectname};
        echo -e "SURVIVOR filter $manta_path/raw/$name.vcf $software_path/0.bed 50 10000000 -1 -1 $manta_path/lowmq/$name.lowmq.vcf" | bsub -n 1 -e $reports_path/logs/$name-err.txt -J ${uuid}filter_lowmq_manta_$name -o $reports_path/logs/$name-out.txt  -q qgp -P ${projectname};
        # echo -e "SURVIVOR filter $speedseq_path/raw/$name.vcf $software_path/0.bed 50 10000000 -1 -1 $speedseq_path/lowmq/$name.lowmq.vcf" | bsub -n 1 -e $reports_path/logs/$name-err.txt -o $reports_path/logs/$name-out.txt  -q qgp -P ${projectname};
         echo -e "SURVIVOR filter $svaba_path/raw/$name.vcf $software_path/0.bed 50 10000000 -1 -1 $svaba_path/lowmq/$name.lowmq.vcf" | bsub -n 1 -e $reports_path/logs/$name-err.txt -o $reports_path/logs/$name-out.txt  -q qgp -P ${projectname};
        # echo -e "SURVIVOR filter $wham_path/raw/$name.vcf $software_path/0.bed 50 10000000 -1 -1 $wham_path/lowmq/$name.lowmq.vcf" | bsub -n 1 -e $reports_path/logs/$name-err.txt -o $reports_path/logs/$name-out.txt  -q qgp -P ${projectname};

        done < "$1"
}

echo "Starting SIDRA SV Discovery";

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

echo $DIR;

uuid=$(cat /dev/urandom | tr -cd 'a-f0-9' | head -c 16);

echo $uuid;

bam_list=$1
work_dir=$2
projectname=$3

ref_data="/gpfs/data_jrnas1/ref_data/Homo_sapiens/hs37d5/Sequences/WholeGenomeSequence/hs37d5.fa"


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

filter_calls_lowmq $bam_list