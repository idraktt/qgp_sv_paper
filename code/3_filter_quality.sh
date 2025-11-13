#!/usr/bin/bash
module load bcftools/1.9

filter_calls_quality() {
        while IFS='\n' read -r line || [[ -n "$line" ]]; do
        fname=$(basename "$line");
        name=$(echo $fname | cut -f 1 -d '.');
        IFS='- ' read -r -a array <<< "$name";
        family="${array[0]}";
        echo -e "bcftools view -i '(SVTYPE = \"DEL\" & (INFO/num_Reads[0] > 10) )' $breakdancer_path/lowmq/$name.lowmq.vcf > $breakdancer_path/filtered/$name.filtered.vcf" | bsub -n 1 -e $reports_path/logs/$name-err.txt -o $reports_path/logs/$name-out.txt -q qgp -P ${projectname};

        echo -e "bcftools view -i '(SVTYPE = \"DEL\" & (FORMAT/PE[0] > 0) )' $breakseq_path/lowmq/$name.lowmq.vcf > $breakseq_path/filtered/$name.filtered.vcf" | bsub -n 1 -e $reports_path/logs/$name-err.txt -o $reports_path/logs/$name-out.txt -q qgp -P ${projectname};

        echo -e "bcftools view -i 'INFO/natorQ0[0] < 0.5 & INFO/natorP1[0] < 0.01' $cnvnator_path/lowmq/$name.lowmq.vcf > $cnvnator_path/filtered/$name.filtered.vcf" | bsub -n 1 -e $reports_path/logs/$name-err.txt -o $reports_path/logs/$name-out.txt -q qgp -P ${projectname};

        echo -e "bcftools view -i '(SVTYPE = \"DEL\" & (INFO/PE[0] > 12 | INFO/SR[0] > 12) ) | (SVTYPE = \"DUP\" & (INFO/PE[0] > 4 | INFO/SR[0] > 4)) | (SVTYPE = \"INV\"& (INFO/PE[0] > 4 | INFO/SR[0] > 4)) ' $delly_path/lowmq/$name.lowmq.vcf > $delly_path/filtered/$name.filtered.vcf" | bsub -n 1 -e $reports_path/logs/$name-err.txt -o $reports_path/logs/$name-out.txt -q qgp -P ${projectname};

        echo -e "bcftools view $erds_path/lowmq/$name.lowmq.vcf > $erds_path/filtered/$name.filtered.vcf" | bsub -n 1 -e $reports_path/logs/$name-err.txt -o $reports_path/logs/$name-out.txt -q qgp -P ${projectname};

        echo -e "bcftools view  $genomestrip_path/lowmq/$name.lowmq.vcf > $genomestrip_path/filtered/$name.filtered.vcf" | bsub -n 1 -e $reports_path/logs/$name-err.txt -o $reports_path/logs/$name-out.txt -q qgp -P ${projectname};

        echo -e "bcftools view -i '(SVTYPE = \"DEL\" & (FORMAT/PR[0:1] > 16 | FORMAT/SR[0:1] > 16) ) | (SVTYPE = \"DUP\" & (FORMAT/PR[0:1] > 16 | FORMAT/SR[0:1] > 16)) | (SVTYPE = \"INV\" & (FORMAT/PR[0:1] > 12 | FORMAT/SR[0:1] > 12))' $manta_path/lowmq/$name.lowmq.vcf > $manta_path/filtered/$name.filtered.vcf" | bsub -n 1 -e $reports_path/logs/$name-err.txt -o $reports_path/logs/$name-out.txt -q qgp -P ${projectname};
        
        echo -e "bcftools view -i '(SVTYPE = \"DEL\" & (INFO/PE[0] > 4 | INFO/SR[0] > 4)) | (SVTYPE = \"DUP\" & (INFO/PE[0] > 5 | INFO/SR[0] > 5)) | (SVTYPE = \"INV\"& (INFO/PE[0] > 10 | INFO/SR[0] > 10))' $speedseq_path/lowmq/$name.lowmq.vcf > $speedseq_path/filtered/$name.filtered.vcf" | bsub -n 1 -e $reports_path/logs/$name-err.txt -o $reports_path/logs/$name-out.txt -q qgp -P ${projectname};
        
        echo -e "bcftools view -i '(SVTYPE = \"DEL\" & (FORMAT/AD[0] > 12) ) | (SVTYPE = \"DUP\" & (FORMAT/AD[0] > 12)) | (SVTYPE = \"INV\"& (FORMAT/AD[0] > 12)) ' $svaba_path/lowmq/$name.lowmq.vcf > $svaba_path/filtered/$name.filtered.vcf" | bsub -n 1 -e $reports_path/logs/$name-err.txt -o $reports_path/logs/$name-out.txt -q qgp -P ${projectname};
        
        echo -e "bcftools view -i '(SVTYPE = \"DEL\" & (INFO/A[0] > 8) ) | (SVTYPE = \"DUP\" & (INFO/A[0] > 8))| (SVTYPE = \"INV\" & (INFO/A[0] > 6 ))' $wham_path/lowmq/$name.lowmq.vcf > $wham_path/filtered/$name.filtered.vcf" | bsub -n 1 -e $reports_path/logs/$name-err.txt -o $reports_path/logs/$name-out.txt -q qgp -P ${projectname};

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

filter_calls_quality $bam_list
