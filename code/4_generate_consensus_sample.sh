#!/usr/bin/bash
module load SURVIVOR/v1.0.7
module load bcftools/1.9
module load bedtools/2.28

merge_consensuns_sample() {
	while IFS='' read -r line || [[ -n "$line" ]]; do
	echo "Text read from file: $line"
        fname=$(basename "$line");
        echo "$fname";
        name=$(echo $fname | cut -f 1 -d '.');
        echo "$name";
        IFS='- ' read -r -a array <<< "$name";
        family="${array[0]}";

        echo "$breakdancer_path/filtered/$name.filtered.vcf" > $reports_path/sample_lists/$name.filtered.list
		echo "$breakseq_path/filtered/$name.filtered.vcf" >> $reports_path/sample_lists/$name.filtered.list
		echo "$cnvnator_path/filtered/$name.filtered.vcf" >> $reports_path/sample_lists/$name.filtered.list
		echo "$delly_path/filtered/$name.filtered.vcf" >> $reports_path/sample_lists/$name.filtered.list
		echo "$erds_path/filtered/$name.filtered.vcf" >> $reports_path/sample_lists/$name.filtered.list
		echo "$genomestrip_path/filtered/$name.filtered.vcf" >> $reports_path/sample_lists/$name.filtered.list
		echo "$manta_path/filtered/$name.filtered.vcf" >> $reports_path/sample_lists/$name.filtered.list
		echo "$speedseq_path/filtered/$name.filtered.vcf" >> $reports_path/sample_lists/$name.filtered.list
		echo "$svaba_path/filtered/$name.filtered.vcf" >> $reports_path/sample_lists/$name.filtered.list
		echo "$wham_path/filtered/$name.filtered.vcf" >> $reports_path/sample_lists/$name.filtered.list

		echo "SURVIVOR merge $reports_path/sample_lists/$name.filtered.list 500 2 1 0 0 50 $reports_path/sample_final/$name.vcf;bedtools intersect -a $reports_path/sample_final/$name.vcf -b $software_path/blacklist_regions/segdup.bed $software_path/blacklist_regions/hg19-blacklist.bed -nonamecheck -wb -f 0.90 -r -v -header | bcftools annotate --set-id '%CHROM\_%POS\_%INFO/END\_%INFO/SVTYPE' > $reports_path/sample_final/$name.id.vcf ;rm -rf $reports_path/sample_final/$name.vcf" |bsub -n 1 -e $reports_path/logs/$name-err.txt -o $reports_path/logs/$name-out.txt -q qgp -P ${projectname};

 	done < "$1"
}

run_manta_insertions() {
	while IFS='' read -r row || [[ -n "$row" ]]; do
	echo "Text read from file: $line"
		line=$(echo $row | cut -f 1 -d ",");
        fname=$(basename "$line");
		name=$(echo $fname | cut -f 1 -d '.');
		echo "$name";
		echo "zgrep ^# $manta_path/$name/results/variants/diploidSV.vcf.gz > $reports_path/sample_insertions/$name.vcf;zgrep -v ^# $manta_path/$name/results/variants/diploidSV.vcf.gz | grep -v -E \"hs37d5|GL0|MT\" >> $reports_path/sample_insertions/$name.vcf;bedtools intersect -a $reports_path/sample_insertions/$name.vcf -b $software_path/blacklist_regions/segdup.bed $software_path/blacklist_regions/hg19-blacklist.bed $lowmq_path/$name.bed $software_path/blacklist_regions/hg19-simple_repeats.bed -nonamecheck -wb -f 0.50 -r -v -header | bcftools view -i '(SVTYPE = \"INS\" & (FORMAT/PR[0:1] > 8 | FORMAT/SR[0:1] > 8))' | bcftools annotate --set-id '%CHROM\_%POS\_%INFO/END\_%INFO/SVTYPE'  > $reports_path/sample_insertions/$name.ins.vcf;" | bsub -n 1 -e $reports_path/logs/$name-err.txt -o $reports_path/logs/$name-out.txt -q qgp -P ${projectname}

		done < $1
}

run_manta_translocations() {
	while IFS='' read -r row || [[ -n "$row" ]]; do
	echo "Text read from file: $line"
		line=$(echo $row | cut -f 1 -d ",");
        fname=$(basename "$line");
		name=$(echo $fname | cut -f 1 -d '.');
		echo "$name";
		echo "zgrep ^# $manta_path/$name/results/variants/diploidSV.vcf.gz > $reports_path/sample_translocations/$name.vcf;zgrep -v ^# $manta_path/$name/results/variants/diploidSV.vcf.gz | grep -v -E \"hs37d5|GL0|MT\" >> $reports_path/sample_translocations/$name.vcf;python /gpfs/software/genomics/manta/1.6.0/share/manta-1.6.0-0/libexec/convertInversion.py /gpfs/software/genomics/SAMtools/1.9/bin/samtools $ref_data $reports_path/sample_translocations/$name.vcf > $reports_path/sample_translocations/$name.tmp.vcf;bedtools intersect -a $reports_path/sample_translocations/$name.tmp.vcf -b $software_path/blacklist_regions/segdup.bed $software_path/blacklist_regions/hg19-blacklist.bed $lowmq_path/$name.bed $software_path/blacklist_regions/hg19-simple_repeats.bed -nonamecheck -wb -f 0.50 -r -v -header | bcftools view -i '(SVTYPE = \"BND\" & (FORMAT/PR[0:1] > 16 | FORMAT/SR[0:1] > 16) )' | bcftools annotate --set-id '%CHROM\_%POS\_%INFO/END\_%INFO/SVTYPE' > $reports_path/sample_translocations/$name.tra.vcf;" | bsub -n 1 -e $reports_path/logs/$name-err.txt -o $reports_path/logs/$name-out.txt -q qgp -P ${projectname}

		done < $1
}

run_manta_inversions() {
	while IFS='' read -r row || [[ -n "$row" ]]; do
	echo "Text read from file: $line"
		line=$(echo $row | cut -f 1 -d ",");
        fname=$(basename "$line");
		name=$(echo $fname | cut -f 1 -d '.');
		echo "$name";
		echo "zgrep ^# $manta_path/$name/results/variants/diploidSV.vcf.gz > $reports_path/sample_inversions/$name.vcf;zgrep -v ^# $manta_path/$name/results/variants/diploidSV.vcf.gz | grep -v -E \"hs37d5|GL0|MT\" >> $reports_path/sample_inversions/$name.vcf;python /gpfs/software/genomics/manta/1.6.0/share/manta-1.6.0-0/libexec/convertInversion.py /gpfs/software/genomics/SAMtools/1.9/bin/samtools $ref_data $reports_path/sample_inversions/$name.vcf > $reports_path/sample_inversions/$name.tmp.vcf;bedtools intersect -a $reports_path/sample_inversions/$name.tmp.vcf -b $software_path/blacklist_regions/segdup.bed $software_path/blacklist_regions/hg19-blacklist.bed $lowmq_path/$name.bed $software_path/blacklist_regions/hg19-simple_repeats.bed -nonamecheck -wb -f 0.50 -r -v -header | bcftools view -i '(SVTYPE = \"INV\" & SVLEN > 50 & SVLEN < 10000000 & (FORMAT/PR[0:1] > 12 | FORMAT/SR[0:1] > 12) )' | bcftools annotate --set-id '%CHROM\_%POS\_%INFO/END\_%INFO/SVTYPE' > $reports_path/sample_inversions/$name.inv.vcf;" | bsub -n 1 -e $reports_path/logs/$name-err.txt -o $reports_path/logs/$name-out.txt -q qgp -P ${projectname}

		done < $1
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

merge_consensuns_sample $bam_list

#run_manta_insertions $bam_list
#run_manta_translocations $bam_list
#run_manta_inversions $bam_list