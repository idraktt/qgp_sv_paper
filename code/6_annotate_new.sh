#!/bin/bash
module load svtyper/v0.7.0
module load ANNOTSV/2.2
module load bcftools/1.9
module load java/1.8.121
module load SURVIVOR/v1.0.7
module load vcftools/0.1.16_source
module load tabix/0.2.6
module load SVTK/29Apr2019


generate_del_dup_inv_annotated_vcf() {
	# ls $reports_path/sample_final/*hq.vcf > $reports_path/sample_final/sample.list;
	# SURVIVOR merge $reports_path/sample_final/sample.list 500 1 1 0 0 50 $reports_path/population_final/$projectname.survivor.vcf;
	# sed -i 's/\.\/\.:/0\/0:/g' $reports_path/population_final/$projectname.survivor.vcf;
	# bcftools annotate --set-id '%CHROM\_%POS\_%INFO/END\_%INFO/SVTYPE' $reports_path/population_final/$projectname.survivor.vcf | bcftools annotate -x "INFO/SUPP_VEC"  > $reports_path/population_final/$projectname.id.vcf;
	module load python/2.7.10
	python2.7 /gpfs/software/genomics/svtyper/v0.7.0/env/bin/vcf_allele_freq.py $reports_path/population_final/$projectname.id.vcf | grep -v -f $reports_path/population_final/blacklist_ids.list > $reports_path/population_final/$projectname.svtyper.vcf;
	java -jar -Xmx128G $software_path/htsSidra-1.0-jar-with-dependencies.jar annotatewithSubpops $reports_path/population_final/$projectname.svtyper.vcf /home/ealiyev_qgp/QF-QBB-RES-ACC-0032/pop_structure.tsv;
	cat $software_path/headers/subpop_header $reports_path/population_final/$projectname.svtyper.subpop.vcf | bcftools view -i 'POS < INFO/END[0]' | bcftools view -i 'INFO/NSAMP[0] > 0'> $reports_path/population_final/$projectname.AF.vcf;
	cat $software_path/headers/bed_header.txt > $reports_path/population_final/$projectname.bed;
	bcftools query -f "%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\t%INFO/SVLEN\t%ID\n" $reports_path/population_final/$projectname.AF.vcf>> $reports_path/population_final/$projectname.bed;
	cut -d$'\t' -f1-9 $reports_path/population_final/$projectname.AF.vcf > $reports_path/population_final/$projectname.sites.vcf
}

generate_tra_annotated_vcf() {
	ls $reports_path/sample_translocations/*tra.vcf > $reports_path/sample_translocations/translocations.list;
	SURVIVOR merge $reports_path/sample_translocations/translocations.list 500 1 1 0 0 0 $reports_path/population_final/$projectname.tra.vcf;
	sed -i 's/\.\/\.:/0\/0:/g' $reports_path/population_final/$projectname.tra.vcf;
	bcftools annotate --set-id '%CHROM\_%POS\_%INFO/END\_%INFO/SVTYPE' $reports_path/population_final/$projectname.tra.vcf | bcftools annotate -x "INFO/SUPP_VEC" > $reports_path/population_final/$projectname.tra.id.vcf;
	java -jar -Xmx128G $software_path/htsSidra-1.0-jar-with-dependencies.jar annotatewithSubpops $reports_path/population_final/$projectname.tra.id.vcf /home/ealiyev_qgp/QF-QBB-RES-ACC-0032/pop_structure.tsv;
	cat $software_path/headers/subpop_header $reports_path/population_final/$projectname.tra.id.subpop.vcf  > $reports_path/population_final/$projectname.tra.filtered.vcf;
	python2.7 /gpfs/software/genomics/svtyper/v0.7.0/env/bin/vcf_allele_freq.py $reports_path/population_final/$projectname.tra.filtered.vcf | bcftools view -i 'INFO/NSAMP[0] > 0' > $reports_path/population_final/$projectname.tra.AF.vcf;
	cat $software_path/headers/bed_header.txt > $reports_path/population_final/$projectname.tra.bed;
	bcftools query -f "%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\t%INFO/SVLEN\t%ID\n" $reports_path/population_final/$projectname.tra.AF.vcf >> $reports_path/population_final/$projectname.tra.bed;
	cut -d$'\t' -f1-9 $reports_path/population_final/$projectname.tra.AF.vcf > $reports_path/population_final/$projectname.tra.sites.vcf
}

generate_ins_annotated_vcf() {
	ls $reports_path/sample_insertions/*ins.vcf > $reports_path/sample_insertions/insertions.list;
	SURVIVOR merge $reports_path/sample_insertions/insertions.list 500 1 1 0 0 0 $reports_path/population_final/$projectname.ins.vcf;
	sed -i 's/\.\/\.:/0\/0:/g' $reports_path/population_final/$projectname.ins.vcf;
	bcftools annotate --set-id '%CHROM\_%POS\_%INFO/END\_%INFO/SVTYPE' $reports_path/population_final/$projectname.ins.vcf | bcftools annotate -x "INFO/SUPP_VEC" > $reports_path/population_final/$projectname.ins.id.vcf;
	java -jar -Xmx128G $software_path/htsSidra-1.0-jar-with-dependencies.jar annotatewithSubpops $reports_path/population_final/$projectname.ins.id.vcf /home/ealiyev_qgp/QF-QBB-RES-ACC-0032/pop_structure.tsv;
	cat $software_path/headers/subpop_header $reports_path/population_final/$projectname.ins.id.subpop.vcf  > $reports_path/population_final/$projectname.ins.filtered.vcf;
	python2.7 /gpfs/software/genomics/svtyper/v0.7.0/env/bin/vcf_allele_freq.py $reports_path/population_final/$projectname.ins.filtered.vcf | bcftools view -i 'POS < INFO/END[0]' | bcftools view -i 'INFO/NSAMP[0] > 0' > $reports_path/population_final/$projectname.ins.AF.vcf;

	cat $software_path/headers/bed_header.txt > $reports_path/population_final/$projectname.ins.bed;
	bcftools query -f "%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\t%INFO/SVLEN\t%ID\n" $reports_path/population_final/$projectname.ins.AF.vcf >> $reports_path/population_final/$projectname.ins.bed;
	cut -d$'\t' -f1-9 $reports_path/population_final/$projectname.ins.AF.vcf > $reports_path/population_final/$projectname.ins.sites.vcf
}


start_annotation() { 
	grep -v "#" $reports_path/population_final/$projectname.tra.sites.vcf >> $reports_path/population_final/$projectname.sites.vcf;
	grep -v "#" $reports_path/population_final/$projectname.ins.sites.vcf >> $reports_path/population_final/$projectname.sites.vcf;
	export LC_ALL=C ;
	AnnotSV -SVinputFile $reports_path/population_final/$projectname.sites.vcf -bedtools /gpfs/software/genomics/bedtools/2.26/bin/bedtools -SVinputInfo 1 -outputFile $projectname -outputDir $reports_path/population_final/;
	java -jar -Xmx64G $software_path/htsSidra-1.0-jar-with-dependencies.jar processAnnotatedFile $reports_path/population_final/$projectname.tsv
}


generate__annotation_counts() {
	cat $software_path/headers/zygosity_header.txt > $reports_path/population_final/$projectname.zygosity_counts.tsv;cat $reports_path/population_final/$projectname.AF.vcf | grep -v ^# | awk '{ unknown=0; homref=0; het=0; homalt=0; for(i=10;i<=NF;i++) { split($i,a,":"); split(a[1],GT,"[/|]");if(GT[1]=="."&&GT[2]==".") {unknown++} else if(GT[1]==0&&GT[2]==0) {homref++} else if(GT[1]==GT[2]) {homalt++} else {het++}};print $3,"\t",unknown,"\t",homref,"\t",het,"\t",homalt}' >> $reports_path/population_final/$projectname.zygosity_counts.tsv;
	cat $software_path/headers/zygosity_header.txt > $reports_path/population_final/zygosity_counts.tsv;cat $reports_path/population_final/$projectname.tra.AF.vcf | grep -v ^# | awk '{ unknown=0; homref=0; het=0; homalt=0; for(i=10;i<=NF;i++) { split($i,a,":"); split(a[1],GT,"[/|]");if(GT[1]=="."&&GT[2]==".") {unknown++} else if(GT[1]==0&&GT[2]==0) {homref++} else if(GT[1]==GT[2]) {homalt++} else {het++}};print $1,"\t",$2,"\t",$3,"\t",$4,"\t",$5,unknown,homref,het,homalt}' >> $reports_path/population_final/$projectname.tra.zygosity_counts.tsv;
	cat $software_path/headers/zygosity_header.txt > $reports_path/population_final/zygosity_counts.tsv;cat $reports_path/population_final/$projectname.ins.AF.vcf | grep -v ^# | awk '{ unknown=0; homref=0; het=0; homalt=0; for(i=10;i<=NF;i++) { split($i,a,":"); split(a[1],GT,"[/|]");if(GT[1]=="."&&GT[2]==".") {unknown++} else if(GT[1]==0&&GT[2]==0) {homref++} else if(GT[1]==GT[2]) {homalt++} else {het++}};print $1,"\t",$2,"\t",$3,"\t",$4,"\t",$5,unknown,homref,het,homalt}' >> $reports_path/population_final/$projectname.ins.zygosity_counts.tsv;

	bcftools query -f "%ID\t%INFO/AF\t%INFO/NSAMP\n" $reports_path/population_final/$projectname.AF.vcf > $reports_path/population_final/$projectname.allele_freq.tsv;
	bcftools query -f "%ID\t%INFO/AF\t%INFO/NSAMP\n" $reports_path/population_final/$projectname.ins.AF.vcf > $reports_path/population_final/$projectname.ins.allele_freq.tsv;
	bcftools query -f "%ID\t%INFO/AF\t%INFO/NSAMP\n" $reports_path/population_final/$projectname.tra.AF.vcf > $reports_path/population_final/$projectname.tra.allele_freq.tsv;

	bcftools query -f "%ID\t%INFO/AF\t%INFO/NSAMP\t%INFO/QGP_ADM_AC\t%INFO/QGP_ADM_AF\t%INFO/QGP_ADM_AN\t%INFO/QGP_ADM_HOMO\t%INFO/QGP_AFR_AC\t%INFO/QGP_AFR_AF\t%INFO/QGP_AFR_AN\t%INFO/QGP_AFR_HOMO\t%INFO/QGP_GAR_AC\t%INFO/QGP_GAR_AF\t%INFO/QGP_GAR_AN\t%INFO/QGP_GAR_HOMO\t%INFO/QGP_PAR_AC\t%INFO/QGP_PAR_AF\t%INFO/QGP_PAR_AN\t%INFO/QGP_PAR_HOMO\t%INFO/QGP_SAS_AC\t%INFO/QGP_SAS_AF\t%INFO/QGP_SAS_AN\t%INFO/QGP_SAS_HOMO\t%INFO/QGP_WEP_AC\t%INFO/QGP_WEP_AF\t%INFO/QGP_WEP_AN\t%INFO/QGP_WEP_HOMO\n" $reports_path/population_final/$projectname.AF.vcf > $reports_path/population_final/${projectname}_allele_freq_subpop.tsv;

	echo -E 'ID	AF	NSAMP	QGP_ADM_AC	QGP_ADM_AF	QGP_ADM_AN	QGP_ADM_HOMO	QGP_AFR_AC	QGP_AFR_AF	QGP_AFR_AN	QGP_AFR_HOMO	QGP_GAR_AC	QGP_GAR_AF	QGP_GAR_AN	QGP_GAR_HOMO	QGP_PAR_AC	QGP_PAR_AF	QGP_PAR_AN	QGP_PAR_HOMO	QGP_SAS_AC	QGP_SAS_AF	QGP_SAS_AN	QGP_SAS_HOMO	QGP_WEP_AC	QGP_WEP_AF	QGP_WEP_AN	QGP_WEP_HOMO\n' | cat - $reports_path/population_final/${projectname}_allele_freq_subpop.tsv > temp && mv temp $reports_path/population_final/${projectname}_allele_freq_subpop.tsv;
	
	bcftools query -f "%ID\t%INFO/AF\t%INFO/NSAMP\t%INFO/QGP_ADM_AC\t%INFO/QGP_ADM_AF\t%INFO/QGP_ADM_AN\t%INFO/QGP_ADM_HOMO\t%INFO/QGP_AFR_AC\t%INFO/QGP_AFR_AF\t%INFO/QGP_AFR_AN\t%INFO/QGP_AFR_HOMO\t%INFO/QGP_GAR_AC\t%INFO/QGP_GAR_AF\t%INFO/QGP_GAR_AN\t%INFO/QGP_GAR_HOMO\t%INFO/QGP_PAR_AC\t%INFO/QGP_PAR_AF\t%INFO/QGP_PAR_AN\t%INFO/QGP_PAR_HOMO\t%INFO/QGP_SAS_AC\t%INFO/QGP_SAS_AF\t%INFO/QGP_SAS_AN\t%INFO/QGP_SAS_HOMO\t%INFO/QGP_WEP_AC\t%INFO/QGP_WEP_AF\t%INFO/QGP_WEP_AN\t%INFO/QGP_WEP_HOMO\n" $reports_path/population_final/$projectname.ins.AF.vcf > $reports_path/population_final/${projectname}_allele_freq_subpop.ins.tsv;

	echo -E 'ID	AF	NSAMP	QGP_ADM_AC	QGP_ADM_AF	QGP_ADM_AN	QGP_ADM_HOMO	QGP_AFR_AC	QGP_AFR_AF	QGP_AFR_AN	QGP_AFR_HOMO	QGP_GAR_AC	QGP_GAR_AF	QGP_GAR_AN	QGP_GAR_HOMO	QGP_PAR_AC	QGP_PAR_AF	QGP_PAR_AN	QGP_PAR_HOMO	QGP_SAS_AC	QGP_SAS_AF	QGP_SAS_AN	QGP_SAS_HOMO	QGP_WEP_AC	QGP_WEP_AF	QGP_WEP_AN	QGP_WEP_HOMO\n' | cat - $reports_path/population_final/${projectname}_allele_freq_subpop.ins.tsv > temp && mv temp $reports_path/population_final/${projectname}_allele_freq_subpop.ins.tsv;
	
	bcftools query -f "%ID\t%INFO/AF\t%INFO/NSAMP\t%INFO/QGP_ADM_AC\t%INFO/QGP_ADM_AF\t%INFO/QGP_ADM_AN\t%INFO/QGP_ADM_HOMO\t%INFO/QGP_AFR_AC\t%INFO/QGP_AFR_AF\t%INFO/QGP_AFR_AN\t%INFO/QGP_AFR_HOMO\t%INFO/QGP_GAR_AC\t%INFO/QGP_GAR_AF\t%INFO/QGP_GAR_AN\t%INFO/QGP_GAR_HOMO\t%INFO/QGP_PAR_AC\t%INFO/QGP_PAR_AF\t%INFO/QGP_PAR_AN\t%INFO/QGP_PAR_HOMO\t%INFO/QGP_SAS_AC\t%INFO/QGP_SAS_AF\t%INFO/QGP_SAS_AN\t%INFO/QGP_SAS_HOMO\t%INFO/QGP_WEP_AC\t%INFO/QGP_WEP_AF\t%INFO/QGP_WEP_AN\t%INFO/QGP_WEP_HOMO\n" $reports_path/population_final/$projectname.tra.AF.vcf > $reports_path/population_final/${projectname}_allele_freq_subpop.tra.tsv;

	echo -E 'ID	AF	NSAMP	QGP_ADM_AC	QGP_ADM_AF	QGP_ADM_AN	QGP_ADM_HOMO	QGP_AFR_AC	QGP_AFR_AF	QGP_AFR_AN	QGP_AFR_HOMO	QGP_GAR_AC	QGP_GAR_AF	QGP_GAR_AN	QGP_GAR_HOMO	QGP_PAR_AC	QGP_PAR_AF	QGP_PAR_AN	QGP_PAR_HOMO	QGP_SAS_AC	QGP_SAS_AF	QGP_SAS_AN	QGP_SAS_HOMO	QGP_WEP_AC	QGP_WEP_AF	QGP_WEP_AN	QGP_WEP_HOMO\n' | cat - $reports_path/population_final/${projectname}_allele_freq_subpop.tra.tsv > temp && mv temp $reports_path/population_final/${projectname}_allele_freq_subpop.tra.tsv;

}

extra_experimental_shit() {

	cat $software_path/headers/report_header > $reports_path/population_final/$projectname.persample_counts.all.tsv;java -jar -Xmx128G $software_path/htsSidra-1.0-jar-with-dependencies.jar countgtmultisampleVCF $reports_path/population_final/$projectname.AF.vcf $reports_path/population_final/$projectname.tsv.processed.tsv ALL >> $reports_path/population_final/$projectname.persample_counts.all.tsv;

	cat $software_path/headers/report_header > $reports_path/population_final/$projectname.persample_counts.genic.tsv;java -jar -Xmx128G $software_path/htsSidra-1.0-jar-with-dependencies.jar countgtmultisampleVCF $reports_path/population_final/$projectname.AF.vcf $reports_path/population_final/$projectname.tsv.processed.tsv GENIC >> $reports_path/population_final/$projectname.persample_counts.genic.tsv;

	cat $software_path/headers/report_header > $reports_path/population_final/$projectname.persample_counts.exonic.tsv;java -jar -Xmx128G $software_path/htsSidra-1.0-jar-with-dependencies.jar countgtmultisampleVCF $reports_path/population_final/$projectname.AF.vcf $reports_path/population_final/$projectname.tsv.processed.tsv EXONIC >> $reports_path/population_final/$projectname.persample_counts.exonic.tsv;

	cat $software_path/headers/report_header > $reports_path/population_final/$projectname.persample_counts.novel.tsv;java -jar -Xmx128G $software_path/htsSidra-1.0-jar-with-dependencies.jar countgtmultisampleVCF $reports_path/population_final/$projectname.AF.vcf $reports_path/population_final/$projectname.tsv.processed.tsv NOVEL >> $reports_path/population_final/$projectname.persample_counts.novel.tsv;

	cat $software_path/headers/report_header > $reports_path/population_final/$projectname.persample_counts.novelgenic.tsv;java -jar -Xmx128G $software_path/htsSidra-1.0-jar-with-dependencies.jar countgtmultisampleVCF $reports_path/population_final/$projectname.AF.vcf $reports_path/population_final/$projectname.tsv.processed.tsv NOVEL_GENIC >> $reports_path/population_final/$projectname.persample_counts.novelgenic.tsv;

	cat $software_path/headers/report_header > $reports_path/population_final/$projectname.persample_counts.novelexonic.tsv;java -jar -Xmx128G $software_path/htsSidra-1.0-jar-with-dependencies.jar countgtmultisampleVCF $reports_path/population_final/$projectname.AF.vcf $reports_path/population_final/$projectname.tsv.processed.tsv NOVEL_EXONIC >> $reports_path/population_final/$projectname.persample_counts.novelexonic.tsv;

	module unload bedtools/2.28
	module load bedtools/2.26

	vcf-sort $reports_path/population_final/$projectname.AF.vcf -t $reports_path/population_final/ > $reports_path/population_final/$projectname.sorted.no_supp_vec.sorted.vcf;
	
	# bgzip $reports_path/population_final/$projectname.sorted.no_supp_vec.sorted.vcf;
	# tabix -p vcf $reports_path/population_final/$projectname.sorted.no_supp_vec.sorted.vcf.gz;

	sed -i 's/SVTYPE=TRA/SVTYPE=BND/g' $reports_path/population_final/$projectname.sorted.no_supp_vec.sorted.vcf;

	svtk count-svtypes $reports_path/population_final/$projectname.sorted.no_supp_vec.sorted.vcf $reports_path/population_final/$projectname.total_counts.txt

	svtk count-svtypes $reports_path/population_final/$projectname.ins.AF.vcf $reports_path/population_final/$projectname.total_counts.ins.txt

	svtk count-svtypes $reports_path/population_final/$projectname.tra.AF.vcf $reports_path/population_final/$projectname.total_counts.tra.txt

	svtk count-svtypes $reports_path/population_final/$projectname.sorted.no_supp_vec.sorted.vcf $reports_path/population_final/$projectname.total_obs.txt --total-obs

	svtk count-svtypes $reports_path/population_final/$projectname.sorted.no_supp_vec.sorted.vcf $reports_path/population_final/$projectname.total_variants.txt --total-variants

	bcftools view -i 'INFO/NSAMP[0] == 1' $reports_path/population_final/$projectname.sorted.no_supp_vec.sorted.vcf > $reports_path/population_final/${projectname}.singleton.vcf;svtk count-svtypes $reports_path/population_final/${projectname}.singleton.vcf $reports_path/population_final/${projectname}_total_counts.singleton.txt;rm -rf $reports_path/population_final/${projectname}.singleton.vcf;

	bcftools view -i 'INFO/AF[0] < 0.01' $reports_path/population_final/$projectname.sorted.no_supp_vec.sorted.vcf > $reports_path/population_final/${projectname}.rare.vcf;svtk count-svtypes $reports_path/population_final/${projectname}.rare.vcf $reports_path/population_final/${projectname}_total_counts.rare.txt;rm -rf $reports_path/population_final/${projectname}.rare.vcf;

	bcftools view -i 'INFO/QGP_ADM_HOMO[0]>1 | INFO/QGP_AFR_HOMO[0]>1 | INFO/QGP_GAR_HOMO[0]>1 | INFO/QGP_PAR_HOMO[0]>1 | INFO/QGP_SAS_HOMO[0]>1 | INFO/QGP_WEP_HOMO[0]>1' $reports_path/population_final/$projectname.sorted.no_supp_vec.sorted.vcf > $reports_path/population_final/${projectname}.homo.vcf;svtk count-svtypes $reports_path/population_final/${projectname}.homo.vcf $reports_path/population_final/${projectname}_total_counts.homo.txt;rm -rf $reports_path/population_final/${projectname}.homo.vcf;

	bcftools view -i 'INFO/SVLEN[0] < -2000000 | INFO/SVLEN[0] > 2000000' $reports_path/population_final/$projectname.sorted.no_supp_vec.sorted.vcf > $reports_path/population_final/${projectname}.bigevents.vcf;

	#svtk standardize $reports_path/population_final/$projectname.sorted.no_supp_vec.sorted.vcf $reports_path/population_final/$projectname.sorted.no_supp_vec.std.vcf manta;

	#svtk annotate $reports_path/population_final/$projectname.sites.vcf $reports_path/population_final/${projectname}.svtk.vcf --gencode $software_path/gencode.v19.annotation.gtf;svtk vcf2bed --info ALL $reports_path/population_final/${projectname}.svtk.vcf $reports_path/population_final/${projectname}.svtk.bed;

	# grep "#CHROM" $reports_path/population_final/$projectname.sorted.no_supp_vec.sorted.vcf | cut -f3,10- > $reports_path/population_final/$projectname.genotypes.vcf;bcftools query -f'%ID[\t%GT]\n' $$reports_path/population_final/$projectname.sorted.no_supp_vec.sorted.vcf >> $reports_path/population_final/$projectname.genotypes.vcf;

	# sed -i 's/0\/1/0|1/g' $reports_path/population_final/$projectname.genotypes.vcf;sed -i 's/1\/1/1|1/g' $reports_path/population_final/$projectname.genotypes.vcf;


	# join -t $'\t' --header -1 8 -2 1 <(head -n 1 $reports_path/population_final/$projectname.tsv.processed.tsv && tail -n 2 $reports_path/population_final/$projectname.tsv.processed.tsv | sort -k8 $reports_path/population_final/$projectname.tsv.processed.tsv) <(head -n 1 $reports_path/population_final/$projectname.genotypes.vcf && tail -n 2 $reports_path/population_final/$projectname.genotypes.vcf | sort -k1 $reports_path/population_final/$projectname.genotypes.vcf) > $reports_path/population_final/$projectname.merged.tsv;

	#java -jar -Xmx64G $software_path/htsSidra-1.0-jar-with-dependencies.jar createDenovoHomozygousReport $reports_path/population_final/$projectname.AF.vcf $software_path/mnd_ped_merged.ped $software_path/mnd_sampleid_merged.txt > $reports_path/population_final/$projectname.mutations.txt
}

echo "Starting SIDRA SV Annotation";

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

# generate_del_dup_inv_annotated_vcf
#generate_ins_annotated_vcf
#generate_tra_annotated_vcf
# start_annotation
# generate__annotation_counts
extra_experimental_shit



