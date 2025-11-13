#!/usr/bin/bash
module load SURVIVOR/v1.0.3
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



generate_lowmq_bed() {
	while IFS='\n' read -r line || [[ -n "$line" ]]; do
        fname=$(basename "$line");
        name=$(echo $fname | cut -f 1 -d '.');
        IFS='- ' read -r -a array <<< "$name";
        family="${array[0]}";
        echo -e "samtools view -H $line > $lowmq_path/$name.sam;samtools view $line | awk '\$5<5 {print \$0}' >>  $lowmq_path/$name.sam;samtools view -S -b -h $lowmq_path/$name.sam > $lowmq_path/$name.bam;samtools depth $lowmq_path/$name.bam >  $lowmq_path/$name.cov;SURVIVOR bincov $lowmq_path/$name.cov 10 2 > $lowmq_path/$name.bed;" | bsub -n 2 -J ${uuid}generate_lowmq_${name} -e $reports_path/logs/$name-err.txt -o $reports_path/logs/$name-out.txt  -P ${projectname}
	done < "$1"
}

filter_calls_lowmq() {
	while IFS='\n' read -r line || [[ -n "$line" ]]; do
        fname=$(basename "$line");
        name=$(echo $fname | cut -f 1 -d '.');
        IFS='- ' read -r -a array <<< "$name";
        family="${array[0]}";
        echo -e "SURVIVOR filter $breakdancer_path/raw/$name.vcf $lowmq_path/$name.bed 300 2000000 -1 -1 $breakdancer_path/lowmq/$name.lowmq.vcf" | bsub -n 1 -J ${uuid}filter_lowmq_breakdancer_${name} -w 'done('$uuid'breakdancer_run_'$name') && done('$uuid'generate_lowmq_'$name')' -e $reports_path/logs/$name-err.txt -o $reports_path/logs/$name-out.txt  -P ${projectname};
        echo -e "SURVIVOR filter $breakseq_path/raw/$name.vcf $lowmq_path/$name.bed 300 2000000 -1 -1 $breakseq_path/lowmq/$name.lowmq.vcf" | bsub -n 1 -J ${uuid}filter_lowmq_breakseq_${name} -w 'done('$uuid'breakseq_run_'$name') && done('$uuid'generate_lowmq_'$name')' -e $reports_path/logs/$name-err.txt -o $reports_path/logs/$name-out.txt  -P ${projectname};
        echo -e "SURVIVOR filter $cnvnator_path/raw/$name.vcf $lowmq_path/$name.bed 300 2000000 -1 -1 $cnvnator_path/lowmq/$name.lowmq.vcf" | bsub -n 1 -J ${uuid}filter_lowmq_cnvnator_${name} -w 'done('$uuid'cnvnator_run_'$name') && done('$uuid'generate_lowmq_'$name')' -e $reports_path/logs/$name-err.txt -o $reports_path/logs/$name-out.txt  -P ${projectname};
        echo -e "SURVIVOR filter $delly_path/raw/$name.vcf $lowmq_path/$name.bed 300 2000000 -1 -1 $delly_path/lowmq/$name.lowmq.vcf" | bsub -n 1 -J ${uuid}filter_lowmq_delly_${name} -w 'done('$uuid'delly_run_'$name') && done('$uuid'generate_lowmq_'$name')' -e $reports_path/logs/$name-err.txt -o $reports_path/logs/$name-out.txt  -P ${projectname};
        #echo -e "SURVIVOR filter $erds_path/raw/$name.vcf $lowmq_path/$name.bed 300 2000000 -1 -1 $erds_path/lowmq/$name.lowmq.vcf" | bsub -n 1 -J ${uuid}filter_lowmq_erds_${name} -w 'done(erds_run_'$name') && done(generate_lowmq_'$name')' -e $reports_path/logs/$name-err.txt -o $reports_path/logs/$name-out.txt  -P ${projectname};
        #echo -e "SURVIVOR filter $genomestrip_path/raw/$name.vcf $lowmq_path/$name.bed 300 2000000 -1 -1 $genomestrip_path/lowmq/$name.lowmq.vcf" | bsub -n 1 -J ${uuid}filter_lowmq_genomestrip_${name} -w 'done(genomestrip_run_'$name') && done(generate_lowmq_'$name')' -e $reports_path/logs/$name-err.txt -o $reports_path/logs/$name-out.txt  -P ${projectname};
        echo -e "SURVIVOR filter $manta_path/raw/$name.vcf $lowmq_path/$name.bed 300 2000000 -1 -1 $manta_path/lowmq/$name.lowmq.vcf" | bsub -n 1 -J ${uuid}filter_lowmq_manta_${name} -w 'done('$uuid'manta_run_'$name') && done('$uuid'generate_lowmq_'$name')' -e $reports_path/logs/$name-err.txt -J ${uuid}filter_lowmq_manta_$name -o $reports_path/logs/$name-out.txt  -P ${projectname};
        echo -e "SURVIVOR filter $speedseq_path/raw/$name.vcf $lowmq_path/$name.bed 300 2000000 -1 -1 $speedseq_path/lowmq/$name.lowmq.vcf" | bsub -n 1 -J ${uuid}filter_lowmq_speedseq_${name} -w 'done('$uuid'speedseq_run_'$name') && done('$uuid'generate_lowmq_'$name')' -e $reports_path/logs/$name-err.txt -o $reports_path/logs/$name-out.txt  -P ${projectname};
        echo -e "SURVIVOR filter $svaba_path/raw/$name.vcf $lowmq_path/$name.bed 300 2000000 -1 -1 $svaba_path/lowmq/$name.lowmq.vcf" | bsub -n 1 -J ${uuid}filter_lowmq_svaba_${name} -w 'done('$uuid'svaba_run_'$name') && done('$uuid'generate_lowmq_'$name')' -e $reports_path/logs/$name-err.txt -o $reports_path/logs/$name-out.txt  -P ${projectname};
        echo -e "SURVIVOR filter $wham_path/raw/$name.vcf $lowmq_path/$name.bed 300 2000000 -1 -1 $wham_path/lowmq/$name.lowmq.vcf" | bsub -n 1 -J ${uuid}filter_lowmq_wham_${name} -w 'done('$uuid'wham_run_'$name') && done('$uuid'generate_lowmq_'$name')' -e $reports_path/logs/$name-err.txt -o $reports_path/logs/$name-out.txt  -P ${projectname};

	done < "$1"
}

filter_calls_quality() {
	while IFS='\n' read -r line || [[ -n "$line" ]]; do
        fname=$(basename "$line");
        name=$(echo $fname | cut -f 1 -d '.');
        IFS='- ' read -r -a array <<< "$name";
        family="${array[0]}";
        echo -e "bcftools view -i 'INFO/num_Reads[0] > 16' $breakdancer_path/lowmq/$name.lowmq.vcf > $breakdancer_path/filtered/$name.filtered.vcf" | bsub -n 1 -J ${uuid}filter_quality_breakdancer_${name} -w 'done('$uuid'filter_lowmq_breakdancer_'$name')' -e $reports_path/logs/$name-err.txt -o $reports_path/logs/$name-out.txt  -P ${projectname};

        echo -e "bcftools view -i 'FORMAT/PE[0] > 3' $breakseq_path/lowmq/$name.lowmq.vcf > $breakseq_path/filtered/$name.filtered.vcf" | bsub -n 1 -J ${uuid}filter_quality_breakseq_${name} -w 'done('$uuid'filter_lowmq_breakseq_'$name')' -e $reports_path/logs/$name-err.txt -o $reports_path/logs/$name-out.txt  -P ${projectname};

        echo -e "bcftools view -i 'INFO/natorQ0[0] > 0.5' $cnvnator_path/lowmq/$name.lowmq.vcf > $cnvnator_path/filtered/$name.filtered.vcf" | bsub -n 1 -J ${uuid}filter_quality_cnvnator_${name} -w 'done('$uuid'filter_lowmq_cnvnator_'$name')' -e $reports_path/logs/$name-err.txt -o $reports_path/logs/$name-out.txt  -P ${projectname};

        echo -e "bcftools view -i 'INFO/PE[0] > 10 | INFO/SR[0] > 13' $delly_path/lowmq/$name.lowmq.vcf > $delly_path/filtered/$name.filtered.vcf" | bsub -n 1 -J ${uuid}filter_quality_delly_${name} -w 'done('$uuid'filter_lowmq_delly_'$name')' -e $reports_path/logs/$name-err.txt -o $reports_path/logs/$name-out.txt  -P ${projectname};

        #echo -e "bcftools view $erds_path/lowmq/$name.lowmq.vcf > $erds_path/filtered/$name.filtered.vcf" | bsub -n 1 -J ${uuid}filter_quality_erds_${name} -w 'done(filter_lowmq_erds_'$name')' -e $reports_path/logs/$name-err.txt -o $reports_path/logs/$name-out.txt  -P ${projectname};

        #echo -e "bcftools view  $genomestrip_path/lowmq/$name.lowmq.vcf > $genomestrip_path/filtered/$name.filtered.vcf" | bsub -n 1 -J ${uuid}filter_quality_genomestrip_${name} -w 'done(filter_lowmq_genomestrip_'$name')' -e $reports_path/logs/$name-err.txt -o $reports_path/logs/$name-out.txt  -P ${projectname};

        echo -e "bcftools view -i 'FORMAT/PR[0] > 18 | FORMAT/SR[0] > 18' $manta_path/lowmq/$name.lowmq.vcf > $manta_path/filtered/$name.filtered.vcf" | bsub -n 1 -J ${uuid}filter_quality_manta_${name} -w 'done('$uuid'filter_lowmq_manta_'$name')' -e $reports_path/logs/$name-err.txt -o $reports_path/logs/$name-out.txt  -P ${projectname};
        
        echo -e "bcftools view -i 'INFO/PE[0] > 12 | INFO/SR[0] > 12' $speedseq_path/lowmq/$name.lowmq.vcf > $speedseq_path/filtered/$name.filtered.vcf" | bsub -n 1 -J ${uuid}filter_quality_speedseq_${name} -w 'done('$uuid'filter_lowmq_speedseq_'$name')' -e $reports_path/logs/$name-err.txt -o $reports_path/logs/$name-out.txt -P ${projectname};
        
        echo -e "bcftools view -i 'FORMAT/AD[0] > 12' $svaba_path/lowmq/$name.lowmq.vcf > $svaba_path/filtered/$name.filtered.vcf" | bsub -n 1 -J ${uuid}filter_quality_svaba_${name} -w 'done('$uuid'filter_lowmq_svaba_'$name')' -e $reports_path/logs/$name-err.txt -o $reports_path/logs/$name-out.txt  -P ${projectname};
        
        echo -e "bcftools view -i 'INFO/A[0] > 8' $wham_path/lowmq/$name.lowmq.vcf > $wham_path/filtered/$name.filtered.vcf" | bsub -n 1 -J ${uuid}filter_quality_wham_${name} -w 'done('$uuid'filter_lowmq_wham_'$name')' -e $reports_path/logs/$name-err.txt -o $reports_path/logs/$name-out.txt  -P ${projectname};

	done < "$1"	
}

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
		#echo "$erds_path/filtered/$name.filtered.vcf" >> $reports_path/sample_lists/$name.filtered.list
		#echo "$genomestrip_path/filtered/$name.filtered.vcf" >> $reports_path/sample_lists/$name.filter#$ed.list
		echo "$manta_path/filtered/$name.filtered.vcf" >> $reports_path/sample_lists/$name.filtered.list
		echo "$speedseq_path/filtered/$name.filtered.vcf" >> $reports_path/sample_lists/$name.filtered.list
		echo "$svaba_path/filtered/$name.filtered.vcf" >> $reports_path/sample_lists/$name.filtered.list
		echo "$wham_path/filtered/$name.filtered.vcf" >> $reports_path/sample_lists/$name.filtered.list

		echo "SURVIVOR merge $reports_path/sample_lists/$name.filtered.list 200 2 1 0 0 300 $reports_path/sample_final/$name.vcf" |bsub -n 1 -J ${uuid}merge_samples_${name} -w 'done('$uuid'filter_quality*)' -e $reports_path/logs/$name-err.txt -o $reports_path/logs/$name-out.txt -P ${projectname};

 	done < "$1"
}

genotype_samples() {
	while IFS='' read -r line || [[ -n "$line" ]]; do
	echo "Text read from file: $line"
        fname=$(basename "$line");
        echo "$fname";
        name=$(echo $fname | cut -f 1 -d '.');
        IFS='- ' read -r -a array <<< "-t /";
        family="${array[0]}";
		echo "svtyper -i $reports_path/population_final/$projectname.vcf -B $line > $reports_path/population_final/tmp/$name.genotyped.vcf" | bsub -J ${uuid}genotype_samples_$name -w 'done('$uuid'generate_population_vcf)' -n 2 -e $reports_path/logs/$name-err.txt -o $reports_path/logs/$name-out.txt -P ${projectname};
	done < "$1"
}

generate_population_vcf() {

	echo "ls $reports_path/sample_final/*.vcf > $reports_path/sample_final/sample.list;SURVIVOR merge $reports_path/sample_final/sample.list 1000 1 1 0 0 300 $reports_path/population_final/$projectname.survivor.vcf;cut -d$'\t' -f1-9 $reports_path/population_final/$projectname.survivor.vcf > $reports_path/population_final/$projectname.vcf;rm -rf $reports_path/population_final/tmp;mkdir $reports_path/population_final/tmp;" |bsub -n 1 -J ${uuid}generate_population_vcf  -w 'done('$uuid'merge_samples_*)' -e $reports_path/logs/$name-err.txt -o $reports_path/logs/$name-out.txt -P ${projectname};
}

generate_genotyped_population_vcf() {
	echo "ls $reports_path/population_final/tmp/*.genotyped.vcf > $reports_path/population_final/tmp/genotype.list;java -jar -Xmx64G $software_path/overlapMonster-1.0-jar-with-dependencies.jar -mode annotate -tool genotype -i $reports_path/population_final/$projectname.vcf -g $reports_path/population_final/tmp/genotype.list -outputPath $reports_path/population_final/$projectname.gt.vcf;" | bsub -J ${uuid}generate_genotyped_population_vcf -w 'done('$uuid'genotype_samples_*)'  -n 16 -e $reports_path/logs/$name-err.txt -o $reports_path/logs/$name-out.txt -P ${projectname};
}

generate_annotations() {
	echo "$software_path/annotate.sh $work_dir $projectname" | -J ${uuid}generate_annotions -w 'done('$uuid'generate_genotyped_population_vcf*)' bsub -n 16 -e $reports_path/logs/$name-err.txt -o $reports_path/logs/$name-out.txt -P ${projectname};

}

generate_dirs() {
	mkdir -p $breakdancer_path/{raw,lowmq,filtered}
	mkdir -p $breakseq_path/{raw,lowmq,filtered}
	mkdir -p $cnvnator_path/{raw,lowmq,filtered}
	mkdir -p $delly_path/{raw,lowmq,filtered}
	mkdir -p $erds_path/{raw,lowmq,filtered}
	mkdir -p $genomestrip_path/{raw,lowmq,filtered}
	mkdir -p $lowmq_path
	mkdir -p $manta_path/{raw,lowmq,filtered}
	mkdir -p $speedseq_path/{raw,lowmq,filtered,discordants,splitters}
	mkdir -p $svaba_path/{raw,lowmq,filtered}
	mkdir -p $wham_path/{raw,lowmq,filtered}
	mkdir -p $reports_path/{annotations,logs,population_final,sample_final,sample_lists,statistics}
}

run_breakdancer() {
	while IFS='' read -r line || [[ -n "$line" ]]; do
	echo "Text read from file: $line"
        fname=$(basename "$line");
        echo "$fname";
        name=$(echo $fname | cut -f 1 -d '.');
        echo "$name";
        IFS='- ' read -r -a array <<< "$name";
        family="${array[0]}";
 		echo -ne "perl /gpfs/software/RHEL7/genomics/breakdancer/lib/bam2cfg.pl $line | tr -d '\n\r' > $breakdancer_path/$name.cfg;breakdancer-max $breakdancer_path/$name.cfg > $breakdancer_path/$name.tsv;cat $breakdancer_path/$name.tsv | python $software_path/breakdancer2vcf.py | grep -E 'Type=ITX|Type=CTX|Type=BND|Type=TRA|hs37d5|GL0|MT' -v > $breakdancer_path/raw/$name.vcf" | bsub -n 1 -J ${uuid}breakdancer_run_${name} -e ${breakdancer_path}/${name}-err.txt -o ${breakdancer_path}/${name}-out.txt -R span[hosts=1] -P ${projectname}
	done < "$1"
}

run_breakseq2() {
	while IFS='' read -r line || [[ -n "$line" ]]; do
	echo "Text read from file: $line"
		fname=$(basename "$line");
		echo "$fname";
		name=$(echo $fname | cut -f 1 -d '.');
		echo "$name";

		echo "run_breakseq2.py --reference $ref_data --bams $line --work $breakseq_path/$name --bwa /gpfs/software/genomics/BWA/0.7.15/bwa --samtools /gpfs/software/genomics/SAMtools/1.3/bin/samtools --bplib_gff $software_path/breakseq2_bplib_20150129.gff --nthreads 2;zgrep ^# $breakseq_path/$name/breakseq.vcf.gz > $breakseq_path/raw/$name.vcf;zgrep -v ^# $breakseq_path/$name/breakseq.vcf.gz | grep -v -E \"SVTYPE=TRA|SVTYPE=BND|hs37d5|GL0|MT\" >> $breakseq_path/raw/$name.vcf" | bsub -n 2 -J ${uuid}breakseq_run_${name} -e ${breakseq_path}/${name}-err.txt -o ${breakseq_path}/${name}-out.txt -R span[hosts=1] -P ${projectname}
	done < $1
}

run_cnvnator() {
	while IFS='' read -r line || [[ -n "$line" ]]; do
	echo "Text read from file: $line"
		fname=$(basename "$line");
		echo "$fname";
		name=$(echo $fname | cut -f 1 -d '.');
		echo "$name";
		echo "cnvnator -root $cnvnator_path/$name.root -tree $line -unique;cnvnator -root $cnvnator_path/$name.root -his 100 -d $cnvnator_chromosomes_path;cnvnator -root $cnvnator_path/$name.root -stat 100 -d $cnvnator_chromosomes_path;cnvnator -root $cnvnator_path/$name.root -partition 100 -d $cnvnator_chromosomes_path;cnvnator -root $cnvnator_path/$name.root -call 100 > $cnvnator_path/$name.tsv;cnvnator2VCF.pl $cnvnator_path/$name.tsv > $cnvnator_path/$name.vcf;grep ^# $cnvnator_path/$name.vcf > $cnvnator_path/raw/$name.vcf;grep -v ^# $cnvnator_path/$name.vcf | grep -v -E \"SVTYPE=TRA|SVTYPE=BND|hs37d5|GL0|MT\" >> $cnvnator_path/raw/$name.vcf" | bsub -n 4 -J ${uuid}cnvnator_run_${name} -e $cnvnator_path/$name-err.txt -o $cnvnator_path/$name-out.txt -P ${projectname}
	done < $1
}

run_delly() {
	while IFS='' read -r line || [[ -n "$line" ]]; do
	echo "Text read from file: $line"
		fname=$(basename "$line");
		echo "$fname";
		name=$(echo $fname | cut -f 1 -d '.');
		echo "$name";
		echo "delly call -g $ref_data -o $delly_path/$name.bcf -x $software_path/human.hg19.excl.tsv -n -q 20 $line;bcftools view $delly_path/$name.bcf > $delly_path/$name.vcf;grep ^# $delly_path/$name.vcf > $delly_path/raw/$name.vcf;grep -v ^# $delly_path/$name.vcf | grep -v -E \"SVTYPE=TRA|SVTYPE=BND|hs37d5|GL0|MT\" >> $delly_path/raw/$name.vcf" | bsub -n 1 -J ${uuid}delly_run_${name} -e $delly_path/$name-err.txt -o $delly_path/$name-out.txt -P ${projectname}
	
	done < $1
}

run_erds() {
	while IFS='' read -r line || [[ -n "$line" ]]; do
	echo "Text read from file: $line"
        fname=$(basename "$line");
        echo "$fname";
        name=$(echo $fname | cut -f 1 -d '.');
        echo "$name";
        outdir="$erds_path/$name"
		rm -rf $outdir;
        mkdir $outdir;
        echo $outdir;
        echo "perl $software_path/erds1.1/erds_pipeline.pl  --samtools /gpfs/software/genomics/SAMtools/1.3/bin/samtools -o $outdir -b $line -v /gpfs/projects/tmedicine/KFakhroLAB/MND/Analysis/vcfFiles/$name.vcf.gz -r $ref_data -n $name" | bsub -n 2 -e $erds_path/$name.err -o $erds_path/$name.out -P ${projectname}       
        #grep "##" $outdir/$sample.erds.vcf > $outdir/$name.tmp.vcf
        #cut -d$'\t' -s -f1-8  $outdir/$sample.erds.vcf >> $outdir/$name.tmp.vcf;
        #cat $outdir/$name.tmp.vcf | python /gpfs/projects/bioinfo/najeeb/withKhalid/SVtoolsCompare/survivor/add_ciend.py 100 >> $outdir/$name.tmp.ciend.vcf;
        #echo "svtyper -i $outdir/$name.tmp.ciend.vcf -B $line --verbose > /gpfs/projects/bioinfo/najeeb/withKhalid/SVtoolsCompare/ERDS/$name.gt.vcf;" | bsub -n 1 -e err.txt -o out.txt -P ${projectname}
	done < "$1"
}

run_genomestrip() {
	#bamlist,workdir,gendermap


	echo "sh $software_path/gs_discovery.sh $bam_list $work_dir $name" | bsub -n 16 -e $genomestrip_path/discovery-err.txt -o $genomestrip_path/discovery-out.txt -R span[hosts=1] -P ${projectname}
}

run_manta() {
	while IFS='' read -r line || [[ -n "$line" ]]; do
	echo "Text read from file: $line"
		fname=$(basename "$line");
		echo "$fname";
		name=$(echo $fname | cut -f 1 -d '.');
		echo "$name";
		IFS='- ' read -r -a array <<< "$name";
		family="${array[0]}";
		rm -rf ${manta_path}/${name}/
		mkdir -p ${manta_path}/${name}/

		/gpfs/software/genomics/manta/1.4.0/bin/configManta.py --bam ${line} --referenceFasta $ref_data --runDir ${manta_path}/${name}/

		echo "python ${manta_path}/${name}/runWorkflow.py -m local -j 4;zgrep ^# $manta_path/$name/results/variants/diploidSV.vcf.gz > $manta_path/raw/$name.vcf;zgrep -v ^# $manta_path/$name/results/variants/diploidSV.vcf.gz | grep -v -E \"SVTYPE=TRA|SVTYPE=BND|hs37d5|GL0|MT\" >> $manta_path/raw/$name.vcf" | bsub -n 2 -J ${uuid}manta_run_${name} -e $manta_path/$name-err.txt -o $manta_path/$name-out.txt -R span[hosts=1] -P ${projectname}

		done < $1
}

run_svaba() {
	while IFS='' read -r line || [[ -n "$line" ]]; do
	echo "Text read from file: $line"
		fname=$(basename "$line");
		echo "$fname";
		name=$(echo $fname | cut -f 1 -d '.');
		echo "$name";
		IFS='- ' read -r -a array <<< "$name";
		family="${array[0]}";

		echo "svaba run -t $line -p 4 -L 6 -a $svaba_path/$name -G $ref_data;java -jar $software_path/htsSidra-1.0-jar-with-dependencies.jar convert $svaba_path/$name.svaba.sv.vcf;sed -i 's/\"\"/\"/g' $svaba_path/$name.converted.vcf;grep ^# $svaba_path/$name.converted.vcf > $svaba_path/raw/$name.vcf;grep -v ^# $svaba_path/$name.converted.vcf | grep -v -E \"SVTYPE=TRA|SVTYPE=BND|hs37d5|GL0|MT\" >> $svaba_path/raw/$name.vcf" | bsub -n 2 -J ${uuid}svaba_run_${name} -e ${svaba_path}/${name}-err.txt -o ${svaba_path}/${name}-out.txt -R span[hosts=1] -P ${projectname}

	done < $1
}

run_speedseq() {
	while IFS='' read -r line || [[ -n "$line" ]]; do
	echo "Text read from file: $line"
		fname=$(basename "$line");
		echo "$fname";
		name=$(echo $fname | cut -f 1 -d '.');
		echo "$name";

		##Extract discordants
		echo "samtools view -b -F 1294 $line > $speedseq_path/discordants/$name.bwa.discordants.unsorted.bam;samtools sort $speedseq_path/discordants/$name.bwa.discordants.unsorted.bam -o $speedseq_path/discordants/$name.bwa.discordants.sorted.bam;rm -f $speedseq_path/discordants/$name.bwa.discordants.unsorted.bam" | bsub -J ${uuid}discordants_$name -n 1 -e ${speedseq_path}/${name}-err.txt -o ${speedseq_path}/${name}-out.txt -R span[hosts=1] -P ${projectname}
		##Extract splitters
		echo "samtools view -h $line | /gpfs/software/genomics/lumpy/scripts/extractSplitReads_BwaMem -i stdin | samtools view -Sb - > $speedseq_path/splitters/$name.bwa.splitters.unsorted.bam; samtools sort $speedseq_path/splitters/$name.bwa.splitters.unsorted.bam -o $speedseq_path/splitters/$name.bwa.splitters.sorted.bam" | bsub -J ${uuid}splitters_$name -n 1  -e ${speedseq_path}/${name}-err.txt -o ${speedseq_path}/${name}-out.txt -R span[hosts=1] -P ${projectname}
		#Variant calling
		echo "speedseq sv -g -t 2 -K /gpfs/software/genomics/speedseq/v0.1.2/speedseq/bin/lumpyexpress.config -x $software_path/ceph18.b37.lumpy.exclude.2014-01-15.bed -o $speedseq_path/$name -R $ref_data -B $line -D $speedseq_path/discordants/$name.bwa.discordants.sorted.bam -S $speedseq_path/splitters/$name.bwa.splitters.sorted.bam;zgrep ^# $speedseq_path/$name.sv.vcf.gz > $speedseq_path/raw/$name.vcf;zgrep -v ^# $speedseq_path/$name.sv.vcf.gz | grep -v -E \"SVTYPE=TRA|SVTYPE=BND|hs37d5|GL0|MT\" >> $speedseq_path/raw/$name.vcf" | bsub -n 1  -J ${uuid}speedseq_run_${name} -w 'done('$uuid'discordants_'$name') && done('$uuid'splitters_'$name')' -e ${speedseq_path}/${name}-err.txt -o ${speedseq_path}/${name}-out.txt -R span[hosts=1] -P ${projectname}

	done < $1
}

run_wham() {
	while IFS='' read -r line || [[ -n "$line" ]]; do
	echo "Text read from file: $line"
		fname=$(basename "$line");
		echo "$fname";
		name=$(echo $fname | cut -f 1 -d '.');
		echo "$name";
		IFS='- ' read -r -a array <<< "$name";
		family="${array[0]}";

		echo "whamg -x 4 -e GL000207.1,GL000226.1,GL000229.1,GL000231.1,GL000210.1,GL000239.1,GL000235.1,GL000201.1,GL000247.1,GL000245.1,GL000197.1,GL000203.1,GL000246.1,GL000249.1,GL000196.1,GL000248.1,GL000244.1,GL000238.1,GL000202.1,GL000234.1,GL000232.1,GL000206.1,GL000240.1,GL000236.1,GL000241.1,GL000243.1,GL000242.1,GL000230.1,GL000237.1,GL000233.1,GL000204.1,GL000198.1,GL000208.1,GL000191.1,GL000227.1,GL000228.1,GL000214.1,GL000221.1,GL000209.1,GL000218.1,GL000220.1,GL000213.1,GL000211.1,GL000199.1,GL000217.1,GL000216.1,GL000215.1,GL000205.1,GL000219.1,GL000224.1,GL000223.1,GL000195.1,GL000212.1,GL000222.1,GL000200.1,GL000193.1,GL000194.1,GL000225.1,GL000192.1,NC_007605,hs37d5 -a $ref_data -f $line > $wham_path/$name.vcf;grep ^# $wham_path/$name.vcf > $wham_path/raw/$name.vcf;grep -v ^# $wham_path/$name.vcf | grep -v -E \"SVTYPE=TRA|SVTYPE=BND|hs37d5|GL0|MT\" >> $wham_path/raw/$name.vcf" | bsub -n 4 -J ${uuid}wham_run_${name} -e $wham_path/$name-err.txt -o $wham_path/$name-out.txt -P ${projectname}
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

generate_dirs

generate_lowmq_bed $bam_list


run_breakdancer $bam_list
run_breakseq2 $bam_list
run_cnvnator $bam_list
run_delly $bam_list
#run_erds $bam_list
#run_genomestrip bam_list
run_manta $bam_list
run_svaba $bam_list
module unload python/2.7
run_speedseq $bam_list
module load python/2.7
run_wham $bam_list

filter_calls_lowmq $bam_list
filter_calls_quality $bam_list

merge_consensuns_sample $bam_list

module unload speedseq/0.1.2
module unload python/2.7

module load svtyper/v0.7.0

generate_population_vcf $bam_list

genotype_samples $bam_list

generate_genotyped_population_vcf 

generate_annotions

echo "SIDRA CV CALLING PIPELINE FINISHED";



