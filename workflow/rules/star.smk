rule star_align:
    input:
        r1 = rules.fastp.output.fq1,
        r2 = rules.fastp.output.fq2,
    output:
        bam = os.path.join(
            star_outpath, "{accession}", 
            "Aligned.sortedByCoord.out.bam"
        ),
        star_log = os.path.join(
            star_outpath, "{accession}", "Log.final.out"
        ),
        sj_out = os.path.join(
            star_outpath, "{accession}", "SJ.out.tab"
        ),
        temp_log = os.path.join(
            star_outpath, "{accession}", "Log.out"
        )
    conda: "../envs/star.yml"
    log: os.path.join(log_path, "star_align", "{accession}.log") #star_align defined where? 
    params:
        index = config['star']['star_ref_path'],
        out_prefix = os.path.join(star_outpath, "{accession}", ""),
        out_samtype = config['star']['out_samtype'],
        out_samattributes = config['star']['out_samattributes'],
        out_samunmapped = config['star']['out_samunmapped'],
        multimap_nmax = config['star']['multimap_nmax'],
        mismatch_nmax = config['star']['mismatch_nmax'],
        extra = config['star']['extra'],
        temp_dir = os.path.join(
            star_outpath, "{accession}", "_STARtmp"
        )
    threads: 8
    resources:
        runtime = "3h",
        mem_mb = 32768
    shell:
        """
        STAR \
          --runThreadN {threads} \
          --genomeDir {params.index} \
          ### no params in below arg atm, clarify what was there with SP. 
          ### ok to leave blank
          --{params.extra} \
          --readFilesIn {input.r1} {input.r2} \
          --readFilesCommand zcat \
          --outSAMtype {params.out_samtype} \
          --outSAMattributes {params.out_samattributes} \
          --outSAMunmapped {params.out_samunmapped} \
          --outFilterMultimapNmax {params.multimap_nmax} \
          --outFilterMismatchNmax {params.mismatch_nmax} \
          --outFileNamePrefix {params.out_prefix} \
          --outStd Log &>> {log}
        
        ### Also to do with SP bug fix, clairfy. 
        ## Deleting this here ensures it is **retained on failure**, but
        ## removed on success as it is entirely redundant at that point
        echo -e "Deleting Log.progress.out"
        rm -f {params.out_prefix}/Log.progress.out

        ## Cleanup the temp_dir
        if [[ -d {params.temp_dir} ]]; then
            echo -e "Deleting {params.temp_dir}" >> {log}
            rm -rf {params.temp_dir}
        fi
        """

### Again, a SP addition to do with git logs, 
### clairfy what they meant by this (should have written it down, 
### remember to write down what stevie says/keep a notebook hand when given instruction) 
#### useful to have with git, more notes in physical notes
rule copy_star_logs:
    input:
        star_log = rules.star_align.output.star_log,
        sj_out = rules.star_align.output.sj_out
    output:
        star_log = os.path.join(
            starlog_path, "{build}", "{accession}_Log.final.out"
        ),
        sj_out = os.path.join(
            starlog_path, "{build}", "{accession}_SJ.out.tab"
        )
    params:
        dir = os.path.dirname(
            os.path.join(starlog_path, "{build}", "{accession}")
        )
    threads: 1
    resources:
        runtime = "1m"
    shell:
        """
        if [[ ! -d {params.dir} ]]; then
            mkdir -p {params.dir}
        fi
        cp {input.star_log} {output.star_log}
        cp {input.sj_out} {output.sj_out}
        """

### Gives stats on primary mapped, secondary mapped, duplicates etc
### May not need this either
### may as well keep, not very intense computationally, may as well haved
rule star_flagstat:
    input: 
        bam = os.path.join(
            star_outpath, "{build}", "{accession}",
            "Aligned.sortedByCoord.out.bam"
        ),
        bai = os.path.join(
            star_outpath, "{build}", "{accession}", 
            "Aligned.sortedByCoord.out.bam.bai"
        )
    output:
        os.path.join(
            starlog_path, "{build}", "{accession}_samtools.flagstat"
        )
    conda: "../envs/samtools.yml"
    threads: 4
    resources:
        runtime = "30m"
    shell:
        """
        samtools flagstat -@ {threads} {input.bam} > {output}
        """

rule star_md5sums:
    input:
        expand(
            os.path.join(
                star_outpath, "{{build}}", "{accession}",
                "Aligned.sortedByCoord.out.bam"
            ),
            accession = accession_table['accession']
        )
    output:
        os.path.join(star_outpath, "{build}", "md5sums.txt")
    threads: 1
    resources:
        runtime="10m"
    shell:
        """
        md5sum {input} > {output}
        """