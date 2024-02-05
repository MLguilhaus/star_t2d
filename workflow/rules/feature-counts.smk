def get_gtf(wildcards):
    pat = os.path.join(ref_path, "*.gtf.gz")
    gtf = glob.glob(pat)[0]
    return(gtf)

rule feature_counts_star:
    input:
        bam = expand(
            os.path.join(star_outpath, "{sample}", "{bam}"),
            bam = ["Aligned.sortedByCoord.out.bam"],
            # SP said something here about lane samples, change to acession?
            sample = accessions
        ),
        bai = expand(
            os.path.join(star_outpath, "{sample}", "{bai}"),
            bai = ["Aligned.sortedByCoord.out.bam.bai"],
            sample = accessions
        ),
        gtf = get_gtf
    output:
        counts = os.path.join(fcount_path, "star_counts.out"),
        summary = os.path.join(fcount_path, "star_counts.out.summary")
    conda: "../envs/subread.yml"
    log: os.path.join(log_path, "feature_counts", "star_feature_counts.log")
    threads: 12
    resources:
        runtime = "6:00:00",
        mem_mb = 8192
    params:
        minOverlap = config['featureCounts']['minOverlap'],
        fracOverlap = config['featureCounts']['fracOverlap'],
        q = config['featureCounts']['minQual'],
        s = config['featureCounts']['strandedness'],
        extra = config['featureCounts']['extra']
    shell:
       """
       featureCounts \
         {params.extra} \
         -Q {params.q} \
         -s {params.s} \
         --minOverlap {params.minOverlap} \
         --fracOverlap {params.fracOverlap} \
         -T {threads} \
         -a {input.gtf} \
         -o {output.counts} \
         {input.bam} &>> {log}
       """