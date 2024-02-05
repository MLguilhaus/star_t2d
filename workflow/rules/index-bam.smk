rule index_bam:
    input:
        bam = expand(
            os.path.join(star_outpath, "{sample}", "{bam}"),
            bam = ["Aligned.sortedByCoord.out.bam"],
            sample = accessions
            ),
    output:
        bai = expand(
            os.path.join(star_outpath, "{sample}", "{bai}"),
            bai = ["Aligned.sortedByCoord.out.bam.bai"],
            sample = accessions
        ),
    params:
        extra = config['samtools']['extra']
    shell:
        """samtools index {input.bam} {params.extra} -o {output.bai}"""