rule index_bam:
    input:
         bam = os.path.join(star_outpath, "{sample}", 
            "Aligned.sortedByCoord.out.bam"),
    output:
        bai = os.path.join(star_outpath, "{sample}",
            "Aligned.sortedByCoord.out.bam.bai"),  
    params:
        extra = config['samtools']['extra']
    shell:
        """samtools index {input.bam} {params.extra} -o {output.bai}"""