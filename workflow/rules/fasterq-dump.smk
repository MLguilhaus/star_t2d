rule get_fastq:
    output:
        os.path.join(raw_path, "{accession}_1.fastq.gz"),
        os.path.join(raw_path, "{accession}_2.fastq.gz"),
    log: "workflow/logs/fasterq-dump/{accession}.log"
    conda: "../envs/fasterq-dump.yml"
    params:
        extra="--skip-technical",
    threads: 2
    retries: 3
    resources:
        runtime = "2h"    
    script:
        "../scripts/fasterq-dump.py" # Taken from the wrapper

rule raw_md5sums:
    input: 
        expand(
            os.path.join(raw_path, "{sample}{tag}.fastq.gz"),
            sample = accessions, 
            tag = ['_1', '_2']
        )
    output: 
        os.path.join(raw_path, "md5sums.txt")
    threads: 1
    resources:
        runtime="10m"
    shell:
        """
        md5sum {input} > {output}
        """