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