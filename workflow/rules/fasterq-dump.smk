rule get_fastq:
    output:
        os.path.join(raw_path, "{accession}.fastq.gz")
    log: "workflow/logs/fasterq-dump/{accession}.fastq.gz.log"
    conda: "../envs/fasterq-dump.yml"
    params:
        extra="--skip-technical",
    threads: 2
    retries: 3
    resources:
        runtime = "2h"    
    script:
        "../scripts/fasterq-dump.py" # Taken from the wrapper