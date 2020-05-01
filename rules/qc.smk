rule download_adapters:
    output:
        temp(directory("BBMap")),
        "data/qc/adapters.fa"
    shell: "git clone --depth 1 https://github.com/BioInfoTools/BBMap.git && cp BBMap/resources/adapters.fa data/qc"

rule skewer:
    input:
        fq=FQ_FOLDER + "{fq}.fastq.gz",
        adapters="data/qc/adapters.fa"
    output:
        fq1="data/trimmed/{fq}.trim_1.fastq.gz",
        fq2="data/trimmed/{fq}.trim_2.fastq.gz",
    threads: 12
    log: "data/trimmed/{fq}-trimmed.status"
    params:
        quality=20,
        ext="data/trimmed/{fq}"
    shell:
        "skewer -q {params.quality} -o {params.ext}"
        " -t {threads} -x {input.adapters} {input.fq1} {input.fq2} &> {log} && "
        "mv {params.ext}-trimmed-pair1.fastq {params.ext}.trim_1.fastq && "
        "mv {params.ext}-trimmed-pair2.fastq {params.ext}.trim_2.fastq && "
        "gzip {params.ext}.trim_*.fastq"

rule fastqc_analysis:
    input:
        fq1=["data/{fq}_1.fastq.gz", "data/trimmed/{fq}.trim_1.fastq.gz"],
        fq2=["data/{fq}_2.fastq.gz", "data/trimmed/{fq}.trim_2.fastq.gz"],
    output:
        fq1=["data/{fq}_1_fastqc.html", "data/{fq}_1_fastqc.zip", "data/trimmed/{fq}.trim_1_fastqc.html", "data/trimmed/{fq}.trim_1_fastqc.zip"],
        fq2=["data/{fq}_2_fastqc.html", "data/{fq}_2_fastqc.zip", "data/trimmed/{fq}.trim_2_fastqc.html", "data/trimmed/{fq}.trim_2_fastqc.zip"],
    log: "data/{fq}.fastqc.log",
    threads: 6
    shell:
        "fastqc -t {threads} {input.fq1} {input.fq2} 2> {log}"
