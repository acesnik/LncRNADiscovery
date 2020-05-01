RAM_MB_REQ = 50000 #mb
RAM_B_REQ = RAM_MB_REQ * 1000 # BYTES

rule star_genome_generate:
    input:
        efa="data/ERCC.fa",
        gfa=FA,
        gff=GFF3 + ".fix.gff3"
    output: STAR_REF_FOLDER + "/SA"
    params: genomedir=STAR_REF_FOLDER
    threads: 99
    shell:
        "STAR --runMode genomeGenerate --runThreadN {threads} --genomeDir {params.genomedir} "
        "--genomeFastaFiles {input.gfa} {input.efa} --sjdbGTFfile {input.gff} --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 100"

rule load_star_genome_firstpass:
    input: STAR_REF_FOLDER + "/SA"
    output: temp("output/loaded_firstpass")
    params: genomedir=STAR_REF_FOLDER
    resources: mem_mb = RAM_MB_REQ
    shell: "STAR --genomeLoad LoadAndExit --genomeDir {params} && touch {output}"

rule fastqc_analysis:
    input: FQ_FOLDER + "{fq}.fastq.gz", #single end only
    output: FQ_FOLDER + "{fq}_fastqc.html"
    log: "output/{fq}.fastqc.log"
    threads: 6
    shell: "fastqc -t {threads} {input} 2> {log}"

rule star_firstpass:
    input:
        temp("output/loaded_firstpass"),
        suffix=STAR_REF_FOLDER + "/SA",
        fastq=FQ_FOLDER + "{fq}.fastq.gz", #single end only
        fastqc=FQ_FOLDER + "{fq}_fastqc.html" # trigger qc analysis
    output: "output/SJ1st/{fq}SJ.out.tab"
    threads: 6
    params:
        junctions="--outFilterIntronMotifs RemoveNoncanonical", # adds XS tag to all alignments that contain a splice junction
        bam="--outSAMtype None",
        sjfilter=" --outSJfilterReads Unique", # for JUM
        genomedir=STAR_REF_FOLDER
    shell:
        "STAR --genomeLoad LoadAndKeep --genomeDir {params.genomedir}"
        " --runThreadN {threads} {params.junctions} {params.bam} {params.sjfilter}"
        " --outFileNamePrefix output/SJ1st/{wildcards.fq}"
        " --readFilesIn <(zcat " + FQ_FOLDER + "{wildcards.fq}.fastq.gz)"

rule unload_firstpass_genome:
    input:
        STAR_REF_FOLDER + "/SA",
        jj=expand("output/SJ1st/{fq}SJ.out.tab", fq=FQ_PREFIXES)
    output: temp("temp/unloaded_firstpass")
    params: genomedir=STAR_REF_FOLDER
    shell:
        "STAR --genomeLoad Remove --genomeDir {params.genomedir} && "
        "touch {output}"

rule star_genome_generate_secondpass:
    input:
        efa="data/ERCC.fa",
        gfa=FA,
        gff=GFF3 + ".fix.gff3",
        jj=expand("output/SJ1st/{fq}SJ.out.tab", fq=FQ_PREFIXES)
    output: STAR_REF_FOLDER + unique_tag() + "/SA"
    params:
        genomedir=STAR_REF_FOLDER + unique_tag()
    threads: 99
    shell:
        "STAR --runMode genomeGenerate --runThreadN {threads} --genomeDir {params.genomedir} "
        "--genomeFastaFiles {input.gfa} {input.efa} --sjdbGTFfile {input.gff} --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 100"
        "--limitSjdbInsertNsj 1200000 --sjdbFileChrStartEnd {input.jj}"

rule load_star_genome_2pass:
    input: STAR_REF_FOLDER + unique_tag() + "/SA"
    output: temp("output/loaded_2pass")
    params: genomedir=STAR_REF_FOLDER + unique_tag()
    resources: mem_mb = RAM_MB_REQ
    shell: "STAR --genomeLoad LoadAndExit --genomeDir {params} && touch {output}"

rule star_2pass:
    input:
        temp("output/loaded_2pass"),
        genomedir=STAR_REF_FOLDER + unique_tag() + "/SA",
        fastq=FQ_FOLDER + "{fq}.fastq.gz" #single end only
    output:
        sj="output/{fq}" + unique_tag() + STAR_OUT_SJ_SUFFIX_2ND,
        bam="output/{fq}" + unique_tag() + STAR_OUT_BAM_SUFFIX_2ND,
        log="output/{fq}" + unique_tag() + STAR_OUT_LOG_SUFFIX_2ND,
        final="output/{fq}" + unique_tag() + STAR_OUT_FIN_SUFFIX_2ND,
        progress="output/{fq}" + unique_tag() + STAR_OUT_PROG_SUFFIX_2ND,
    threads: 6
    params:
        junctions="--outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical", # adds XS tag to all alignments that contain a splice junction
        bam="--outSAMtype BAM SortedByCoordinate --outBAMcompression 10 --limitBAMsortRAM " + str(RAM_B_REQ),
        gatk="--outSAMattrRGline ID:1 PU:platform  PL:illumina SM:sample LB:library --outSAMmapqUnique 60",
        sjfilter=" --outSJfilterReads Unique", # for JUM
        genomedir=STAR_REF_FOLDER + unique_tag()
    shell:
        "STAR --runMode alignReads --genomeLoad LoadAndKeep --genomeDir {params.genomedir} {params.sjfilter}"
        " --runThreadN {threads} {params.junctions} {params.bam} {params.gatk}"
        " --outFileNamePrefix output/{wildcards.fq}" + unique_tag() +
        " --readFilesIn <(zcat " + FQ_FOLDER + "{wildcards.fq}.fastq.gz)"

rule unload_2pass_genome:
    input:
        STAR_REF_FOLDER + unique_tag() + "/SA",
        bam="output/{fq}" + unique_tag() + STAR_OUT_BAM_SUFFIX_2ND
    output: "temp/unloaded_2pass"
    params: genomedir=STAR_REF_FOLDER + unique_tag()
    shell:
        "STAR --genomeLoad Remove --genomeDir {params.genomedir} && "
        "touch {output}"

rule merge_bams:
    input:
        bams=expand("output/{fq}" + unique_tag() + STAR_OUT_BAM_SUFFIX_2ND, fq=FQ_PREFIXES)
    output:
        sorted="output/combined.sorted.bam",
        stats="output/combined.sorted.stats"
    params:
        compression="9",
        tempprefix="output/combined.sorted"
    log: "output/combined.sorted.log"
    threads: 12
    resources: mem_mb=16000
    shell:
        "(ls {input.bams} | "
        "{{ read firstbam; "
        "samtools view -h ""$firstbam""; "
        "while read bam; do samtools view ""$bam""; done; }} | "
        "samtools view -ubS - | "
        "samtools sort -@ {threads} -l {params.compression} -T {params.tempprefix} -o {output.sorted} - && "
        "samtools index {output.sorted} && "
        "samtools flagstat -@ {threads} {output.sorted} > {output.stats}) 2> {log}"
