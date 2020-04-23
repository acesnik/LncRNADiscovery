REFSTAR_PREFIX = f"ensembl/Homo_sapiens.{GENEMODEL_VERSION}RsemStar/RsemStarReference"
REFSTAR_FOLDER = f"ensembl/Homo_sapiens.{GENEMODEL_VERSION}RsemStar/"

rule rsem_star_genome:
    '''Create an RSEM reference with STAR indices'''
    input:
        gfa=FA,
        customGtf = "data/combined.sorted.filtered.withcds.gtf"
    output:
        REFSTAR_PREFIX + ".gtf",
        suffix = REFSTAR_FOLDER + "SA"
    threads: 99
    resources: mem_mb=60000
    log: "data/ensembl/prepare-reference.log"
    shell:
        "(rsem-prepare-reference --num-threads {threads} --star --gtf {input.customGtf} \"{input.gfa}\" " + REFSTAR_PREFIX +
        ") 2> {log}"

rule rsem_star_align:
    '''Align to transcripts with STAR and quantify with RSEM'''
    input:
        suffix=REFSTAR_FOLDER + "SA",
        gtf=REFSTAR_PREFIX + ".gtf",
        # fq=expand("output/{fq}" + unique_tag() + STAR_OUT_BAM_SUFFIX_2ND, fq=FQ_PREFIXES)
        fq1="data/trimmed/{sra}.trim_1.fastq.gz" if check_sra() else "data/{fq}_1.fastq.gz",
        fq2="D/trimmed/{sra}.trim_2.fastq.gz" if check_sra() else "data/{fq}_2.fastq.gz",
    output:
        "data/{sra}.isoforms.results",
        "data/{sra}.genes.results",
        "data/{sra}.time",
        directory("data/{sra}.stat"),
    resources: mem_mb=50000
    threads: 12
    log: "data/{sra}calculate-expression.log"
    shell:
        "(rsem-calculate-expression --no-bam-output --time --star" # --calc-ci" not doing credibility intervals for now; they take a long time to calc.
        " --num-threads {threads} --paired-end <(zcat {input.fq1}) <(zcat {input.fq2}) " + REFSTAR_PREFIX + " data/{wildcards.sra}) &> {log}"

rule make_rsem_dataframe:
    '''Take the results from RSEM and put them in a usable dataframe'''
    input:
        expand("data/{sra}.genes.results", sra=config["sra"], dir=config["analysisDirectory"]),
        gff="data/ensembl/" + REF + "." + config["release"] + ".gff3" + ".fix.gff3"
    output:
        counts="data/Counts.csv",
        names="data/IdsToNames.csv",
        tpms="data/Tpms.csv"
    shell:
        "python scripts/make_rsem_dataframe.py {input.gff} {output.counts} {output.tpms} {output.names}"
