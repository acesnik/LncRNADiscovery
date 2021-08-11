REFSTAR_PREFIX = f"ensembl/Homo_sapiens.{GENEMODEL_VERSION}RsemStar/RsemStarReference"
REFSTAR_FOLDER = f"ensembl/Homo_sapiens.{GENEMODEL_VERSION}RsemStar/"

rule rsem_star_genome:
    '''Create an RSEM reference with STAR indices'''
    input:
        efa="data/ERCC.fa",
        gfa=FA,
        customGtf="output/combined.filtered.withcds.gtf"
    output:
        suffix = REFSTAR_FOLDER + "SA"
    threads: 99
    resources: mem_mb=60000
    log: "output/ensembl/prepare-reference.log"
    shell:
        "(rsem-prepare-reference --num-threads {threads} --star --gtf {input.customGtf} \"{input.efa}\",\"{input.gfa}\" " + REFSTAR_PREFIX +
        ") 2> {log}"

rule rsem_star_align:
    '''Align to transcripts with STAR and quantify with RSEM'''
    input:
        suffix=REFSTAR_FOLDER + "SA",
        fastq=FQ_FOLDER + "{fq}.fastq.gz", #single end only
    output:
        "output/{fq}.isoforms.results",
        "output/{fq}.genes.results",
        "output/{fq}.time",
        directory("output/{fq}.stat"),
    resources: mem_mb=50000
    threads: 12
    log: "output/{fq}calculate-expression.log"
    shell:
        "(rsem-calculate-expression --no-bam-output --time --star" # --calc-ci" not doing credibility intervals for now; they take a long time to calc.
        " --num-threads {threads} <(zcat " + FQ_FOLDER + "{wildcards.fq}.fastq.gz) " + REFSTAR_PREFIX + " output/{wildcards.fq}) &> {log}"

rule make_rsem_dataframe:
    '''Take the results from RSEM and put them in a usable dataframe'''
    input:
        expand("output/{fq}.genes.results", fq=FQ_PREFIXES),
        gff=GFF3 + ".fix.gff3"
    output:
        counts="output/Counts.csv",
        tpms="output/Tpms.csv"
    params:
        names="output/IdsToNames.csv",
    shell:
        "python scripts/make_rsem_dataframe.py {input.gff} {output.counts} {output.tpms} {params.names}"
