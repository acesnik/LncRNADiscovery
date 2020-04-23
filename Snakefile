GENOME_VERSION = "GRCh38"
ENSEMBL_VERSION = "96"
GENEMODEL_VERSION = GENOME_VERSION + "." + ENSEMBL_VERSION
GENOME_FA = f"ensembl/Homo_sapiens.{GENOME_VERSION}.dna.primary_assembly.fa"
ENSEMBL_GFF = f"ensembl/Homo_sapiens.{GENEMODEL_VERSION}.gff3"
ENSEMBL_GTF = f"ensembl/Homo_sapiens.{GENEMODEL_VERSION}.gtf"
TEST_GENOME_FA = f"ensembl/202122.fa"
TEST_ENSEMBL_GFF = f"ensembl/202122.gff3"
FA=GENOME_FA # for analysis; can also be TEST_GENOME_FA
GFF3=ENSEMBL_GFF # for analysis; can also be TEST_ENSEMBL_GFF
REFSTAR_PREFIX = f"ensembl/Homo_sapiens.{GENEMODEL_VERSION}RsemStar/RsemStarReference"
REFSTAR_FOLDER = f"ensembl/Homo_sapiens.{GENEMODEL_VERSION}RsemStar/"
REF_PREFIX = f"ensembl/Homo_sapiens.{GENEMODEL_VERSION}Rsem/RsemReference"
REF_FOLDER = f"ensembl/Homo_sapiens.{GENEMODEL_VERSION}Rsem/"
STAR_REF_FOLDER = f"ensembl/Homo_sapiens.{GENEMODEL_VERSION}Star" # no slash

# FQ_FOLDER = "/mnt/e/ProjectsActive/tonsil/data_spritzsnake/trimmed/"
FQ_FOLDER = "../CellCycleQuant/ESCG_data/"
(FQ_PREFIXES,) = glob_wildcards(FQ_FOLDER + "{fq}.fastq.gz")

STAR_OUT_SJ_SUFFIX_2ND = "SJ.out.tab"
STAR_OUT_BAM_SUFFIX_2ND = "Aligned.sortedByCoord.out.bam"
STAR_OUT_FIN_SUFFIX_2ND = "Log.final.out"
STAR_OUT_LOG_SUFFIX_2ND = "Log.out"
STAR_OUT_PROG_SUFFIX_2ND = "Log.progress.out"

THREADS=24
TEST="pvalue"
CUTOFF="0.05"

def unique_tag():
    return str(len(str(expand("{fq}", fq=FQ_PREFIXES))))

rule all:
    input:
        "output/{fq}" + unique_tag() + STAR_OUT_BAM_SUFFIX_2ND + "slncky/annotated.lncs.info.txt", # lncRNA discovery
        "data/Tpms.csv", "data/Counts.csv", "data/IdsToNames.csv" # quantification

# rule clean:
#     shell: "rm -rf ensembl output"

include: "rules/align.smk"
include: "rules/isoforms.smk"
include: "rules/lncRNAs.smk"
include: "rules/qc.smk"
include: "rules/quant.smk"
include: "rules/isoforms.smk"

rule download_ensembl_genome:
    '''Download the genome with the version specified above'''
    output: gfa=GENOME_FA,
    log: "ensembl/downloads_fa.log"
    shell:
        "(wget -O - ftp://ftp.ensembl.org/pub/release-81//fasta/homo_sapiens/dna/Homo_sapiens." + GENOME_VERSION + ".dna.primary_assembly.fa.gz | "
        "gunzip -c > {output.gfa}) 2> {log}"

rule download_ensembl_genemodel:
    '''Download the gene model with the version specified above'''
    output:
        gff=ENSEMBL_GFF,
        gtf=ENSEMBL_GTF
    log: "ensembl/download_gff.log"
    shell:
        "(wget -O - ftp://ftp.ensembl.org/pub/release-" + ENSEMBL_VERSION + "/gff3/homo_sapiens/Homo_sapiens." + GENEMODEL_VERSION + ".gff3.gz | "
        "gunzip -c > {output.gff} && "
        "wget -O - ftp://ftp.ensembl.org/pub/release-" + ENSEMBL_VERSION + "/gtf/homo_sapiens/Homo_sapiens." + GENEMODEL_VERSION + ".gtf.gz | "
        "gunzip -c > {output.gtf}) 2> {log}"

rule fix_gff3_for_rsem:
    '''This script changes descriptive notes in column 4 to "gene" if a gene row, and it also adds ERCCs to the gene model'''
    input: ENSEMBL_GFF
    output: ENSEMBL_GFF + ".fix.gff3"
    shell: "python scripts/fix_gff3_for_rsem.py {input} {output}"

rule filter_gff3:
    '''For testing, make a smaller gene model'''
    input: ENSEMBL_GFF + ".fix.gff3"
    output: TEST_ENSEMBL_GFF + ".fix.gff3"
    shell: "grep \"^ERCC\\|^#\\|^20\\|^21\\|^22\" {input} > {output}"

rule filter_fa:
    '''For testing, make a smaller genome'''
    input: GENOME_FA
    output: TEST_GENOME_FA
    script: "scripts/filter_fasta.py"
