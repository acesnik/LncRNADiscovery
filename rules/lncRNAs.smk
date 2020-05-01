SLNCKY_ANNOTATIONS = "https://www.dropbox.com/s/jifnrr3kjq95w64/annotations.tar.gz"

rule install_slncky:
    output: "tools/slncky/slncky.v1.0"
    log: "logs/install_slncky.log"
    shell: "git clone https://github.com/slncky/slncky.git tools/slncky &> {log}"

rule download_annotations:
    input: "tools/slncky/slncky.v1.0"
    output: "tools/slncky/annotations/hg38.fa"
    log: "logs/download_annotations.log"
    shell: "(cd tools/slncky && wget -c {SLNCKY_ANNOTATIONS} -O - | tar -xz) &> {log}"

rule index_fa:
    input: GENOME_FA
    output: GENOME_FA + ".fai"
    shell: "samtools faidx {input}"

rule convert_predicted_gtf_to_bed12:
    '''Use tools from BEDOPS to convert custom GTF to sorted BED12 file'''
    input:
        gtf="output/combined.sorted.filtered.withcds.gtf",
        fai=GENOME_FA + ".fai"
    output:
        temp("output/combined.sorted.filtered.withcds.gtf.genePred"),
        bed12=temp("output/combined.sorted.filtered.withcds.filtered.bed12"),
        bed12sorted="output/combined.sorted.filtered.withcds.filtered.sorted.bed12"
    log: "output/combined.sorted.filtered.withcds.filtered.bed12.log"
    shell:
        "(gtfToGenePred {input.gtf} {input.gtf}.genePred && "
        "genePredToBed {input.gtf}.genePred {output.bed12} && "
        "bedtools sort -faidx {input.fai} -i {output.bed12} > {output.bed12sorted}) &> {log}"

rule download_chromosome_mappings:
    output: "ChromosomeMappings/GRCh38_ensembl2UCSC.txt"
    shell: "git clone https://github.com/dpryan79/ChromosomeMappings.git"

rule convert_firstBed12_column_to_ucsc:
    '''The first column of the BED12 needs to be converted from Ensembl format (1, 2, MT, etc.) to UCSC format (chr1, chr2, chrM, etc.)'''
    input:
        "ChromosomeMappings/GRCh38_ensembl2UCSC.txt",
        bed12="output/combined.sorted.filtered.withcds.filtered.sorted.bed12"
    output: "output/combined.sorted.filtered.withcds.filtered.sorted.ucsc.bed12"
    shell: "python scripts/convert_ensembl2ucsc.py {input.bed12} {output}"

rule annotate_lncrnas:
    input:
        "tools/slncky/annotations/hg38.fa",
        slncky="tools/slncky/slncky.v1.0",
        bed12="output/combined.sorted.filtered.withcds.filtered.sorted.ucsc.bed12"
    output:
        "output/combined.slncky/annotated.canonical_to_lncs.txt",
        "output/combined.slncky/annotated.cluster_info.txt",
        "output/combined.slncky/annotated.filtered_info.txt",
        "output/combined.slncky/annotated.lncs.bed",
        "output/combined.slncky/annotated.lncs.info.txt",
        "output/combined.slncky/annotated.orfs.txt",
        "output/combined.slncky/annotated.orthologs.top.txt",
        "output/combined.slncky/annotated.orthologs.txt",
    params: reference="hg38"
    threads: 12
    shell:
        "python2 {input.slncky} --threads {threads} {input.bed12} {params.reference} output/combined.slncky/annotated"
