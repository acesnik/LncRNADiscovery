SLNCKY_ANNOTATIONS = "https://www.dropbox.com/s/jifnrr3kjq95w64/annotations.tar.gz"

rule install_slncky:
    output: dir("tools/slncky")
    log: "logs/install_slncky.log"
    shell: "git clone https://github.com/slncky/slncky.git tools/ &> {log}"

rule download_annotations:
    input: dir("tools/slncky")
    output: dir("tools/slncky/annotations")
    log: "logs/download_annotations.log"
    shell: "(cd tools/slncky && wget -c {SLNCKY_ANNOTATIONS} -O - | tar -xz) &> {log}"

rule index_fa:
    input: GENOME_FA
    output: GENOME_FA + ".fai"
    shell: "samtools faidx {input}"

rule convert_predicted_gtf_to_bed12:
    '''Use tools from BEDOPS to convert custom GTF to sorted BED12 file'''
    input:
        "output/{fq}" + unique_tag() + STAR_OUT_BAM_SUFFIX_2ND + ".filtered.gtf",
        fai=GENOME_FA + ".fai"
    output:
        temp("output/{fq}" + unique_tag() + STAR_OUT_BAM_SUFFIX_2ND + ".filtered.gtf.genePred"),
        bed12=temp("output/{fq}" + unique_tag() + STAR_OUT_BAM_SUFFIX_2ND + ".filtered.bed12"),
        bed12sorted="output/{fq}" + unique_tag() + STAR_OUT_BAM_SUFFIX_2ND + ".filtered.sorted.bed12"
    log: "output/{fq}" + unique_tag() + STAR_OUT_BAM_SUFFIX_2ND + ".filtered.bed12.log"
    shell:
        "(gtfToGenePred {input} {input}.genePred && "
        "genePredToBed {input}.genePred {output.bed12} && "
        "bedtools sort -faidx {input.fai} -i {output.bed12sorted}) &> {log}"

rule convert_firstBed12_column_to_ucsc:
    '''The first column of the BED12 needs to be converted from Ensembl format (1, 2, MT, etc.) to UCSC format (chr1, chr2, chrM, etc.)'''
    input: "output/{fq}" + unique_tag() + STAR_OUT_BAM_SUFFIX_2ND + ".filtered.sorted.bed12"
    output: "output/{fq}" + unique_tag() + STAR_OUT_BAM_SUFFIX_2ND + ".filtered.sorted.ucsc.bed12"
    shell: # need to write this

rule annotate_lncrnas:
    input:
        dir("tools/slncky"),
        dir("data/annotations"),
        bed12= "output/{fq}" + unique_tag() + STAR_OUT_BAM_SUFFIX_2ND + ".filtered.sorted.ucsc.bed12"
    output:
        "output/{fq}" + unique_tag() + STAR_OUT_BAM_SUFFIX_2ND + "slncky/annotated.canonical_to_lncs.txt",
        "output/{fq}" + unique_tag() + STAR_OUT_BAM_SUFFIX_2ND + "slncky/annotated.cluster_info.txt",
        "output/{fq}" + unique_tag() + STAR_OUT_BAM_SUFFIX_2ND + "slncky/annotated.filtered_info.txt",
        "output/{fq}" + unique_tag() + STAR_OUT_BAM_SUFFIX_2ND + "slncky/annotated.lncs.bed",
        "output/{fq}" + unique_tag() + STAR_OUT_BAM_SUFFIX_2ND + "slncky/annotated.lncs.info.txt",
        "output/{fq}" + unique_tag() + STAR_OUT_BAM_SUFFIX_2ND + "slncky/annotated.orfs.txt",
        "output/{fq}" + unique_tag() + STAR_OUT_BAM_SUFFIX_2ND + "slncky/annotated.orthologs.top.txt",
        "output/{fq}" + unique_tag() + STAR_OUT_BAM_SUFFIX_2ND + "slncky/annotated.orthologs.txt",
    params: reference="hg38"
    threads: 6
    shell:
        "tools/slncky/slncky.v1.0 --threads {threads} {input.bed12} {params.reference} output/{fq}" + unique_tag() + STAR_OUT_BAM_SUFFIX_2ND + "slncky/annotated"
