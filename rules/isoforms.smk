rule assemble_transcripts:
    input:
        bam="output/{fq}" + unique_tag() + STAR_OUT_BAM_SUFFIX_2ND,
        gff=GFF3 + ".fix.gff3"
    output: "output/{fq}.gtf"
    threads: 6
    log: "output/{fq}.gtf.log"
    shell:
        "stringtie {input.bam} -p {threads} -G {input.gff} -o {output} 2> {log}" # strandedness: --fr for forwared or --rf for reverse

rule merge_transcripts:
    input:
        custom_gtfs=expand("output/{fq}.gtf", fq=FQ_PREFIXES),
        gff=GFF3 + ".fix.gff3"
    output: "output/combined.gtf"
    threads: 12
    log: "output/combined.gtf.log"
    shell:
        "stringtie --merge -o {output} -G {input.gff} -p {threads} -i {input.custom_gtfs} 2> {log}"

rule build_gtf_sharp:
    output: "GtfSharp/GtfSharp/bin/Release/netcoreapp2.1/GtfSharp.dll"
    log: "output/GtfSharp.build.log"
    shell:
        "(cd GtfSharp && "
        "dotnet restore && "
        "dotnet build -c Release GtfSharp.sln) &> {log}"

rule filter_transcripts_add_cds:
    input:
        gtfsharp="GtfSharp/GtfSharp/bin/Release/netcoreapp2.1/GtfSharp.dll",
        gtf="output/combined.gtf",
        fa=FA,
        efa="data/ERCC.fa",
        refg=GFF3 + ".fix.gff3"
    output:
        temp("output/combined.filtered.gtf"),
        "output/combined.filtered.withcds.gtf",
        efa=temp("ensembl/efa.fa")
    log: "output/filter_add_cds.log"
    shell:
        "cat {input.fa} <(grep -v '>ERCC_ID' {input.efa}) > {output.efa} &&"
        "dotnet {input.gtfsharp} -f {output.efa} -g {input.gtf} -r {input.refg} &> {log}"
