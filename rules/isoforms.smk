rule assemble_transcripts:
    input:
        bam="data/combined.sorted.bam",
        gff=GFF3 + ".fix.gff3"
    output: "data/combined.sorted.gtf"
    threads: 12
    log: "data/combined.sorted.gtf.log"
    shell:
        "stringtie {input.bam} -p {threads} -G {input.gff} -o {output} 2> {log}" # strandedness: --fr for forwared or --rf for reverse

rule build_gtf_sharp:
    output: "GtfSharp/GtfSharp/bin/Release/netcoreapp2.1/GtfSharp.dll"
    log: "data/GtfSharp.build.log"
    shell:
        "(cd GtfSharp && "
        "dotnet restore && "
        "dotnet build -c Release GtfSharp.sln) &> {log}"

rule filter_transcripts_add_cds:
    input:
        gtfsharp="GtfSharp/GtfSharp/bin/Release/netcoreapp2.1/GtfSharp.dll",
        gtf="data/combined.sorted.gtf",
        fa=FA,
        refg=GFF3 + ".fix.gff3"
    output:
        "data/combined.sorted.filtered.gtf",
        "data/combined.sorted.filtered.withcds.gtf",
    shell:
        "mv {input.gtf} temporary/combined.sorted.gtf && "
        "dotnet {input.gtfsharp} -f {input.fa} -g temporary/combined.sorted.gtf -r {input.refg} && "
        "mv temporary/* data"
