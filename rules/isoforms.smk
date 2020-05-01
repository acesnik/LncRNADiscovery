rule assemble_transcripts:
    input:
        bam="output/combined.sorted.bam",
        gff=GFF3 + ".fix.gff3"
    output: "output/combined.sorted.gtf"
    threads: 12
    log: "output/combined.sorted.gtf.log"
    shell:
        "stringtie {input.bam} -p {threads} -G {input.gff} -o {output} 2> {log}" # strandedness: --fr for forwared or --rf for reverse

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
        gtf="output/combined.sorted.gtf",
        fa=FA,
        refg=GFF3 + ".fix.gff3"
    output:
        temp("output/combined.sorted.filtered.gtf"),
        "output/combined.sorted.filtered.withcds.gtf",
    log: "output/filter_add_cds.log"
    shell:
        "dotnet {input.gtfsharp} -f {input.fa} -g {input.gtf} -r {input.refg} &> {log}"
