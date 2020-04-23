# Snakemake analysis pipeline for lncRNA discovery and quantification

Use conda to create and activate the conda environment:
1. Download and install a conda distributions, such as `miniconda`: https://docs.conda.io/en/latest/miniconda.html
2. From this directory run `conda env create -n lncRNAs -f environment.yaml; conda activate lncRNAs`. This will create an environment with all the tools necessary to run this pipeline.
3. Recommended command to run this pipeline: `snakemake --cores 24 --resources mem_mb=100000`

If you use this pipeline, please cite:
* Cesnik, A. J.; Yang, B.; Truong, A.; Spiniello, M.; Steinbrink, M.; Shortreed, M. R.; Frey, B. L.; Jarrard, D. F.; Smith, L. M. “Long Noncoding RNAs AC009014.3 and Newly Discovered XPLAID Differentiate Aggressive and Indolent Prostate Cancers.” Translational Oncology, 2018, 11, 808–814.

This pipeline uses the following tools:
  * `skewer`: Jiang H, Lei R, Ding S-W, and Zhu S (2014). Skewer: a fast and accurate adapter trimmer for next-generation sequencing paired-end reads. BMC Bioinformatics 15(1), 182.
  * `RSEM`: Li B and Dewey CN (2011). RSEM: accurate transcript quantification from RNA-Seq data with or without a reference genome. BMC Bioinformatics 12 (1), 323.
  * `STAR`: Dobin A, Davis CA, Schlesinger F, Drenkow J, Zaleski C, Jha S, Batut P, Chaisson M, and Gingeras TR (2013). STAR: ultrafast universal RNA-seq aligner. Bioinformatics 29(1), 15–21.
  * `stringtie`: Pertea M, Pertea GM, Antonescu CM, Chang TC, Mendell JT & Salzberg SL. StringTie enables improved reconstruction of a transcriptome from RNA-seq reads Nature Biotechnology 2015, doi:10.1038/nbt.3122.
  * `slncky`: Chen J, Shishkin AA, Zhu X, Kadri S, Maza I, Guttman M, Hanna JH, Regev A, and Garber M (2016). Evolutionary analysis across mammals reveals distinct classes of long non-coding RNAs. Genome Biol 17(1), 19.
