# MetaLoRe
A Modular and Reproducible Pipeline for Long-Read Metagenomic Analysis Using ONT and PacBio Sequencing Data


## Easy Installation (Linux) (recommended)
Download apptainer image (XXX.gb) (apptainer is similar to docker but safe for clusters)

```
apptainer pull --arch amd64 library://wheaton5/souporcell/souporcell:release

```

If you are running on a scientific cluster, they will likely have apptainer, contact your sysadmin for more details. If you are running on your own linux box you may need to install  [apptainer](https://apptainer.org/docs/user/main/quick_start.html)
However, if you find the installation process difficult, you can use conda to install it.
```
conda install -y apptainer
```
requires singularity >= 3.0 or apptainer version >= 1.0.0 (singularity rebranded as apptainer and changed its version numbers)

```
which apptainer
apptainer --version

```
## Usage

You should now be able to run XXXX.py and XXXX.py through the singularity container. Singularity automatically mounts the current working directory and directories downstream from where you run it, otherwise you would need to manually mount those directories. 

### Single-Sample Analysis Overview 
MetaLoRe’s single-sample module automates the full processing of long-read metagenomic data (ONT or PacBio).This module is run using MetaLoRe-Base.py. Starting from raw FASTQ, it performs quality filtering, optional host-read removal, and generates interactive QC reports. The pipeline then assembles contigs (Flye), converts to FASTA, annotates genes (Prokka/Abricate), classifies taxonomy (Centrifuge/Krona), and bins genomes (MetaBAT2, MaxBin2, SemiBin2), producing standardized output files and publication-ready plots for each sample.
The options for using MetaLoRe-Base.py are:

```
singularity exec souporcell_latest.sif souporcell_pipeline.py -h
usage: souporcell_pipeline.py [-h] -i BAM -b BARCODES -f FASTA -t THREADS -o
                              OUT_DIR -k CLUSTERS [-p PLOIDY]
                              [--min_alt MIN_ALT] [--min_ref MIN_REF]
                              [--max_loci MAX_LOCI] [--restarts RESTARTS]
                              [--common_variants COMMON_VARIANTS]
                              [--known_genotypes KNOWN_GENOTYPES]
                              [--known_genotypes_sample_names KNOWN_GENOTYPES_SAMPLE_NAMES [KNOWN_GENOTYPES_SAMPLE_NAMES ...]]
                              [--skip_remap SKIP_REMAP] [--ignore IGNORE]

single cell RNAseq mixed genotype clustering using sparse mixture model
clustering with tensorflow.

optional arguments:
  -h, --help            show this help message and exit
  -i BAM, --bam BAM     cellranger bam
  -b BARCODES, --barcodes BARCODES
                        barcodes.tsv from cellranger
  -f FASTA, --fasta FASTA
                        reference fasta file
  -t THREADS, --threads THREADS
                        max threads to use
  -o OUT_DIR, --out_dir OUT_DIR
                        name of directory to place souporcell files
  -k CLUSTERS, --clusters CLUSTERS
                        number cluster, tbd add easy way to run on a range of
                        k
  -p PLOIDY, --ploidy PLOIDY
                        ploidy, must be 1 or 2, default = 2
  --min_alt MIN_ALT     min alt to use locus, default = 10.
  --min_ref MIN_REF     min ref to use locus, default = 10.
  --max_loci MAX_LOCI   max loci per cell, affects speed, default = 2048.
  --restarts RESTARTS   number of restarts in clustering, when there are > 12
                        clusters we recommend increasing this to avoid local
                        minima
                         --common_variants COMMON_VARIANTS
                        common variant loci or known variant loci vcf, must be
                        vs same reference fasta
  --known_genotypes KNOWN_GENOTYPES
                        known variants per clone in population vcf mode, must
                        be .vcf right now we dont accept gzip or bcf sorry
  --known_genotypes_sample_names KNOWN_GENOTYPES_SAMPLE_NAMES [KNOWN_GENOTYPES_SAMPLE_NAMES ...]
                        which samples in population vcf from known genotypes
                        option represent the donors in your sample
  --skip_remap SKIP_REMAP
                        don't remap with minimap2 (not recommended unless in
                        conjunction with --common_variants
  --ignore IGNORE       set to True to ignore data error assertions

```

A typical command looks like

```
singularity exec /path/to/souporcell_latest.sif souporcell_pipeline.py -i /path/to/possorted_genome_bam.bam -b /path/to/barcodes.tsv -f /path/to/reference.fasta -t num_threads_to_use -o output_dir_name -k num_clusters

```

### Differential Analysis Overview 
Once all samples have been processed, MetaLoRe’s comparative module merges taxonomic and functional abundance tables and conducts group-wise statistics and diversity analyses. It supports LEfSe and STAMP for differential feature discovery, computes α-diversity (Shannon, Simpson) and β-diversity (Bray–Curtis) metrics, and performs ordination (PCoA/NMDS) with PERMANOVA. The result is a comprehensive set of tables and figures highlighting biomarkers, group differences, and overall community structure across sample groups.
The options for using  MetaLoRe-Compare.py are:

```
singularity exec souporcell_latest.sif souporcell_pipeline.py -h
usage: souporcell_pipeline.py [-h] -i BAM -b BARCODES -f FASTA -t THREADS -o
                              OUT_DIR -k CLUSTERS [-p PLOIDY]
                              [--min_alt MIN_ALT] [--min_ref MIN_REF]
                              [--max_loci MAX_LOCI] [--restarts RESTARTS]
                              [--common_variants COMMON_VARIANTS]
                              [--known_genotypes KNOWN_GENOTYPES]
                              [--known_genotypes_sample_names KNOWN_GENOTYPES_SAMPLE_NAMES [KNOWN_GENOTYPES_SAMPLE_NAMES ...]]
                              [--skip_remap SKIP_REMAP] [--ignore IGNORE]

single cell RNAseq mixed genotype clustering using sparse mixture model
clustering with tensorflow.

optional arguments:
  -h, --help            show this help message and exit
  -i BAM, --bam BAM     cellranger bam
  -b BARCODES, --barcodes BARCODES
                        barcodes.tsv from cellranger
  -f FASTA, --fasta FASTA
                        reference fasta file
  -t THREADS, --threads THREADS
                        max threads to use
  -o OUT_DIR, --out_dir OUT_DIR
                        name of directory to place souporcell files
  -k CLUSTERS, --clusters CLUSTERS
                        number cluster, tbd add easy way to run on a range of
                        k
  -p PLOIDY, --ploidy PLOIDY
                        ploidy, must be 1 or 2, default = 2
  --min_alt MIN_ALT     min alt to use locus, default = 10.
  --min_ref MIN_REF     min ref to use locus, default = 10.
  --max_loci MAX_LOCI   max loci per cell, affects speed, default = 2048.
  --restarts RESTARTS   number of restarts in clustering, when there are > 12
                        clusters we recommend increasing this to avoid local
                        minima
                         --common_variants COMMON_VARIANTS
                        common variant loci or known variant loci vcf, must be
                        vs same reference fasta
  --known_genotypes KNOWN_GENOTYPES
                        known variants per clone in population vcf mode, must
                        be .vcf right now we dont accept gzip or bcf sorry
  --known_genotypes_sample_names KNOWN_GENOTYPES_SAMPLE_NAMES [KNOWN_GENOTYPES_SAMPLE_NAMES ...]
                        which samples in population vcf from known genotypes
                        option represent the donors in your sample
  --skip_remap SKIP_REMAP
                        don't remap with minimap2 (not recommended unless in
                        conjunction with --common_variants
  --ignore IGNORE       set to True to ignore data error assertions

```

A typical command looks like

```
singularity exec /path/to/souporcell_latest.sif souporcell_pipeline.py -i /path/to/possorted_genome_bam.bam -b /path/to/barcodes.tsv -f /path/to/reference.fasta -t num_threads_to_use -o output_dir_name -k num_clusters

```
