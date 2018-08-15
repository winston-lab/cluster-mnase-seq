
# *spt6-1004* intragenic TSS MNase-seq clustering

## description

Code used to cluster *spt6-1004* intragenic TSSs based on the surrounding MNase-seq signal, as described in [our preprint](https://www.biorxiv.org/content/early/2018/06/15/347575). To reproduce the clustering exactly as in the publication, use the archived version of this analysis from [Zenodo](https://doi.org/10.5281/zenodo.1325930).

## requirements

### required software

- Unix-like operating system (tested on CentOS 7.2.1511)
- Git
- [conda](https://conda.io/docs/user-guide/install/index.html)
- [snakemake](https://snakemake.readthedocs.io/en/stable/)

### required files

Again, to reproduce the clustering from the publication, use the archived version of this analysis in the `cluster_intragenic_mnase` directory from [Zenodo](https://doi.org/10.5281/zenodo.1325930).

- annotation file of the summit positions of *spt6-1004* upregulated intragenic TSSs, generated by the [TSS-seq pipeline](https://github.com/winston-lab/tss-seq) in the `tss-seq-publication` directory of the Zenodo archive.
- matrix of wild-type and *spt6-1004* MNase-seq coverage flanking *spt6-1004* upregulated intragenic TSSs, generated by the [MNase-seq pipeline](https://github.com/winston-lab/mnase-seq) in the `mnase-seq-publication` directory of the Zenodo archive.

## instructions

**0**.  If running the analysis in the Zenodo archive, the pipelines in the `tss-seq-publication` and `mnase-seq-publication` directories need to be run before running this pipeline.

**1**. Run the pipeline using [snakemake](https://snakemake.readthedocs.io/en/stable/).

```bash
snakemake -p --use-conda
```


