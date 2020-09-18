# DEPENDENCIES


Everything hecatomb requires can be installed with conda (or mamba)!

### Notes

- Almost everything here is going to be in [bioconda](https://bioconda.github.io/), and so you will need to add `-c conda-forge -c bioconda` to all installs. Alternatively, you can add those channels to your conda config: `conda --add channels conda-forge; conda --add channels bioconda` before you begin
- We use [mamba](https://anaconda.org/conda-forge/mamba) instead of conda, but you can use conda if you don't mind waiting for the environment resolution!
- Therefore, when you see a command like `mamba install cd-hit` you can also consider it equivalent to `conda install -c conda-forge -c bioconda cd-hit`

# Install everything

This conda command should install all dependencies for you

```
mamba install -c conda-forge -c bioconda cd-hit prinseq-plus-plus bowtie2 bbmap samtools
```

For **download_databases.snakefile** you will need:
- [cURL](https://en.wikipedia.org/wiki/CURL) - This should already be on your system
- [cd-hit](http://weizhongli-lab.org/cd-hit/) - install with conda: `mamba install cd-hit`

For **hecatomb_alt.snakefile**
- [prinseq++](https://github.com/Adrian-Cantu/PRINSEQ-plus-plus) - install with conda: `mamba install prinseq-plus-plus`
- [bbtools](https://jgi.doe.gov/data-and-tools/bbtools/) - install with conda: `mamba install bbmap`
- [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) - install with conda: `mamba install bowtie2`
- [samtools](http://www.htslib.org/) - install with conda: `mamba install samtools`
