
# Installing hecatomb from source

hecatomb has a few requirements:

- [R](https://www.r-project.org/)
    - [tidyverse](https://www.tidyverse.org/packages/)
- [snakemake](https://snakemake.readthedocs.io/)
- [bbtools](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/installation-guide/)


Please note that `tidyverse` requires libxml2-dev or libxml2-devel depending on your operating system.

To install the required `R` packages, you can run the accessory [install_packages.R](accessory/install_packages.R) like so: `Rscript install_packages.R`. Note that you do not need to be root to install these packages.

# Conda installation

Please check the [dependencies](DEPENDENCIES.md) for instructions on how to install all the dependencies using CONDA.

# The config file

You should make a copy of the config file. We typically make a copy of that file into each directory where we are working. Then if you make any changes to that file they reside with the data. 
There are [example config files](configs/) in both [JSON](configs/sample_config.json) and [YAML](configs/sample_config.yaml), and of course snakemake can use either. (If you are not sure, YAML is probably easier to start with than JSON).

The key things in the config file are:

1. The database file location. You can set that in the config file and then create the database as described below
2. The directory name where your raw reads (`fastq files`) reside. 

You can adjust almost everything else as needed, and the directories will be created when the data is generated.

# Resources

The snakemake conda applications require use of the `resources` options available in snakemake to request additional resources (memory/time/cpus) on an as-needed basis (i.e. not every application requests these by default). We strongly recommend setting up default configurations for these.

Create a directory `~/.config/snakemake/slurm/` and then create the file `~/.config/snakemake/slurm/config.yaml` with these contents:


```
jobs: 10
cluster: "sbatch -t {resources.time_min} --mem={resources.mem_mb} -c {resources.cpus} -o logs_slurm/{rule}_{jobid}.out -e logs_slurm/{rule}_{jobid}.err "
default-resources: [cpus=1, mem_mb=2000, time_min=60]
latency-wait: 60
local-cores: 32
```

Now run snakemake as `snakemake --profile slurm --configfile config.yaml -s snakefile` to use those default settings.


# Setting up the databases

Before you begin, you need to set up the databases. We have several different databases that we screen the data against:

- bacterial genomes
- primer and vector contamination
- host (we typically screen against human, but you can substitute or append to this).

You can easily download and compile the databases as described in the [databases/](databases/) directory. This will take a few minutes but you will only need to do it once.

*Note:* The database download is 1.6 GB, and the uncompressed databases require 32 GB of disk space after extraction and compilation.

# Testing hecatomb

Once you have the databases installed you can run hecatomb on the test data that we have provided.

```bash
cd test_data
snakemake --snakefile ~/GitHubs/hecatomb/snakemake/hecatomb.snakefile --configfile config.yaml
```


