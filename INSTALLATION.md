
# Hecatomb is now dependent on conda and snakemake alone.

You only need to install `conda` and `snakemake` and `hecatomb` will do the rest for you. No more software installs!

## Step 1. Install `conda`. 

This is undoubtedly the hardest part, and it is quite easy. Install [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) following their instructions.

## Step 2. Install `snakemake` and `mamba`

We highly recommend you install `mamba` as it is a _lot_ faster than conda. But you don't have to (if you don't, skip the mamba step below)

```
conda install -c conda-forge -c bioconda snakemake mamba
```


### Step 3. Clone this repo

We will move this into a lightweight repo soon. At the moment you have all the unnecessary development cruft here

```
git clone https://github.com/linsalrob/hecatomb.git
```


### Step 4. Set up your snakemake environment

We recommend using profiles (see these [great](https://www.sichong.site/2020/02/25/snakemake-and-slurm-how-to-manage-workflow-with-resource-constraint-on-hpc/) and [great](http://bluegenes.github.io/Using-Snakemake_Profiles/) blogs for more information). 

Our default profile encompasses both snakemake and slurm details. If you are not using slurm, then feel free to leave those parts out.

You will need to make a directory in you `.config` directory and add this information:

```bash
mkdir ~/.config/snakemake/slurm/
vi ~/.config/snakemake/slurm/config.yaml
```

Put this in that file:

```yaml
# non-slurm settings

jobs: 10
use-conda: True
conda-frontend: mamba
default-resources: [cpus=1, mem_mb=2000, time_min=60]
keep-going: True


# slurm settings

cluster: "sbatch -t {resources.time_min} --mem={resources.mem_mb} -c {resources.cpus} -o logs_slurm/{rule}_{jobid}.out -e logs_slurm/{rule}_{jobid}.err "
latency-wait: 60
local-cores: 32
```




### Step 5. Run the code

With $HECATOMB as the path to this directory you can run:


```bash
snakemake --profile slurm --configfile config.yaml -s $HECATOMB/snakemake/hecatomb_alt.snakefile 
```

# The config file

You should make a copy of the config file. We typically make a copy of that file into each directory where we are working. Then if you make any changes to that file they reside with the data. 
There are [example config files](configs/) in both [JSON](configs/sample_config.json) and [YAML](configs/sample_config.yaml), and of course snakemake can use either. (If you are not sure, YAML is probably easier to start with than JSON).

The key things in the config file are:

1. The database file location. You can set that in the config file and then create the database as described below
2. The directory name where your raw reads (`fastq files`) reside. 

You can adjust almost everything else as needed, and the directories will be created when the data is generated.


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
snakemake --snakefile $HECATOMB/snakemake/hecatomb.snakefile --configfile config.yaml
```


