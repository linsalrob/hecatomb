"""

The hecatomb snakefile has all the parts of the hecatomb pipeline, however
you will need to download the databases before we can begin!

This snakefile also uses conda. To run this, use


```
snakemake --profile slurm --configfile config.yaml -s $HECATOMB/snakemake/hecatomb_alt.snakefile --use-conda --conda-frontend mamba
```


Rob Edwards, Jan-Sep 2020
"""

import os
import sys

if not config:
    sys.stderr.write("FATAL: Please define a config file using the --configfile command line option.\n")
    sys.stderr.write("examples are provided in the Git repo\n")
    sys.exit()


###################################################################
#                                                                 #
# The output directories where results and analyses are written   #
#                                                                 #
###################################################################

# paths for our data. This is where we will read and put things
READDIR = config['Paths']['Reads']
CLUMPED = config['Output']["Clumped"]
QC = config['Output']['QC']
if not os.path.exists(QC):
    os.makedirs(QC, exist_ok=True)

RESULTS = config['Output']['Results']

###################################################################
#                                                                 #
# We create subdirectories in this temp directory as needed       #
#                                                                 #
###################################################################

TMPDIR = config['Paths']['Temp']
if not os.path.exists(TMPDIR):
    os.makedirs(TMPDIR, exist_ok=True)


###################################################################
#                                                                 #
# The database and path structure. Note that below we define the  #
# databases that we are actaully looking for, and test            #
# to see if each of them exist. You should be able                #
# to automatically download them.                                 #
#                                                                 #
###################################################################

DBDIR = config['Paths']['Databases']
BACBT2 = os.path.join(DBDIR, "bac_giant_unique_species", "bac_uniquespecies_giant.masked_Ns_removed")
HOSTBT2 = os.path.join(DBDIR, "human_masked", "human_virus_masked")
CONPATH = os.path.join(DBDIR, "contaminants")
PROTPATH = os.path.join(DBDIR, "proteins")
NUCLPATH = os.path.join(DBDIR, "nucleotides")
TAXPATH  = os.path.join(DBDIR, "taxonomy")

###################################################################
#                                                                 #
# The bacterial and host bowtie2 indexes                          #
#                                                                 #
###################################################################

for bti in [BACBT2, HOSTBT2]:
    if not os.path.exists(f"{bti}.1.bt2l") and not os.path.exists(f"{bti}.1.bt2"):
        sys.stderr.write(f"FATAL: You do not appear to have the bowtie2 indexes for {bti}\n")
        sys.stderr.write("This version of hecatomb uses bowtie to remove host, bacteria, and contaminants\n")
        sys.stderr.write("Please re-run the alt version of download databases\n")
        sys.exit()

if not os.path.exists(os.path.join(CONPATH, "line_sine.1.bt2")):
    sys.stderr.write("FATAL: You do not appear to have the bowtie2 indexes for the line/sine database\n")
    sys.stderr.write("This version of hecatomb uses bowtie to remove contaminants\n")
    sys.stderr.write("Please re-run the alt version of download databases\n")
    sys.exit()

if not os.path.exists(PROTPATH):
    sys.stderr.write("FATAL: You appear not to have the protein databases. Please download the databases using the download_databases.snakefile\n")
    sys.exit()

###################################################################
#                                                                 #
# Amino acid definitions and paths for searches                   #
#                                                                 #
###################################################################

AA_OUT  = os.path.join(RESULTS, "mmseqs_aa_out")
if not os.path.exists(AA_OUT):
    os.makedirs(AA_OUT, exist_ok=True)

AA_OUT_CHECKED  = os.path.join(RESULTS, "mmseqs_aa_checked_out")
if not os.path.exists(AA_OUT_CHECKED):
    os.makedirs(AA_OUT_CHECKED, exist_ok=True)

###################################################################
#                                                                 #
# nucleotide definitions and paths for searches                   #
#                                                                 #
###################################################################

NTDB = os.path.join(NUCLPATH, "refseq_virus_nt_UniVec_masked", "nt.fnaDB")
if not os.path.exists(NTDB):
    sys.stderr.write(f"FATAL: You appear not to have the nucleotide ")
    sys.stderr.write(f"database {NTDB} installed.\n")
    sys.stderr.write(f"Please download the databases using the download_databases.snakefile\n")
    sys.exit()


NT_OUT = os.path.join(RESULTS, "mmseqs_nt_out")
if not os.path.exists(NT_OUT):
    os.makedirs(NT_OUT)

NT_CHECKED_OUT = os.path.join(RESULTS, "mmseqs_nt_checked_out")
if not os.path.exists(NT_CHECKED_OUT):
    os.makedirs(NT_CHECKED_OUT)

# note that we have two databases called "nt.fnaDB". Sorry.
BVMDB = os.path.join(NUCLPATH, "bac_virus_masked", "nt.fnaDB")
if not os.path.exists(BVMDB):
    sys.stderr.write(f"FATAL: You appear not to have the nucleotide ")
    sys.stderr.write(f"database {BVMDB} installed.\n")
    sys.stderr.write(f"Please download the databases using the download_databases.snakefile\n")
    sys.exit()



###################################################################
#                                                                 #
# Taxonomy databases and related information                      #
#                                                                 #
###################################################################

TAXTAX = os.path.join(TAXPATH, "taxonomizr_accessionTaxa.sql")
if not os.path.exists(TAXTAX):
    sys.stderr.write(f"FATAL: You appear not to have the taxonomizr ")
    sys.stderr.write(f"database {TAXTAX} installed.\n")
    sys.stderr.write(f"Please download the databases using the download_databases.snakefile\n")
    sys.exit()


PHAGE_LINEAGES = os.path.join(DBDIR, "phages", "phage_taxonomic_lineages.txt")
if not os.path.exists(PHAGE_LINEAGES):
    sys.stderr.write("FATAL: phages/phage_taxonomic_lineages.txt not ")
    sys.stderr.write("found in the databases directory. Please check ")
    sys.stderr.write("you have the latest version of the databases\n")
    sys.exit()

###################################################################
#                                                                 #
# Uniprot databases and related information                       #
#                                                                 #
###################################################################

URVPATH = os.path.join(PROTPATH, "uniref_plus_virus")
URVDB = os.path.join(URVPATH, "uniref50_virus.db") # uniref50 + viruses database
if not os.path.exists(URVDB):
    sys.stderr.write("FATAL: {URVDB} not found.\n")
    sys.stderr.write("Please make sure that you have run ")
    sys.stderr.write("download_databases.snakefile before commencing\n")
    sys.exit()

VIRDB = os.path.join(PROTPATH, "uniprot_virus_c99.db")
if not os.path.exists(VIRDB):
    sys.stderr.write(f"FATAL: {VIRDB} does not exist. Please ensure you")
    sys.stderr.write(" have installed the databases\n")
    sys.exit()



# how much memory we have
XMX = config['System']['Memory']


###################################################################
#                                                                 #
# Read the sequence files and parse the file names.               #
#                                                                 #
###################################################################

SAMPLES,EXTENSIONS = glob_wildcards(os.path.join(READDIR, '{sample}_R1{extentions}'))

if not EXTENSIONS:
    sys.stderr.write("""
        FATAL: We could not parse the sequence file names.
        We are expecting {sample}_R1{extension}, and so your files
        should contain the characters '_R1' in the fwd reads
        and '_R2' in the rev reads
        """)
    sys.exit()
# we just get the generic extension. This is changed in Step 1

file_extension = EXTENSIONS[0]
# a convenience so we don't need to use '{sample}_R1' all the time
PATTERN_R1 = '{sample}_R1'
PATTERN_R2 = '{sample}_R2'

if len(SAMPLES) == 0:
    sys.stderr.write("FATAL: We could not detect any samples at all.\n")
    sys.stderr.write("You should complain to Rob\n")
    sys.exit()

rule all:
    input:
        os.path.join(AA_OUT, "phage_tax_table.tsv"),
        os.path.join(AA_OUT, "viruses_tax_table.tsv"),
        os.path.join(AA_OUT, "pviral_aa_unclassified_seqs.fasta"),
        os.path.join(AA_OUT_CHECKED, "taxonomyResult.firsthit.m8"),
        os.path.join(AA_OUT_CHECKED, "taxonomyResult.report"),
        os.path.join(AA_OUT_CHECKED, "viruses_checked_aa_table.tsv"),
        os.path.join(AA_OUT_CHECKED, "viruses_checked_aa_tax_table.tsv"),
        os.path.join(AA_OUT_CHECKED, "unclassified_checked_aa_seqs.fasta"),
        os.path.join(NT_OUT, "resultDB.firsthit.m8"),
        os.path.join(NT_CHECKED_OUT, "mmseqs_pviral_nt_lineage.tsv"),
        os.path.join(NT_CHECKED_OUT, "phage_nt_seqs.fasta"),
        os.path.join(NT_CHECKED_OUT, "pviral_virus_nt_seqs.fasta"),
        os.path.join(NT_CHECKED_OUT, "mmseqs_pviral_nt_checked_lineage.tsv"),
        os.path.join(RESULTS, "viruses_tax_table.tsv"),
        os.path.join(RESULTS, "phage_tax_table.tsv"),
        os.path.join(RESULTS, "aa.aln.m8"),
        os.path.join(RESULTS, "nt.aln.m8"),
        "family_reads"

"""
Clean the data.

This is from contaminant_removal.snakefile

"""

rule trim_low_quality:
    """
    Step 0: Trim low-quality bases.

    Here we use prinseq++ to trim and dereplicate the sequences by quaity
    """


    input:
        r1 = os.path.join(READDIR, PATTERN_R1 + file_extension),
        r2 = os.path.join(READDIR, PATTERN_R2 + file_extension)
    output:
        r1 = os.path.join(QC, "step_0", PATTERN_R1 + ".good_out.s0.fastq"),
        r2 = os.path.join(QC, "step_0", PATTERN_R2 + ".good_out.s0.fastq"),
        s1 = os.path.join(QC, "step_0", PATTERN_R1 + ".single_out.s0.fastq"),
        s2 = os.path.join(QC, "step_0", PATTERN_R2 + ".single_out.s0.fastq"),
        b1 = temporary(os.path.join(QC, "step_0", PATTERN_R1 + ".bad_out_R1.fastq")),
        b2 = temporary(os.path.join(QC, "step_0", PATTERN_R2 + ".bad_out_R2.fastq"))
    benchmark:
        "benchmarks/trim_low_quality_{sample}.txt"
    resources:
        mem_mb=20000,
        cpus=8
    params:
        o = os.path.join(QC, "step_0")
    conda:
        "envs/prinseq.yaml"
    shell:
        """
            prinseq++ -min_len 60 -min_qual_mean 25 -ns_max_n 1 -derep 1 \
                    -out_format 0 -trim_tail_left 5 -trim_tail_right 5 \
                    -ns_max_n 5  -trim_qual_type min -trim_qual_left 30 \
                    -trim_qual_right 30 -trim_qual_window 10 \
                    -threads {resources.cpus} \
                    -out_good {output.r1} -out_single {output.s1} -out_bad {output.b1} \
                    -out_good2 {output.r2} -out_single2 {output.s2} -out_bad2 {output.b2} \
                    -fastq {input.r1} \
                    -fastq2 {input.r2};
        """


rule remove_leftmost_primerB:
    """
    Step 1: Remove leftmost primerB. Not the reverse complements
    """
    input:
        r1 = os.path.join(QC, "step_0", PATTERN_R1 + ".good_out.s0.fastq"),
        r2 = os.path.join(QC, "step_0", PATTERN_R2 + ".good_out.s0.fastq"),
        primers = os.path.join(CONPATH, "primerB.fa")
    output:
        r1 = os.path.join(QC, "step_1", PATTERN_R1 + ".s1.out.fastq"),
        r2 = os.path.join(QC, "step_1", PATTERN_R2 + ".s1.out.fastq"),
        stats = os.path.join(QC, "step_1", "{sample}.s1.stats.txt")
    benchmark:
        "benchmarks/removeprimerB_{sample}.txt"
    resources:
        mem_mb=20000,
        cpus=8
    conda:
        "envs/bbmap.yaml"
    shell:
        """
        bbduk.sh in={input.r1} in2={input.r2} \
            ref={input.primers} \
            out={output.r1} out2={output.r2} \
            stats={output.stats} \
            k=16 hdist=1 mink=11 ktrim=l restrictleft=20 \
            removeifeitherbad=f trimpolya=10 ordered=t \
            rcomp=f ow=t {XMX}
        """

rule remove_3prime_contaminant:
    """
    Step 2: Remove 3' read through contaminant
    """
    input:
        r1 = os.path.join(QC, "step_1", PATTERN_R1 + ".s1.out.fastq"),
        r2 = os.path.join(QC, "step_1", PATTERN_R2 + ".s1.out.fastq"),
        primers = os.path.join(CONPATH, "rc_primerB_ad6.fa")
    output:
        r1 = os.path.join(QC, "step_2", PATTERN_R1 + ".s2.out.fastq"),
        r2 = os.path.join(QC, "step_2", PATTERN_R2 + ".s2.out.fastq"),
        stats = os.path.join(QC, "step_2", "{sample}.s2.stats.txt")
    benchmark:
        "benchmarks/remove_3prime_contaminant_{sample}.txt"
    resources:
        mem_mb=20000,
        cpus=8
    conda:
        "envs/bbmap.yaml"
    shell:
        """
        bbduk.sh in={input.r1} in2={input.r2} \
            ref={input.primers} \
            out={output.r1} out2={output.r2} \
            stats={output.stats} \
            k=16 hdist=1 mink=11 ktrim=r removeifeitherbad=f \
            ordered=t rcomp=f ow=t {XMX}
        """

rule remove_primer_free_adapter:
    """
    Step 3: Remove primer free adapter (both orientations)
    """
    input:
        r1 = os.path.join(QC, "step_2", PATTERN_R1 + ".s2.out.fastq"),
        r2 = os.path.join(QC, "step_2", PATTERN_R2 + ".s2.out.fastq"),
        primers = os.path.join(CONPATH, "nebnext_adapters.fa")
    output:
        r1 = os.path.join(QC, "step_3", PATTERN_R1 + ".s3.out.fastq"),
        r2 = os.path.join(QC, "step_3", PATTERN_R2 + ".s3.out.fastq"),
        stats = os.path.join(QC, "step_3", "{sample}.s3.stats.txt")
    benchmark:
        "benchmarks/remove_primer_free_adapter_{sample}.txt"
    resources:
        mem_mb=20000,
        cpus=8
    conda:
        "envs/bbmap.yaml"
    shell:
        """
        bbduk.sh in={input.r1} in2={input.r2} \
            ref={input.primers} \
            out={output.r1} out2={output.r2} \
            stats={output.stats} \
            k=16 hdist=1 mink=10 ktrim=r removeifeitherbad=f \
            ordered=t rcomp=t ow=t {XMX}
        """

rule remove_adapter_free_primer:
    """
    Step 4: Remove adapter free primer (both orientations)
    """
    input:
        r1 = os.path.join(QC, "step_3", PATTERN_R1 + ".s3.out.fastq"),
        r2 = os.path.join(QC, "step_3", PATTERN_R2 + ".s3.out.fastq"),
        primers = os.path.join(CONPATH, "rc_primerB_ad6.fa")
    output:
        r1 = os.path.join(QC, "step_4", PATTERN_R1 + ".s4.out.fastq"),
        r2 = os.path.join(QC, "step_4", PATTERN_R2 + ".s4.out.fastq"),
        stats = os.path.join(QC, "step_4", "{sample}.s4.stats.txt")
    benchmark:
        "benchmarks/remove_adapter_free_primer_{sample}.txt"
    resources:
        mem_mb=20000,
        cpus=8
    conda:
        "envs/bbmap.yaml"
    shell:
        """
        bbduk.sh in={input.r1} in2={input.r2} \
            ref={input.primers} \
            out={output.r1} out2={output.r2} \
            stats={output.stats} \
            k=16 hdist=0 removeifeitherbad=f ordered=t rcomp=t ow=t {XMX}
        """

rule remove_vector_contamination:
    """
    Step 5: Vector contamination removal (PhiX + NCBI UniVecDB)
    """
    input:
        r1 = os.path.join(QC, "step_4", PATTERN_R1 + ".s4.out.fastq"),
        r2 = os.path.join(QC, "step_4", PATTERN_R2 + ".s4.out.fastq"),
        primers = os.path.join(CONPATH, config['DatabaseFiles']['contaminants'])
    output:
        r1 = os.path.join(QC, "step_5", PATTERN_R1 + ".s5.out.fastq"),
        r2 = os.path.join(QC, "step_5", PATTERN_R2 + ".s5.out.fastq"),
        stats = os.path.join(QC, "step_5", "{sample}.s5.stats.txt")
    benchmark:
        "benchmarks/remove_vector_contamination_{sample}.txt"
    resources:
        mem_mb=20000,
        cpus=8
    conda:
        "envs/bbmap.yaml"
    shell:
        """
        bbduk.sh in={input.r1} in2={input.r2} \
            ref={input.primers} \
            out={output.r1} out2={output.r2} \
            stats={output.stats} \
            k=31 hammingdistance=1 ordered=t ow=t {XMX}
        """

rule bowtie2_host_removal:
    """
    Step 6. Remove host and LINE/SINES
    These steps have been converted to bowtie2 from bbmap as it is more efficient
    """
    input:
        r1 = os.path.join(QC, "step_5", PATTERN_R1 + ".s5.out.fastq"),
        r2 = os.path.join(QC, "step_5", PATTERN_R2 + ".s5.out.fastq"),
    output:
        os.path.join(QC, "step_6", '{sample}.host.bam')
    benchmark:
        "benchmarks/{rules.myrule.name}_{sample}.txt"
    resources:
        mem_mb=20000,
        cpus=8
    params:
        idx = HOSTBT2
    conda:
        "envs/bowtie2.yaml"
    shell:
        """
        bowtie2 -p {resources.cpus} -x {params.idx} -1 {input.r1} -2 {input.r2} | \
        samtools view --threads {resources.cpus} -bh | \
        samtools sort --threads {resources.cpus} -o {output} -
        """

rule reads_host_mapped:
    """
    Note, in deconseq.snakefile these are all separate rules,
    but I merge them here so that they run on the same cluster node
    and it will reduce the size of the dag
    """
    input:
        os.path.join(QC, "step_6", '{sample}.host.bam')
    output:
        r1 = os.path.join(QC, "step_6", PATTERN_R1 + '_host.mapped.fastq'),
        r2 = os.path.join(QC, "step_6", PATTERN_R2 + '_host.mapped.fastq'),
        s = os.path.join(QC, "step_6", '{sample}_singletons_host.mapped.fastq')
    benchmark:
        "benchmarks/{rules.myrule.name}_{sample}.txt"
    resources:
        mem_mb=20000,
        cpus=8
    conda:
        "envs/bowtie2.yaml"
    shell:
        """
        samtools fastq --threads {resources.cpus} -G 12 -f 65 {input} > {output.r1} &&
        samtools fastq --threads {resources.cpus} -G 12 -f 129 {input} > {output.r2} &&
        samtools fastq --threads {resources.cpus} -F 5 {input} > {output.s}
        """

rule reads_host_unmapped:
    input:
        os.path.join(QC, "step_6", '{sample}.host.bam')
    output:
        r1 = os.path.join(QC, "step_6", PATTERN_R1 + '_host.unmapped.fastq'),
        r2 = os.path.join(QC, "step_6", PATTERN_R2 + '_host.unmapped.fastq'),
        s = os.path.join(QC, "step_6", '{sample}_singletons_host.unmapped.fastq')
    benchmark:
        "benchmarks/{rules.myrule.name}_{sample}.txt"
    resources:
        mem_mb=20000,
        cpus=8
    conda:
        "envs/samtools.yaml"
    shell:
        """
        samtools fastq --threads {resources.cpus} -f 77  {input} > {output.r1} && 
        samtools fastq --threads {resources.cpus} -f 141 {input} > {output.r2} &&
        samtools fastq --threads {resources.cpus} -f 4 -F 1  {input} > {output.s}
        """

rule line_sine_bam:
    input:
        r1 = os.path.join(QC, "step_6", PATTERN_R1 + '_host.unmapped.fastq'),
        r2 = os.path.join(QC, "step_6", PATTERN_R2 + '_host.unmapped.fastq'),
        s = os.path.join(QC, "step_6", '{sample}_singletons_host.unmapped.fastq')
    output:
        os.path.join(QC, "step_6", '{sample}.linesine.bam')
    benchmark:
        "benchmarks/{rules.myrule.name}_{sample}.txt"
    resources:
        mem_mb=20000,
        cpus=8
    params:
        idx = os.path.join(CONPATH, "line_sine")
    conda:
        "envs/bowtie2.yaml"
    shell:
        """
        bowtie2 -p {resources.cpus} -x {params.idx} -1 {input.r1} -2 {input.r2} \
        -U {input.s} | samtools view --threads {resources.cpus} -bh | \
        samtools sort --threads {resources.cpus} -o {output} -
        """

rule reads_linesine_mapped:
    input:
        os.path.join(QC, "step_6", '{sample}.linesine.bam')
    output:
        r1 = os.path.join(QC, "step_6", PATTERN_R1 + '_linesine.mapped.fastq'),
        r2 = os.path.join(QC, "step_6", PATTERN_R2 + '_linesine.mapped.fastq'),
        s = os.path.join(QC, "step_6", '{sample}_singletons_linesine.mapped.fastq')
    benchmark:
        "benchmarks/{rules.myrule.name}_{sample}.txt"
    resources:
        mem_mb=20000,
        cpus=8
    conda:
        "envs/samtools.yaml"
    shell:
        """
        samtools fastq --threads {resources.cpus} -G 12 -f 65 {input} > {output.r1} &&
        samtools fastq --threads {resources.cpus} -G 12 -f 129 {input} > {output.r2} &&
        samtools fastq --threads {resources.cpus} -F 5 {input} > {output.s}
        """

rule reads_linesine_unmapped:
    input:
        os.path.join(QC, "step_6", '{sample}.linesine.bam')
    output:
        r1 = os.path.join(QC, "step_6", PATTERN_R1 + ".s6.out.fastq"),
        r2 = os.path.join(QC, "step_6", PATTERN_R2 + ".s6.out.fastq"),
        s = os.path.join(QC, "step_6", '{sample}_singletons.s6.out.fastq')
    benchmark:
        "benchmarks/{rules.myrule.name}_{sample}.txt"
    resources:
        mem_mb=20000,
        cpus=8
    conda:
        "envs/samtools.yaml"
    shell:
        """
        samtools fastq --threads {resources.cpus}  -f 77  {input} > {output.r1} && 
        samtools fastq --threads {resources.cpus}  -f 141 {input} > {output.r2} &&
        samtools fastq --threads {resources.cpus}  -f 4 -F 1  {input} > {output.s}
        """

rule remove_bacteria:
    """
    Step 8: Remove bacterial contaminants reserving viral and ambiguous sequences
    """
    input:
        r1 = os.path.join(QC, "step_6", PATTERN_R1 + ".s6.out.fastq"),
        r2 = os.path.join(QC, "step_6", PATTERN_R2 + ".s6.out.fastq"),
        s = os.path.join(QC, "step_6", '{sample}_singletons.s6.out.fastq')
    output:
        os.path.join(QC, "step_7", '{sample}.bacteria.bam')
    params:
        idx = BACBT2
    benchmark:
        "benchmarks/{rules.myrule.name}_{sample}.txt"
    resources:
        mem_mb=20000,
        cpus=8
    conda:
        "envs/bowtie2.yaml"
    shell:
        """
        bowtie2 -p {resources.cpus} -x {params.idx} -1 {input.r1} -2 {input.r2} \
        -U {input.s} | samtools view --threads {resources.cpus} -bh | \
        samtools sort --threads {resources.cpus} -o {output} -
        """

rule reads_bacteria_mapped:
    input:
        os.path.join(QC, "step_7", '{sample}.bacteria.bam')
    output:
        r1 = os.path.join(QC, "step_7", PATTERN_R1 + '_bacteria.mapped.fastq'),
        r2 = os.path.join(QC, "step_7", PATTERN_R2 + '_bacteria.mapped.fastq'),
        s = os.path.join(QC, "step_7", '{sample}_singletons_bacteria.mapped.fastq')
    benchmark:
        "benchmarks/{rules.myrule.name}_{sample}.txt"
    resources:
        mem_mb=20000,
        cpus=8
    conda:
        "envs/samtools.yaml"
    shell:
        """
        samtools fastq --threads {resources.cpus} -G 12 -f 65 {input} > {output.r1} &&
        samtools fastq --threads {resources.cpus}  -G 12 -f 129 {input} > {output.r2} &&
        samtools fastq  --threads {resources.cpus} -F 5 {input} > {output.s}
        """

rule reads_bacteria_unmapped:
    input:
        os.path.join(QC, "step_7", '{sample}.bacteria.bam')
    output:
        r1 = os.path.join(QC, "step_7", PATTERN_R1 + ".s7.out.fastq"),
        r2 = os.path.join(QC, "step_7", PATTERN_R2 + ".s7.out.fastq"),
        s = os.path.join(QC, "step_7", '{sample}_singletons.s7.out.fastq')
    benchmark:
        "benchmarks/{rules.myrule.name}_{sample}.txt"
    resources:
        mem_mb=20000,
        cpus=8
    conda:
        "envs/samtools.yaml"
    shell:
        """
        samtools fastq --threads {resources.cpus} -f 77  {input} > {output.r1} && 
        samtools fastq --threads {resources.cpus} -f 141 {input} > {output.r2} &&
        samtools fastq --threads {resources.cpus} -f 4 -F 1  {input} > {output.s}
        """

"""
End the contaminant_removal.snakefile



This section is from cluster_count.snakefile and counts all the clusters!
"""



"""
Currently from this point forwards we just work with the R1 reads.
In a previous iteration we merged the reads but that has been removed
"""

"""

There were two steps, `Remove exact duplicates` and `Dereplicate` that
have been removed, because prinseq++ does those steps, approximately. 

The dedupliate step allows for additional mismatches, but we are not 
doing that here and sticking more with the OTU approach.

"""

rule deduplicate:
    """
    Step 10: Dereplicate
    """
    input:
        r1 = os.path.join(QC, "step_7", PATTERN_R1 + ".s7.out.fastq"),
    output:
        fa = os.path.join(QC, "step_10", "{sample}_R1.best.fasta"),
        stats = os.path.join(QC, "step_10", "{sample}_R1.stats.txt")
    benchmark:
        "benchmarks/{rules.myrule.name}_{sample}.txt"
    resources:
        mem_mb=20000,
        cpus=8
    conda:
        "envs/bbmap.yaml"
    shell:
        """
        dedupe.sh in={input.r1} \
            csf={output.stats} out={output.fa} \
            ow=t s=4 rnc=t pbr=t {XMX}
        """


rule extract_seq_counts:
    """
    Step 11: Extract sequences and counts for seqtable (count table)
    """
    input:
        os.path.join(QC, "step_10", "{sample}_R1.best.fasta")
    output:
        os.path.join(QC, "step_11", "{sample}_R1.reformated.fasta")
    benchmark:
        "benchmarks/{rules.myrule.name}_{sample}.txt"
    resources:
        mem_mb=20000,
        cpus=8
    conda:
        "envs/bbmap.yaml"
    shell:
        """
        reformat.sh in={input} out={output} \
            deleteinput=f fastawrap=0 \
            ow=t \
            {XMX}
        """

rule extract_counts:
    """
    Step 12. Parse and combine stats and contig files
    """
    input:
        os.path.join(QC, "step_11", "{sample}_R1.reformated.fasta") 
    output:
        os.path.join(QC, "counts", "{sample}_R1.seqs.txt")
    shell:
        """
        grep -v '>' {input} | sed '1i sequence' > {output}
        """

rule extract_counts_ids:
    """
    Step 12b. Extract sequence IDs
    """
    input:
        os.path.join(QC, "step_11", "{sample}_R1.reformated.fasta")
    output:
        os.path.join(QC, "counts", "{sample}_R1.contig_ids.txt")
    shell:
        """
        grep '>' {input} | sed 's|>Cluster_||' | awk -F "," '{{ print$1 }}' | sort -n | sed '1i contig_ids' > {output}
        """

rule exract_count_stats:
    """
    Step 12c. Extract counts
    """
    input:
        os.path.join(QC, "step_10", "{sample}_R1.stats.txt")
    output:
        os.path.join(QC, "counts", "{sample}_R1.counts.txt")
    params:
        s = "{sample}"
    shell:
        """
        cut -f 2 {input} | sed "1s/size/{params.s}/" > {output}
        """

rule create_seq_table:
    """
    Step 12d. Create sequence table
    """
    input:
        seq = os.path.join(QC, "counts", "{sample}_R1.seqs.txt"),
        cnt = os.path.join(QC, "counts", "{sample}_R1.counts.txt")
    output:
        os.path.join(QC, "counts", "{sample}_seqtable.txt")
    shell:
        """
        paste {input.seq} {input.cnt} > {output}
        """

rule merge_seq_table:
    """
    Merge seq counts
    """
    input:
        files = expand(os.path.join(QC, "counts", "{sample}_seqtable.txt"), sample=SAMPLES)
    output:
        seqtable = os.path.join(RESULTS, "seqtable_all.tsv"),
        tab2fx = os.path.join(RESULTS, "seqtable.tab2fx")
    benchmark:
        "benchmarks/{rules.myrule.name}.txt"
    resources:
        mem_mb=20000,
        cpus=8
    params:
        resultsdir = directory(RESULTS),
    conda:
        "envs/R.yaml"
    script:
        "scripts/seqtable_merge.R"




"""
End cluster_count.snakefile

This section is from  mmseqs_pviral_aa.snakefile

"""


rule convert_seqtable_to_fasta:
    input:
        os.path.join(RESULTS, "seqtable.tab2fx")
    output:
        os.path.join(RESULTS, "seqtable.fasta")
    shell:
        "sed -e 's/^/>/; s/\\t/\\n/' {input} > {output}"

rule create_seqtable_db:
    input:
        os.path.join(RESULTS, "seqtable.fasta")
    output:
        os.path.join(AA_OUT, "seqtable_query.db")
    benchmark:
        "benchmarks/{rules.myrule.name}.txt"
    resources:
        mem_mb=20000,
        cpus=16
    conda:
        "envs/mmseqs2.yaml"
    shell:
        "mmseqs createdb --shuffle 0 --dbtype 0 {input} {output}"

rule seqtable_taxsearch:
    input:
        sq = os.path.join(AA_OUT, "seqtable_query.db"),
    output:
        tr = os.path.join(AA_OUT, "taxonomyResult.dbtype")
    params:
        tr = os.path.join(AA_OUT, "taxonomyResult")
    benchmark:
        "benchmarks/{rules.myrule.name}.txt"
    resources:
        mem_mb=20000,
        cpus=16
    conda:
        "envs/mmseqs2.yaml"
    shell:
        """
        mmseqs taxonomy {input.sq} {VIRDB} {params.tr} $(mktemp -d -p {TMPDIR}) \
        -a --start-sens 1 --sens-steps 3 -s 7  --threads {resources.cpus} \
        --search-type 2 --tax-output-mode 1
        """

rule seqtable_convert_alignments:
    input:
        sq = os.path.join(AA_OUT, "seqtable_query.db"),
        tr = os.path.join(AA_OUT, "taxonomyResult.dbtype")
    params:
        tr = os.path.join(AA_OUT, "taxonomyResult")
    benchmark:
        "benchmarks/{rules.myrule.name}.txt"
    resources:
        mem_mb=20000,
        cpus=16
    output:
        os.path.join(AA_OUT, "aln.m8")
    conda:
        "envs/mmseqs2.yaml"
    shell:
        """
        mmseqs convertalis {input.sq} {VIRDB} {params.tr} {output} --threads {resources.cpus} \
        --format-output "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qaln,taln"
        """

rule seqtable_lca:
    input:
        tr = os.path.join(AA_OUT, "taxonomyResult.dbtype")
    output:
        os.path.join(AA_OUT, "lca.db.dbtype")
    params:
        lc = os.path.join(AA_OUT, "lca.db"),
        tr = os.path.join(AA_OUT, "taxonomyResult")
    benchmark:
        "benchmarks/{rules.myrule.name}.txt"
    resources:
        mem_mb=20000,
        cpus=16
    conda:
        "envs/mmseqs2.yaml"
    shell:
        """
        mmseqs lca {VIRDB} {params.tr} {params.lc} --tax-lineage 1 --threads {resources.cpus} \
        --lca-ranks "superkingdom,phylum,class,order,family,genus,species";
        """

rule seqtable_taxtable_tsv:
    input:
        sq = os.path.join(AA_OUT, "seqtable_query.db"),
        lc = os.path.join(AA_OUT, "lca.db.dbtype")
    params:
        lc = os.path.join(AA_OUT, "lca.db")
    output:
        os.path.join(AA_OUT, "taxonomyResult.tsv")
    benchmark:
        "benchmarks/{rules.myrule.name}.txt"
    resources:
        mem_mb=20000,
        cpus=16
    conda:
        "envs/mmseqs2.yaml"
    shell:
        """
        mmseqs createtsv --threads {resources.cpus} {input.sq} {params.lc} {output}
        """

rule seqtable_create_kraken:
    input:
        lc = os.path.join(AA_OUT, "lca.db")
    output:
        os.path.join(AA_OUT, "taxonomyResult.report")
    benchmark:
        "benchmarks/{rules.myrule.name}.txt"
    resources:
        mem_mb=20000,
        cpus=16
    conda:
        "envs/mmseqs2.yaml"
    shell:
        """
        mmseqs taxonomyreport --threads {resources.cpus} {VIRDB} {input.lc} {output}
        """

## Adjust taxonomy table and extract viral lineages
# Extract all (virus + phage) potential viral sequences

rule find_viruses:
    input:
        os.path.join(AA_OUT, "taxonomyResult.tsv")
    output:
        os.path.join(AA_OUT, "all_viruses_table.tsv")
    shell:
        """
        grep 'Viruses;' {input} | cut -f1,5 | sed 's/phi14:2/phi14_2/g' | \
                sed 's/;/\\t/g' | \
                sort -n -k1 > {output}
        """

# Extract phage viral lineages and generate taxonomy table for import into R as PhyloSeq object
rule find_phages:
    input:
        av = os.path.join(AA_OUT, "all_viruses_table.tsv")
    output:
        os.path.join(AA_OUT, "phage_table.tsv")
    shell:
        "grep -f {PHAGE_LINEAGES} {input.av} > {output}"

rule find_phage_seqs:
    input:
        os.path.join(AA_OUT, "phage_table.tsv")
    output:
        os.path.join(AA_OUT, "phage_seqs.list")
    shell:
        "cut -f1 {input} > {output}"

rule pull_phage_seqs:
    input:
        fa = os.path.join(RESULTS, "seqtable.fasta"),
        ls = os.path.join(AA_OUT, "phage_seqs.list")
    output:
        os.path.join(AA_OUT, "phage_seqs.fasta")
    shell:
        """
        grep --no-group-separator -A 1 -Fwf {input.ls} {input.fa} > {output}
        """

rule phage_seqs_to_tab:
    input:
        os.path.join(AA_OUT, "phage_seqs.fasta")
    output:
        os.path.join(AA_OUT, "phage_seqs.tab")
    shell:
        """
        perl -pe 'if (s/^>//) {{chomp; s/$/\t/}}' {input} > {output}
        """

rule phage_to_tax_table:
    input:
        tab = os.path.join(AA_OUT, "phage_seqs.tab"),
        tsv = os.path.join(AA_OUT, "phage_table.tsv")
    output:
        os.path.join(AA_OUT, "phage_tax_table.tsv")
    shell:
        """
        join {input.tab} {input.tsv} | \
        cut -d ' ' --output-delimiter=$'\t' -f 2-9 | \
        sed '1isequence\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies' \
        > {output}
        """

# Extract non-phage viral lineages and generate taxonomy table for import into R as PhyloSeq object
rule find_non_phages:
    input:
        av = os.path.join(AA_OUT, "all_viruses_table.tsv")
    output:
        os.path.join(AA_OUT, "viruses_table.tsv")
    shell:
        "grep -vf {PHAGE_LINEAGES} {input.av} > {output}"

rule find_non_phage_seqs:
    input:
        os.path.join(AA_OUT, "viruses_table.tsv")
    output:
        os.path.join(AA_OUT, "viruses_seqs.list")
    shell:
        "cut -f1 {input} > {output}"

rule pull_non_phage_seqs:
    input:
        fa = os.path.join(RESULTS, "seqtable.fasta"),
        ls = os.path.join(AA_OUT, "viruses_seqs.list")
    output:
        os.path.join(AA_OUT, "viruses_seqs.fasta")
    shell:
        """
        grep --no-group-separator -A 1 -Fwf {input.ls} {input.fa} > {output}
        """

rule non_phage_seqs_to_tab:
    input:
        os.path.join(AA_OUT, "viruses_seqs.fasta")
    output:
        os.path.join(AA_OUT, "viruses_seqs.tab")
    shell:
        """
        perl -pe 'if (s/^>//) {{chomp; s/$/\t/}}' {input} > {output}
        """

rule non_phage_to_tax_table:
    input:
        tab = os.path.join(AA_OUT, "viruses_seqs.tab"),
        tsv = os.path.join(AA_OUT, "viruses_table.tsv")
    output:
        os.path.join(AA_OUT, "viruses_tax_table.tsv")
    shell:
        """
        join {input.tab} {input.tsv} | \
        cut -d ' ' --output-delimiter=$'\t' -f 2-9 | \
        sed '1isequence\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies' \
        > {output}
        """


# Extract unclassified lineages
rule unclassified_lineages:
    input:
        os.path.join(AA_OUT, "taxonomyResult.tsv")
    output:
        os.path.join(AA_OUT, "pviral_unclassified_seqs.list")
    shell:
        """
        grep -v 'Viruses;' {input} | cut -f1 | \
                sort -n -k1 > {output}
        """

rule pull_unclassified_seqs:
    input:
        fa = os.path.join(RESULTS, "seqtable.fasta"),
        ls = os.path.join(AA_OUT, "pviral_unclassified_seqs.list")
    output:
        os.path.join(AA_OUT, "pviral_aa_unclassified_seqs.fasta")
    shell:
        """
        grep --no-group-separator -A 1 -Fwf {input.ls} {input.fa} > {output}
        """



"""
This section is from mmseqs_pviral_aa_check.snakefile

Based upon ../base/mmseqs_pviral_aa_check.sh

Rob Edwards, April 2020

"""



rule create_viral_seqs_db:
    input:
        os.path.join(AA_OUT, "viruses_seqs.fasta")
    output:
        os.path.join(AA_OUT_CHECKED, "viral_seqs_queryDB")
    benchmark:
        "benchmarks/{rules.myrule.name}.txt"
    resources:
        mem_mb=20000,
        cpus=16
    conda:
        "envs/mmseqs2.yaml"
    shell:
        """
        mmseqs createdb {input} {output} --shuffle 0 --dbtype 0
        """

rule viral_seqs_tax_search:
    """
    The taxnomy result maybe several disjoint files
    """
    input:
        vqdb = os.path.join(AA_OUT_CHECKED, "viral_seqs_queryDB"),
    output:
        os.path.join(AA_OUT_CHECKED, "taxonomyResult.dbtype"),
        os.path.join(AA_OUT_CHECKED, "taxonomyResult.index")
    params:
        tr = os.path.join(AA_OUT_CHECKED, "taxonomyResult")
    benchmark:
        "benchmarks/{rules.myrule.name}.txt"
    resources:
        mem_mb=100000,
        cpus=16
    conda:
        "envs/mmseqs2.yaml"
    shell:
        """
        mmseqs taxonomy {input.vqdb} {URVDB} {params.tr} \
            $(mktemp -d -p {TMPDIR}) --threads {resources.cpus} \
            -a -s 7 --search-type 2 --tax-output-mode 1
        """

rule viral_seqs_convertalis:
    input:
        vqdb = os.path.join(AA_OUT_CHECKED, "viral_seqs_queryDB"),
        tr = os.path.join(AA_OUT_CHECKED, "taxonomyResult.dbtype")
    params:
        tr = os.path.join(AA_OUT_CHECKED, "taxonomyResult")
    output:
        os.path.join(AA_OUT_CHECKED, "aln.m8")
    benchmark:
        "benchmarks/{rules.myrule.name}.txt"
    resources:
        mem_mb=20000,
        cpus=16
    conda:
        "envs/mmseqs2.yaml"
    shell:
        """
        mmseqs convertalis {input.vqdb} {URVDB} {params.tr} {output} --threads {resources.cpus} \
        --format-output "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qaln,taln"
        """

rule viral_seqs_lca:
    input:
        vqdb = os.path.join(AA_OUT_CHECKED, "viral_seqs_queryDB"),
        tr = os.path.join(AA_OUT_CHECKED, "taxonomyResult.dbtype")
    params:
        tr = os.path.join(AA_OUT_CHECKED, "taxonomyResult"),
        lca = os.path.join(AA_OUT_CHECKED, "lca.db")
    output:
        os.path.join(AA_OUT_CHECKED, "lca.db.dbtype"),
        os.path.join(AA_OUT_CHECKED, "lca.db.index")
    benchmark:
        "benchmarks/{rules.myrule.name}.txt"
    resources:
        mem_mb=20000,
        cpus=16
    conda:
        "envs/mmseqs2.yaml"
    shell:
        """
        mmseqs lca {URVDB} {params.tr} {params.lca} \
        --tax-lineage 1 --threads {resources.cpus} \
        --lca-ranks "superkingdom,phylum,class,order,family,genus,species"
        """

rule extract_top_hit:
    input:
        os.path.join(AA_OUT_CHECKED, "taxonomyResult.dbtype")
    params:
        tr = os.path.join(AA_OUT_CHECKED, "taxonomyResult"),
        fh = os.path.join(AA_OUT_CHECKED, "taxonomyResult.firsthit")
    output:
        os.path.join(AA_OUT_CHECKED, "taxonomyResult.firsthit.dbtype"),
        os.path.join(AA_OUT_CHECKED, "taxonomyResult.firsthit.index")
    benchmark:
        "benchmarks/{rules.myrule.name}.txt"
    resources:
        mem_mb=20000,
        cpus=16
    conda:
        "envs/mmseqs2.yaml"
    shell:
        """
        mmseqs filterdb {params.tr} {params.fh} --extract-lines 1 --threads {resources.cpus}
        """

rule convertalis_vsqd:
    input:
        vqdb = os.path.join(AA_OUT_CHECKED, "viral_seqs_queryDB"),
        trfhdb = os.path.join(AA_OUT_CHECKED, "taxonomyResult.firsthit.dbtype")
    output:
        os.path.join(AA_OUT_CHECKED, "taxonomyResult.firsthit.m8")
    params:
        trfh = os.path.join(AA_OUT_CHECKED, "taxonomyResult.firsthit")
    benchmark:
        "benchmarks/{rules.myrule.name}.txt"
    resources:
        mem_mb=20000,
        cpus=16
    conda:
        "envs/mmseqs2.yaml"
    shell:
        """
        mmseqs convertalis --threads {resources.cpus} {input.vqdb} {URVDB} {params.trfh} {output} 
        """

rule create_taxtable_vsqd:
    input:
        vqdb = os.path.join(AA_OUT_CHECKED, "viral_seqs_queryDB"),
        lcadb = os.path.join(AA_OUT_CHECKED, "lca.db.dbtype")
    params:
        lcadb = os.path.join(AA_OUT_CHECKED, "lca.db")
    output:
        os.path.join(AA_OUT_CHECKED, "taxonomyResult.tsv")
    benchmark:
        "benchmarks/{rules.myrule.name}.txt"
    resources:
        mem_mb=20000,
        cpus=16
    conda:
        "envs/mmseqs2.yaml"
    shell:
        """
        mmseqs createtsv --threads {resources.cpus} {input.vqdb} {params.lcadb} {output}
        """

rule create_kraken_vsqd:
    input:
        lcadb = os.path.join(AA_OUT_CHECKED, "lca.db.dbtype")
    params:
        lcadb = os.path.join(AA_OUT_CHECKED, "lca.db")
    output:
        os.path.join(AA_OUT_CHECKED, "taxonomyResult.report")
    benchmark:
        "benchmarks/{rules.myrule.name}.txt"
    resources:
        mem_mb=20000,
        cpus=16
    conda:
        "envs/mmseqs2.yaml"
    shell:
        """
        mmseqs taxonomyreport --threads {resources.cpus} {URVDB} {params.lcadb} {output}
        """

rule nonphages_to_pyloseq_table:
    input:
        os.path.join(AA_OUT_CHECKED, "taxonomyResult.tsv")
    output:
        os.path.join(AA_OUT_CHECKED, "viruses_checked_aa_table.tsv")
    shell:
        """
        grep -v 'Bacteria;' {input} | \
            grep 'Viruses;' | \
            grep -v -f  {PHAGE_LINEAGES} | cut -f1,5 | \
            sed 's/;/\t/g' | \
                sort -n -k1 > {output}
        """

rule nonphages_to_pyloseq_list:
    input:
        os.path.join(AA_OUT_CHECKED, "viruses_checked_aa_table.tsv")
    output:
        os.path.join(AA_OUT_CHECKED, "viruses_checked_aa_seqs.list")
    shell:
        """
        cut -f1 {input} > {output}
        """

rule pull_nonphage_seqs:
    input:
        fa = os.path.join(RESULTS, "seqtable.fasta"),
        ls = os.path.join(AA_OUT_CHECKED, "viruses_checked_aa_seqs.list")
    output:
        os.path.join(AA_OUT_CHECKED, "viruses_checked_aa_seqs.fasta")
    shell:
        """
        grep --no-group-separator -A 1 -Fwf {input.ls} {input.fa} > {output}
        """

rule nonphage_seqs_to_tab:
    input:
        os.path.join(AA_OUT_CHECKED, "viruses_checked_aa_seqs.fasta")
    output:
        os.path.join(AA_OUT_CHECKED, "viruses_checked_aa_seqs.tab")
    shell:
        """
        perl -pe 'if (s/^>//) {{chomp; s/$/\t/}}' {input} > {output}
        """

rule nonphage_to_tax_table:
    input:
        tab = os.path.join(AA_OUT_CHECKED, "viruses_checked_aa_seqs.tab"),
        tsv = os.path.join(AA_OUT_CHECKED, "viruses_checked_aa_table.tsv")
    output:
        os.path.join(AA_OUT_CHECKED, "viruses_checked_aa_tax_table.tsv")
    shell:
        """
        join {input.tab} {input.tsv} | \
        cut -d ' ' --output-delimiter=$'\t' -f 2-9 | \
        sed '1isequence\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies' \
        > {output}
        """

rule non_viral_lineages:
    input:
        os.path.join(AA_OUT_CHECKED, "taxonomyResult.tsv")
    output:
        os.path.join(AA_OUT_CHECKED, "unclassified_checked_aa_seqs.list")
    shell:
        """
        grep -v 'Viruses;' {input} | cut -f1 | \
           sort -n -k1 > {output}
        """
 
rule pull_non_viral_lineages:
    input:
        fa = os.path.join(RESULTS, "seqtable.fasta"),
        ls = os.path.join(AA_OUT_CHECKED, "unclassified_checked_aa_seqs.list")
    output:
        os.path.join(AA_OUT_CHECKED, "unclassified_checked_aa_seqs.fasta")
    shell:
        """
        grep --no-group-separator -A 1 -Fwf {input.ls} {input.fa} > {output}
        """



"""

End the code from mmseqs_pviral_aa_check.snakefile

This code comes from mmseqs_pviral_nt.snakefile and is 
based on mmseqs_pviral_nt.sh

"""


rule create_nt_querydb:
    input:
        os.path.join(AA_OUT, "pviral_aa_unclassified_seqs.fasta")
    output:
        idx = os.path.join(NT_OUT, "seqtable_queryDB.index"),
        dbt = os.path.join(NT_OUT, "seqtable_queryDB.dbtype")
    params:
        st = os.path.join(NT_OUT, "seqtable_queryDB")
    benchmark:
        "benchmarks/{rules.myrule.name}.txt"
    resources:
        mem_mb=20000,
        cpus=16
    conda:
        "envs/mmseqs2.yaml"
    shell:
        """
        mmseqs createdb {input} {params.st}  --dbtype 2
        """

rule nt_search:
    input:
        idx = os.path.join(NT_OUT, "seqtable_queryDB.index"),
        dbt = os.path.join(NT_OUT, "seqtable_queryDB.dbtype")
    output:
        idx = os.path.join(NT_OUT, "resultDB.index"),
        dbt = os.path.join(NT_OUT, "resultDB.dbtype")
    params:
        st = os.path.join(NT_OUT, "seqtable_queryDB"),
        rdb = os.path.join(NT_OUT, "resultDB")
    benchmark:
        "benchmarks/{rules.myrule.name}.txt"
    resources:
        mem_mb=20000,
        cpus=16
    conda:
        "envs/mmseqs2.yaml"
    shell:
        """
        mmseqs search {params.st} {NTDB} {params.rdb} $(mktemp -d -p {TMPDIR}) \
        -a -e 0.000001 --search-type 3 --cov-mode 2 -c 0.95 --threads {resources.cpus}
        """

rule nt_top_hit:
    input:
        idx = os.path.join(NT_OUT, "resultDB.index"),
        dbt = os.path.join(NT_OUT, "resultDB.dbtype")
    output:
        os.path.join(NT_OUT, "resultDB.firsthit.dbtype"),
        os.path.join(NT_OUT, "resultDB.firsthit.index")
    params:
        rdb = os.path.join(NT_OUT, "resultDB"),
        rdbfh = os.path.join(NT_OUT, "resultDB.firsthit")
    benchmark:
        "benchmarks/{rules.myrule.name}.txt"
    resources:
        mem_mb=20000,
        cpus=16
    conda:
        "envs/mmseqs2.yaml"
    shell:
        """
        mmseqs filterdb {params.rdb} {params.rdbfh} --extract-lines 1 --threads {resources.cpus}
        """

rule nt_to_m8:
    input:
        sidx = os.path.join(NT_OUT, "seqtable_queryDB.index"),
        sdbt = os.path.join(NT_OUT, "seqtable_queryDB.dbtype"),
        ridx = os.path.join(NT_OUT, "resultDB.firsthit.dbtype"),
        rdbt = os.path.join(NT_OUT, "resultDB.firsthit.index")
    output:
        os.path.join(NT_OUT, "resultDB.firsthit.m8")
    params:
        st = os.path.join(NT_OUT, "seqtable_queryDB"),
        rdbfh = os.path.join(NT_OUT, "resultDB.firsthit")
    benchmark:
        "benchmarks/{rules.myrule.name}.txt"
    resources:
        mem_mb=20000,
        cpus=16
    conda:
        "envs/mmseqs2.yaml"
    shell:
        """
        mmseqs convertalis {params.st} {NTDB} {params.rdbfh} {output} --threads {resources.cpus}
        """

rule nt_annotate:
    input:
        fhtbl = os.path.join(NT_OUT, "resultDB.firsthit.m8")
    output:
        linout = os.path.join(NT_CHECKED_OUT, "mmseqs_pviral_nt_lineage.tsv")
    params:
        taxtax = TAXTAX
    conda:
        "envs/R.yaml"
    script:
        "scripts/mmseqs_pviral_nt_annotate.R"



"""
End mmseqs_pviral_nt.snakefile

Based on mmseqs_pviral_nt_check.sh

"""


rule find_nt_phages:
    input:
        # this input comes from mmseqs_pviral_nt.snakefile
        os.path.join(NT_CHECKED_OUT, "mmseqs_pviral_nt_lineage.tsv")
    output:
        os.path.join(NT_CHECKED_OUT, "phage_nt_table.tsv")
    shell:
        """
        tail -n+2 {input} | grep -f {PHAGE_LINEAGES} |  sort -n -k1 \
                > {output}
        """

rule list_nt_phages:
    input:
        os.path.join(NT_CHECKED_OUT, "phage_nt_table.tsv")
    output:
         os.path.join(NT_CHECKED_OUT, "phage_nt_table.list")
    shell:
        """
        cut -f1 {input} > {output}
        """

rule pull_nt_phage_seqs:
    input:
        fa = os.path.join(RESULTS, "seqtable.fasta"),
        ls = os.path.join(NT_CHECKED_OUT, "phage_nt_table.list")
    output:
        os.path.join(NT_CHECKED_OUT, "phage_nt_seqs.fasta")
    shell:
        """
        grep --no-group-separator -A 1 -Fwf {input.ls} {input.fa} > {output}
        """

rule non_nt_phage_viruses:
    input:
        os.path.join(NT_CHECKED_OUT, "mmseqs_pviral_nt_lineage.tsv")
    output:
        os.path.join(NT_CHECKED_OUT, "pviral_virus_nt_table.tsv")
    shell:
        """
        tail -n+2 {input} | grep -vf {PHAGE_LINEAGES} |  sort -n -k1 \
                > {output}
        """

rule list_non_nt_viruses:
    input:
        os.path.join(NT_CHECKED_OUT, "pviral_virus_nt_table.tsv")
    output:
        os.path.join(NT_CHECKED_OUT, "pviral_virus_nt_table.list")
    shell:
        """
        cut -f1 {input} > {output}
        """

rule pull_nt_non_phage_seqs:
    input:
        fa = os.path.join(RESULTS, "seqtable.fasta"),
        ls = os.path.join(NT_CHECKED_OUT, "pviral_virus_nt_table.list")
    output:
        os.path.join(NT_CHECKED_OUT, "pviral_virus_nt_seqs.fasta")
    shell:
        """
        grep --no-group-separator -A 1 -Fwf {input.ls} {input.fa} > {output}
        """

rule create_nt_db:
    input:
        os.path.join(NT_CHECKED_OUT, "pviral_virus_nt_seqs.fasta")
    output:
         idx = os.path.join(NT_CHECKED_OUT, "seqtable_queryDB.index"),
         dbt = os.path.join(NT_CHECKED_OUT, "seqtable_queryDB.dbtype")
    params:
         st = os.path.join(NT_CHECKED_OUT, "seqtable_queryDB")
    benchmark:
        "benchmarks/{rules.myrule.name}.txt"
    resources:
        mem_mb=20000,
        cpus=16
    conda:
        "envs/mmseqs2.yaml"
    shell:
        """
        mmseqs createdb {input} {params.st} --dbtype 2 
        """

rule nt_search_checked:
    input:
        idx = os.path.join(NT_CHECKED_OUT, "seqtable_queryDB.index"),
        dbt = os.path.join(NT_CHECKED_OUT, "seqtable_queryDB.dbtype")
    output:
        idx = os.path.join(NT_CHECKED_OUT, "resultDB.index"),
        dbt = os.path.join(NT_CHECKED_OUT, "resultDB.dbtype")
    params:
        st = os.path.join(NT_CHECKED_OUT, "seqtable_queryDB"),
        rdb = os.path.join(NT_CHECKED_OUT, "resultDB")
    benchmark:
        "benchmarks/{rules.myrule.name}.txt"
    resources:
        mem_mb=20000,
        cpus=16
    conda:
        "envs/mmseqs2.yaml"
    shell:
        """
        mmseqs search {params.st} {NTDB} {params.rdb} $(mktemp -d -p {TMPDIR}) \
        -a -e 0.000001 --search-type 3 --cov-mode 2 -c 0.95 --threads {resources.cpus}
        """


rule filter_nt_db:
    input:
        idx = os.path.join(NT_CHECKED_OUT, "resultDB.index"),
        dbt = os.path.join(NT_CHECKED_OUT, "resultDB.dbtype")
    output:
        idx = os.path.join(NT_CHECKED_OUT, "resultDB.firsthit.index"),
        dbt = os.path.join(NT_CHECKED_OUT, "resultDB.firsthit.dbtype")
    params:
        rdb = os.path.join(NT_CHECKED_OUT, "resultDB"),
        rfh = os.path.join(NT_CHECKED_OUT, "resultDB.firsthit")
    benchmark:
        "benchmarks/{rules.myrule.name}.txt"
    resources:
        mem_mb=20000,
        cpus=16
    conda:
        "envs/mmseqs2.yaml"
    shell:
        """
        mmseqs filterdb {params.rdb} {params.rfh}  --extract-lines 1 --threads {resources.cpus}
        """

rule convert_nt_alias:
    input:
        sqi = os.path.join(NT_CHECKED_OUT, "seqtable_queryDB.index"),
        sqd = os.path.join(NT_CHECKED_OUT, "seqtable_queryDB.dbtype"),
        idx = os.path.join(NT_CHECKED_OUT, "resultDB.firsthit.index"),
        dbt = os.path.join(NT_CHECKED_OUT, "resultDB.firsthit.dbtype")
    output:
        os.path.join(NT_CHECKED_OUT, "resultDB.firsthit.m8")
    params:
        fh = os.path.join(NT_CHECKED_OUT, "resultDB.firsthit"),
        sq = os.path.join(NT_CHECKED_OUT, "seqtable_queryDB"),
    benchmark:
        "benchmarks/{rules.myrule.name}.txt"
    resources:
        mem_mb=20000,
        cpus=16
    conda:
        "envs/mmseqs2.yaml"
    shell:
        """
        mmseqs convertalis {params.sq} {BVMDB} {params.fh} {output} --threads {resources.cpus}
        """

rule annotate_checked_nt:
    input:
        fhtbl = os.path.join(NT_CHECKED_OUT, "resultDB.firsthit.m8")
    output:
        linout = os.path.join(NT_CHECKED_OUT, "mmseqs_pviral_nt_checked_lineage.tsv")
    params:
        taxtax = TAXTAX
    conda:
        "envs/R.yaml"
    script:
        "scripts/mmseqs_pviral_nt_check_annotate.R"


"""

End mmseqs_pviral_nt_check.snakefile
this section comes from concatenate.snakefile
based on concatenate_results.sh

take all the results and make some cool data sets from them!
"""


rule jive_aa_annotation:
    input:
        os.path.join(AA_OUT_CHECKED, "viruses_checked_aa_table.tsv")
    output:
        os.path.join(AA_OUT_CHECKED, "viruses_checked_aa_table_edited.tsv")
    shell:
        """
        sed 's/uc_//g' {input} > {output}
        """

rule jive_nt_annotation:
    input:
        os.path.join(NT_CHECKED_OUT, "mmseqs_pviral_nt_checked_lineage.tsv")
    output:
        os.path.join(NT_CHECKED_OUT, "mmseqs_pviral_nt_checked_lineage_edited.tsv")
    shell:
        """
        tail -n +2 {input} | sed 's/NA/unknown/g; s/uc_//g' > {output}
        """

# this is the name of the rule in concatenate_results.sh and it is too cute to change
rule happily_marry_outputs:
    input:
        aa = os.path.join(AA_OUT_CHECKED, "viruses_checked_aa_table_edited.tsv"),
        nt = os.path.join(NT_CHECKED_OUT, "mmseqs_pviral_nt_checked_lineage_edited.tsv")
    output:
        temporary(os.path.join(RESULTS, "viruses_tax_table_tmp.tsv"))
    shell:
        """
        cat {input.aa} {input.nt} | sort -n -k 1 > {output}
        """

rule add_crown_to_marriage:
    # not really, its a title
    input:
        os.path.join(RESULTS, "viruses_tax_table_tmp.tsv")
    output:
        os.path.join(RESULTS, "viruses_tax_table.tsv")
    shell:
        """
        sed -e '1iid\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies' {input} > {output}
        """

rule fix_phage_names:
    input:
        tsv = os.path.join(AA_OUT, "phage_table.tsv")
    output:
        temporary(os.path.join(RESULTS, "phage_tax_table_temp.tsv"))
    shell:
        """
        sed 's/uc_//g' {input} > {output}
        """

rule fix_phage_cols:
    input:
        os.path.join(RESULTS, "phage_tax_table_temp.tsv")
    output:
        temporary(os.path.join(RESULTS, "phage_tax_table_temp2.tsv"))
    shell:
        """
        cut -f1-8 {input} > {output}
        """
            
rule add_phage_title:
    input:
        os.path.join(RESULTS, "phage_tax_table_temp2.tsv")
    output:
        os.path.join(RESULTS, "phage_tax_table.tsv")
    shell:
        """
        sed -e '1iid\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies' {input} > {output}
        """

rule add_aa_tax_header:
    input:
        os.path.join(AA_OUT_CHECKED, "taxonomyResult.firsthit.m8")
    output:
        os.path.join(RESULTS, "aa.aln.m8")
    shell:
        """
        sed -e '1iquery\ttarget\tpercent_id\talignment_length\tnum_mismatches\tnumber_gaps\tstart_query\tend_query\tstart_target\tend_target\te_value\tbit_score' \
                {input} > {output}
        """

rule add_nt_tax_header:
    input:
        os.path.join(NT_OUT, "resultDB.firsthit.m8")
    output:
        os.path.join(RESULTS, "nt.aln.m8")
    shell:
        """
        sed -e '1iquery\ttarget\tpercent_id\talignment_length\tnum_mismatches\tnumber_gaps\tstart_query\tend_query\tstart_target\tend_target\te_value\tbit_score' \
                {input} > {output}
        """

"""
I am not very happy with this solution and looking for some alternative ideas
Basically this is very linear, and I can't figure out how to make it non-linear!

I think the solution is to make a list of all virus familes (e.g. read
that from a file) and then use expand on that list to create each of the
fasta files and then finally remove any empty files, but this will work for now
"""

rule make_fasta:
    input:
        os.path.join(RESULTS, "viruses_tax_table.tsv")
    output:
        directory("family_reads")
    shell:
        """
        mkdir -p {output} && 
        for FAM in $(tail -n +2 {input} | cut -f6 | awk '!s[$0]++');
        do 
            grep $FAM {input} | cut -f 1 | \
                grep -A1 -Fwf - results/seqtable.fasta > {output}/$FAM.fasta
        done
        """







