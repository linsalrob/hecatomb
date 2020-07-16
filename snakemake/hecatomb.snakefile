"""

The hecatomb snakefile has all the parts of the hecatomb pipeline, however
you will need to download the databases before we can begin!

Rob Edwards, Jan 2020
"""

import os
import sys

if not config:
    sys.stderr.write("FATAL: Please define a config file using the --configfile command line option.\n")
    sys.stderr.write("examples are provided in the Git repo\n")
    sys.exit()


DBDIR = config['Paths']['Databases']
TMPDIR = config['Paths']['Temp']
if not os.path.exists(TMPDIR):
    os.makedirs(TMPDIR, exist_ok=True)


# paths for our databases
BACPATH = os.path.join(DBDIR, "bac_giant_unique_species")
HOSTPATH = os.path.join(DBDIR, "human_masked")
CONPATH = os.path.join(DBDIR, "contaminants")
PROTPATH = os.path.join(DBDIR, "proteins")

if not os.path.exists(os.path.join(HOSTPATH, "ref")):
    sys.stderr.write("FATAL: You appear not to have the host databases. Please download the databases using the download_databases.snakefile\n")
    sys.exit()

if not os.path.exists(PROTPATH):
    sys.stderr.write("FATAL: You appear not to have the protein databases. Please download the databases using the download_databases.snakefile\n")
    sys.exit()


# paths for our data. This is where we will read and put things
READDIR = config['Paths']['Reads']
CLUMPED = config['Output']["Clumped"]
QC = config['Output']['QC']
RESULTS = config['Output']['Results']
AA_OUT  = os.path.join(RESULTS, "mmseqs_aa_out")
if not os.path.exists(AA_OUT):
    os.makedirs(AA_OUT, exist_ok=True)

AA_OUT_CHECKED  = os.path.join(RESULTS, "mmseqs_aa_checked_out")
if not os.path.exists(AA_OUT_CHECKED):
    os.makedirs(AA_OUT_CHECKED, exist_ok=True)

VIRDB = os.path.join(PROTPATH, "uniprot_virus_c99.db")
if not os.path.exists(VIRDB):
    sys.stderr.write(f"FATAL: {VIRDB} does not exist. Please ensure you")
    sys.stderr.write(" have installed the databases\n")
    sys.exit()

PHAGE_LINEAGES = os.path.join(DBDIR, "phages", "phage_taxonomic_lineages.txt")
if not os.path.exists(PHAGE_LINEAGES):
    sys.stderr.write("FATAL: phages/phage_taxonomic_lineages.txt not ")
    sys.stderr.write("found in the databases directory. Please check ")
    sys.stderr.write("you have the latest version of the databases\n")
    sys.exit()

# uniref50 + viruses
URVPATH = os.path.join(PROTPATH, "uniref_plus_virus")
URVDB = os.path.join(URVPATH, "uniref50_virus.db") # uniref50 + viruses database
if not os.path.exists(URVDB):
    sys.stderr.write("FATAL: {URVDB} not found.\n")
    sys.stderr.write("Please make sure that you have run ")
    sys.stderr.write("download_databases.snakefile before commencing\n")
    sys.exit()


# how much memory we have
XMX = config['System']['Memory']

SAMPLES, = glob_wildcards(os.path.join(READDIR, '{sample}_R1.fastq.gz'))
PATTERN_R1 = '{sample}_R1'
PATTERN_R2 = '{sample}_R2'


rule all:
    input:
        os.path.join(AA_OUT, "phage_tax_table.tsv"),
        os.path.join(AA_OUT, "viruses_tax_table.tsv"),
        os.path.join(AA_OUT, "pviral_aa_unclassified_seqs.fasta"),
        os.path.join(AA_OUT_CHECKED, "taxonomyResult.firsthit.m8"),
        os.path.join(AA_OUT_CHECKED, "taxonomyResult.report"),
        os.path.join(AA_OUT_CHECKED, "viruses_checked_aa_table.tsv"),
        os.path.join(AA_OUT_CHECKED, "viruses_checked_aa_tax_table.tsv"),
        os.path.join(AA_OUT_CHECKED, "unclassified_checked_aa_seqs.fasta")



"""
Clean the data.

This is from contaminant_removal.snakefile

"""


rule clumpify:
    """
    Step 0: Clumpify and deduplicate reads
    """
    input:
        r1 = os.path.join(READDIR, PATTERN_R1 + ".fastq.gz"),
        r2 = os.path.join(READDIR, PATTERN_R2 + ".fastq.gz")
    output:
        r1 = os.path.join(CLUMPED, PATTERN_R1 + ".clumped.fastq.gz"),
        r2 = os.path.join(CLUMPED, PATTERN_R2 + ".clumped.fastq.gz")
    shell:
        """
        clumpify.sh in={input.r1} in2={input.r2} \
            out={output.r1} out2={output.r2} \
            reorder=a \
            ow=t {XMX}
        """

rule remove_leftmost_primerB:
    """
    Step 1: Remove leftmost primerB. Not the reverse complements
    """
    input:
        r1 = os.path.join(CLUMPED, PATTERN_R1 + ".clumped.fastq.gz"),
        r2 = os.path.join(CLUMPED, PATTERN_R2 + ".clumped.fastq.gz"),
        primers = os.path.join(CONPATH, "primerB.fa")
    output:
        r1 = os.path.join(QC, "step_1", PATTERN_R1 + ".s1.out.fastq"),
        r2 = os.path.join(QC, "step_1", PATTERN_R2 + ".s1.out.fastq"),
        stats = os.path.join(QC, "step_1", "{sample}.s1.stats.txt")
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
    shell:
        """
        bbduk.sh in={input.r1} in2={input.r2} \
            ref={input.primers} \
            out={output.r1} out2={output.r2} \
            stats={output.stats} \
            k=31 hammingdistance=1 ordered=t ow=t {XMX}
        """

rule host_removal:
    """
    Step 6: Host removal
    """
    input:
        r1 = os.path.join(QC, "step_5", PATTERN_R1 + ".s5.out.fastq"),
        r2 = os.path.join(QC, "step_5", PATTERN_R2 + ".s5.out.fastq"),
        refpath = os.path.join(HOSTPATH, "ref")
    output:
        unmapped = os.path.join(QC, "step_6", "{sample}.host.unmapped.s6.out.fastq"),
        mapped = os.path.join(QC, "step_6", "{sample}.host.mapped.s6.out.fastq")
    params:
        hostpath = HOSTPATH 
    shell:
        """
        bbmap.sh in={input.r1} in2={input.r2} \
            outu={output.unmapped} outm={output.mapped} \
            path={params.hostpath} \
            semiperfectmode=t quickmatch fast ordered=t ow=t {XMX}
        """

rule line_sine_removal:
    """
    Step 6a. Remove any LINES and SINES in the sequences.
    """
    input:
        unmapped = os.path.join(QC, "step_6", "{sample}.host.unmapped.s6.out.fastq"),
        linesine = os.path.join(CONPATH, "line_sine.fasta")
    output:
        unmapped = os.path.join(QC, "step_6", "{sample}.linesine.unmapped.s6.out.fastq"),
        mapped   = os.path.join(QC, "step_6", "{sample}.linesine.mapped.s6.out.fastq"),
        stats    = os.path.join(QC, "step_6", "{sample}.linesine.stats")
    shell:
        """
        bbduk.sh in={input.unmapped} out={output.unmapped} \
          outm={output.mapped} \
          ref={input.linesine} k=31 hdist=1 stats={output.stats}
        """

rule repair:
    """
    Step 6b. Repair the paired ends
    """
    input:
        unmapped = os.path.join(QC, "step_6", "{sample}.linesine.unmapped.s6.out.fastq")
    output:
        r1 = os.path.join(QC, "step_6", PATTERN_R1 + ".s6.out.fastq"),
        r2 = os.path.join(QC, "step_6", PATTERN_R2 + ".s6.out.fastq")
    shell:
        """   
        repair.sh in={input.unmapped} \
            out={output.r1} out2={output.r2} \
            ow=t {XMX}
        """

rule trim_low_quality:
    """
    Step 7: Trim low-quality bases
    """
    input:
        r1 = os.path.join(QC, "step_6", PATTERN_R1 + ".s6.out.fastq"),
        r2 = os.path.join(QC, "step_6", PATTERN_R2 + ".s6.out.fastq")
    output:
        r1 = os.path.join(QC, "step_7", PATTERN_R1 + ".s7.out.fastq"),
        r2 = os.path.join(QC, "step_7", PATTERN_R2 + ".s7.out.fastq"),
        singletons = os.path.join(QC, "step_7", "{sample}.singletons.s7.out.fastq"),
        stats = os.path.join(QC, "step_7", "{sample}.s7.stats.txt")
    shell:
        """
        bbduk.sh in={input.r1} in2={input.r2} \
            out={output.r1} out2={output.r2} outs={output.singletons} \
            stats={output.stats} \
            qtrim=r trimq=20 maxns=2 minlength=50 ordered=t ow=t {XMX}
        """

"""
Currently from this point forwards we just work with the R1 reads.
In a previous iteration we merged the reads but that has been removed
"""

rule remove_bacteria:
    """
    Step 8: Remove bacterial contaminants reserving viral and ambiguous sequences
    """
    input:
        r1 = os.path.join(QC, "step_7", PATTERN_R1 + ".s7.out.fastq"),
        bacpath = os.path.join(BACPATH, "ref")
    output:
        mapped = os.path.join(QC, "step_8", "{sample}_R1.bacterial.fastq"),
        unmapped = os.path.join(QC, "step_8", "{sample}_R1.viral_amb.fastq"),
        scafstats = os.path.join(QC, "step_8", "{sample}_R1.scafstats.txt")
    params:
        bacpath = BACPATH
    shell:
        """
        bbmap.sh in={input.r1} \
            path={params.bacpath} \
            outm={output.mapped} outu={output.unmapped} \
            scafstats={output.scafstats} \
            semiperfectmode=t quickmatch fast ordered=t ow=t {XMX}
       """


"""
End the contaminant_removal.snakefile



This section is from cluster_count.snakefile and counts all the clusters!
"""


rule remove_exact_dups:
    """
    Step 9: Remove exact duplicates
    """
    input:
        os.path.join(QC, "step_8", "{sample}_R1.viral_amb.fastq")
    output:
        os.path.join(QC, "step_9", "{sample}_R1.s9.deduped.out.fastq")
    shell:
        """
        dedupe.sh in={input} \
                out={output} \
                ac=f  ow=t {XMX}
        """

rule deduplicate:
    """
    Step 10: Dereplicate
    """
    input:
        os.path.join(QC, "step_9", "{sample}_R1.s9.deduped.out.fastq")
    output:
        fa = os.path.join(QC, "step_10", "{sample}_R1.best.fasta"),
        stats = os.path.join(QC, "step_10", "{sample}_R1.stats.txt")
    shell:
        """
        dedupe.sh in={input} \
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
    shell:
        """
        reformat.sh in={input} out={output} \
            deleteinput=t fastawrap=0 \
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
    params:
        resultsdir = directory(RESULTS),
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
    shell:
        "mmseqs createdb --shuffle 0 --dbtype 0 {input} {output}"

rule seqtable_taxsearch:
    input:
        sq = os.path.join(AA_OUT, "seqtable_query.db"),
    output:
        tr = os.path.join(AA_OUT, "taxonomyResult.dbtype")
    params:
        tr = os.path.join(AA_OUT, "taxonomyResult")
    shell:
        """
        mmseqs taxonomy {input.sq} {VIRDB} {params.tr} $(mktemp -d -p {TMPDIR}) \
        -a --start-sens 1 --sens-steps 3 -s 7 \
        --search-type 2 --tax-output-mode 1
        """

rule seqtable_convert_alignments:
    input:
        sq = os.path.join(AA_OUT, "seqtable_query.db"),
        tr = os.path.join(AA_OUT, "taxonomyResult.dbtype")
    params:
        tr = os.path.join(AA_OUT, "taxonomyResult")
    output:
        os.path.join(AA_OUT, "aln.m8")
    shell:
        """
        mmseqs convertalis {input.sq} {VIRDB} {params.tr} {output} \
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
    shell:
        """
        mmseqs lca {VIRDB} {params.tr} {params.lc} --tax-lineage true \
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
    shell:
        """
        mmseqs createtsv {input.sq} {params.lc} {output}
        """

rule seqtable_create_kraken:
    input:
        lc = os.path.join(AA_OUT, "lca.db")
    output:
        os.path.join(AA_OUT, "taxonomyResult.report")
    shell:
        """
        mmseqs taxonomyreport {VIRDB} {input.lc} {output}
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
    shell:
        """
        mmseqs taxonomy {input.vqdb} {URVDB} {params.tr} \
            $(mktemp -d -p {TMPDIR}) \
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
    shell:
        """
        mmseqs convertalis {input.vqdb} {URVDB} {params.tr} {output} \
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
    shell:
        """
        mmseqs lca {URVDB} {params.tr} {params.lca} \
        --tax-lineage true \
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
    shell:
        """
        mmseqs filterdb {params.tr} {params.fh} --extract-lines 1
        """

rule convertalis_vsqd:
    input:
        vqdb = os.path.join(AA_OUT_CHECKED, "viral_seqs_queryDB"),
        trfhdb = os.path.join(AA_OUT_CHECKED, "taxonomyResult.firsthit.dbtype")
    output:
        os.path.join(AA_OUT_CHECKED, "taxonomyResult.firsthit.m8")
    params:
        trfh = os.path.join(AA_OUT_CHECKED, "taxonomyResult.firsthit")
    shell:
        """
        mmseqs convertalis {input.vqdb} {URVDB} {params.trfh} {output} 
        """

rule create_taxtable_vsqd:
    input:
        vqdb = os.path.join(AA_OUT_CHECKED, "viral_seqs_queryDB"),
        lcadb = os.path.join(AA_OUT_CHECKED, "lca.db.dbtype")
    params:
        lcadb = os.path.join(AA_OUT_CHECKED, "lca.db")
    output:
        os.path.join(AA_OUT_CHECKED, "taxonomyResult.tsv")
    shell:
        """
        mmseqs createtsv {input.vqdb} {params.lcadb} {output}
        """

rule create_kraken_vsqd:
    input:
        lcadb = os.path.join(AA_OUT_CHECKED, "lca.db.dbtype")
    params:
        lcadb = os.path.join(AA_OUT_CHECKED, "lca.db")
    output:
        os.path.join(AA_OUT_CHECKED, "taxonomyResult.report")
    shell:
        """
        mmseqs taxonomyreport {URVDB} {params.lcadb} {output}
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
            sed 's/:/\t/g' | \
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

"""
