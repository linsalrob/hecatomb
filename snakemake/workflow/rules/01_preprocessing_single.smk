"""
Snakemake rule file to preprocess Illumina sequence data for virome analysis.

What is accomplished with these rules?
    - Non-biological sequence removal (primers, adapters)
    - Host sequence removal
    - Removal of redundant sequences (duplicates + clustering)
        - Creation of sequence count table
        - Calculation of sequence properties (e.g. GC content, tetramer frequencies)

Rob Edwards, Jan 2020
Updated: Scott Handley, March 2021
Updated: Michael Roach, Q2/3 2021
Updated: Sarah Beecroft Q1 2022
"""

rule fastp_preprocessing:
    """Preprocessing step 01: fastp_preprocessing.
    
    Use fastP to remove adaptors, vector contaminants, low quality sequences, poly-A tails and reads shorts than minimum length, plus deduplicate.
    """
    input:
        r1 = lambda wildcards: sampleReads[wildcards.sample]['R1'],
        contaminants = os.path.join(CONPATH, "vector_contaminants.fa"),
        summ = optionalSummary[0]
    output:
        r1 = temp(os.path.join(TMPDIR, "p01", "{sample}_R1.s1.out.fastq")),
        stats = os.path.join(STATS, "p01", "{sample}.s1.stats.json")
    benchmark:
        os.path.join(BENCH, "fastp_preprocessing.{sample}.txt")
    log:
        os.path.join(STDERR, "fastp_preprocessing.{sample}.log")
    resources:
        mem_mb = FastpMem
    threads:
        FastpCPU
    conda:
        "../envs/fastp.yaml"
    shell:
        """
        fastp -i {input.r1} -o {output.r1} \
            -z {config[COMPRESSION]} \
            -j {output.stats} \
            --qualified_quality_phred {config[QSCORE]} \
            --length_required {config[READ_MINLENGTH]} \
            --adapter_fasta {input.contaminants} \
            --cut_tail --cut_tail_window_size {config[CUTTAIL_WINDOW]} --cut_tail_mean_quality {config[QSCORE]} \
            --dedup --dup_calc_accuracy {config[DEDUP_ACCURACY]} \
            --trim_poly_x \
            --thread {threads} 2> {log}
            rm {log}
        """

rule create_host_index:
    """Step 02. Create the minimap2 index for mapping to the host; this will save time."""
    input:
        HOSTFA,
    output:
        HOSTINDEX
    benchmark:
        os.path.join(BENCH, "create_host_index.txt")
    log:
        os.path.join(STDERR, 'create_host_index.log')
    resources:
        mem_mb = BBToolsMem
    threads:
        BBToolsCPU
    conda:
        "../envs/minimap2.yaml"
    shell:
        """
        minimap2 -t {threads} -d {output} <(cat {input}) 2> {log}
        rm {log}
        """

rule host_removal_mapping:
    """Preprocessing step 02a: Host removal: mapping to host.
    
    Must define host in config file (see Paths: Host: in config.yaml). Host should be masked of viral sequence.
    If your reference is not available you need to add it using 'Hecatomb addHost'
    """
    input:
        r1 = os.path.join(TMPDIR, "p01", "{sample}_R1.s1.out.fastq"),
        host = HOSTINDEX
    output:
        r1 = temp(os.path.join(TMPDIR, "p02", "{sample}_R1.unmapped.fastq")),
        s = temp(os.path.join(TMPDIR, "p02", "{sample}_R1.unmapped.singletons.fastq")),
        o = temp(os.path.join(TMPDIR, "p02", "{sample}_R1.other.singletons.fastq"))
    benchmark:
        os.path.join(BENCH, "host_removal_mapping.{sample}.txt")
    log:
        mm = os.path.join(STDERR, "host_removal_mapping.{sample}.minimap.log"),
        sv = os.path.join(STDERR, "host_removal_mapping.{sample}.samtoolsView.log"),
        fq = os.path.join(STDERR, "host_removal_mapping.{sample}.samtoolsFastq.log")
    resources:
        mem_mb = BBToolsMem
    threads:
        BBToolsCPU
    conda:
        "../envs/minimap2.yaml"
    shell:
        """
        minimap2 -ax sr -t {threads} --secondary=no {input.host} {input.r1} 2> {log.mm} \
            | samtools view -f 4 -h 2> {log.sv} \
            | samtools fastq -NO -1 {output.r1} -0 {output.o} -s {output.s} 2> {log.fq}
        rm {log.mm} {log.sv} {log.fq}
        """

rule nonhost_read_repair:
    """Preprocessing step 03: Parse R1/R2 singletons (if singletons at all)"""
    input:
        s = os.path.join(TMPDIR, "p02", "{sample}_R1.unmapped.singletons.fastq"),
        o = os.path.join(TMPDIR, "p02", "{sample}_R1.other.singletons.fastq")
    output:
        sr1 = temp(os.path.join(TMPDIR, "p03", "{sample}_R1.u.singletons.fastq")),
        or1 = temp(os.path.join(TMPDIR, "p03", "{sample}_R1.o.singletons.fastq")),
    benchmark:
        os.path.join(BENCH, "nonhost_read_repair.{sample}.txt")
    log:
        os.path.join(STDERR, "nonhost_read_repair.{sample}.log")
    resources:
        mem_mb = BBToolsMem,
        javaAlloc = int(0.95 * BBToolsMem)
    threads:
        BBToolsCPU
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        {{ reformat.sh in={input.s} out={output.sr1} \
            -Xmx{resources.javaAlloc}m;
        reformat.sh in={input.o} out={output.or1} \
            -Xmx{resources.javaAlloc}m; }} 2>> {log}
        rm {log}
        """

rule nonhost_read_combine:
    """Preprocessing step 04: Combine paired and singleton reads. """
    input:
        r1 = os.path.join(TMPDIR, "p02", "{PATTERN}_R1.unmapped.fastq"),
        sr1 = os.path.join(TMPDIR, "p03", "{PATTERN}_R1.u.singletons.fastq"),
        or1 = os.path.join(TMPDIR, "p03", "{PATTERN}_R1.o.singletons.fastq")
    output:
        t1 = os.path.join(TMPDIR, "p04", "{PATTERN}_R1.singletons.fastq"),
        r1 = os.path.join(TMPDIR, "p04", "{PATTERN}_R1.all.fastq")
    benchmark:
        os.path.join(BENCH, "nonhost_read_combine.{PATTERN}.txt")
    log:
        os.path.join(STDERR, "nonhost_read_combine.{PATTERN}.log")
    shell:
        """
        {{ cat {input.sr1} {input.or1} > {output.t1};
        cat {input.r1} {output.t1} > {output.r1}; }} 2> {log}
        rm {log}
        """
          
rule cluster_similar_sequences:
    """Preprocessing step 05: Cluster similar sequences.
     
     Sequences clustered at CLUSTERID in config.yaml.
    """
    input:
        fq = os.path.join(TMPDIR, "p04", "{sample}_R1.all.fastq"),
        summ = optionalSummary[1]
    output:
        temp(os.path.join(TMPDIR, "p05", "{sample}_R1_rep_seq.fasta")),
        temp(os.path.join(TMPDIR, "p05", "{sample}_R1_cluster.tsv")),
        temp(os.path.join(TMPDIR, "p05", "{sample}_R1_all_seqs.fasta"))
    params:
        respath = os.path.join(TMPDIR, "p05"),
        tmppath = os.path.join(TMPDIR, "p05", "{sample}_TMP"),
        prefix = '{sample}_R1',
        config = config['linclustParams']
    benchmark:
        os.path.join(BENCH, "cluster_similar_sequences.{sample}.txt")
    log:
        os.path.join(STDERR, "cluster_similar_sequences.{sample}.log")
    resources:
        mem_mb = MMSeqsMem
    threads:
        MMSeqsCPU
    conda:
        "../envs/mmseqs2.yaml"
    shell:
        """ 
        mmseqs easy-linclust {input.fq} {params.respath}/{params.prefix} {params.tmppath} \
            {params.config} \
            --threads {threads} &> {log}
        rm {log}
        """
        
rule create_individual_seqtables:
    """Preprocessing step 06: Create individual seqtables. 
    
    A seqtable is a count table with each sequence as a row, each column as a sample and each cell the counts of each 
    sequence per sample.
    """
    input:
        seqs = os.path.join(TMPDIR, "p05", "{sample}_R1_rep_seq.fasta"),
        counts = os.path.join(TMPDIR, "p05", "{sample}_R1_cluster.tsv"),
        summ = optionalSummary[2]
    output:
        seqs = temp(os.path.join(TMPDIR, "p06", "{sample}_R1.seqs")),
        counts = temp(os.path.join(TMPDIR, "p06", "{sample}_R1.counts")),
        seqtable = temp(os.path.join(TMPDIR, "p06", "{sample}_R1.seqtable"))
    benchmark:
        os.path.join(BENCH, "individual_seqtables.{sample}.txt")
    log:
        os.path.join(STDERR, "individual_seqtables.{sample}.txt")
    resources:
        mem_mb = MMSeqsMem
    threads:
        MMSeqsCPU
    conda:
        "../envs/seqkit.yaml"
    shell:
        """
        {{ seqkit sort {input.seqs} --quiet -j {threads} -w 5000 -t dna \
            | seqkit fx2tab -w 5000 -t dna \
            | sed 's/\\t\\+$//' \
            | cut -f2,3 \
            | sed '1i sequence' > {output.seqs};
        cut -f1 {input.counts} \
            | sort \
            | uniq -c \
            | awk -F ' ' '{{print$2"\\t"$1}}' \
            | cut -f2 \
            | sed "1i {wildcards.sample}" > {output.counts};
        paste {output.seqs} {output.counts} > {output.seqtable}; }} 2> {log}
        rm {log}
        """


rule merge_seq_table:
    """Preprocessing step 07: Merge seq tables
    
    Reads the sequences and counts from each samples' seqtable text file and converts to fasta format for the rest of 
    the pipline.
    """
    input:
        seqtable = expand(os.path.join(TMPDIR, "p06", "{sample}_R1.seqtable"), sample=SAMPLES)
    output:
        fa = os.path.join(RESULTS, "seqtable.fasta"),
        tsv = os.path.join(RESULTS, "sampleSeqCounts.tsv")
    params:
        samples = list(SAMPLES),
        tmpdir = os.path.join(TMPDIR, 'p06')
    #conda:
    #    os.path.join('..', 'envs', 'pysam.yaml')
    benchmark:
        os.path.join(BENCH, "merge_seq_table.txt")
    log:
        os.path.join(STDERR, 'merge_seq_table.log')
    script:
        os.path.join('../', 'scripts', 'mergeSeqTable.py')
