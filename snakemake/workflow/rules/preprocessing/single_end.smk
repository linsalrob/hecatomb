
# Add preprocessing-specific targets
targets.preprocessing += [
        expand(os.path.join(dir.out.assembly,"{sample}_R1.unmapped.fastq.gz"), sample=samples.names),
        expand(os.path.join(dir.out.assembly,"{sample}_R1.singletons.fastq.gz"), sample=samples.names),
        expand(os.path.join(dir.out.assembly,"{sample}_R1.all.fastq.gz"), sample=samples.names),
    ]


# rules
rule prinseq_trim:
    """Preprocessing step 01: fastp_preprocessing.

    Use fastP to remove adaptors, vector contaminants, low quality sequences, poly-A tails and reads shorts than minimum length, plus deduplicate.
    """
    input:
        r1=lambda wildcards: samples.reads[wildcards.sample]['R1'],
    output:
        r1=temp(os.path.join(dir.out.temp,"p01","{sample}_good_out.fastq")),
        b1=temp(os.path.join(dir.out.temp,"p01","{sample}_bad_out.fastq")),
    benchmark:
        os.path.join(dir.out.bench,"prinseq_trim.{sample}.txt")
    log:
        os.path.join(dir.out.stderr,"prinseq_trim.{sample}.log")
    resources:
        mem_mb=config.resources.sml.mem
    threads:
        config.resources.sml.cpu
    conda:
        os.path.join(dir.env,"prinseqpp.yaml")
    params:
        params=config.prinseq,
        prefix=os.path.join(dir.out.temp,"p01","{sample}")
    shell:
        """
        prinseq++ {params.params} \
          -threads {threads} \
          -out_name {params.prefix} \
          -fastq {input.r1} &> {log}
        rm {log}
        """


rule create_host_index:
    """Step 02. Create the minimap2 index for mapping to the host; this will save time."""
    input:
        dir.dbs.host.fasta,
    output:
        dir.dbs.host.index
    benchmark:
        os.path.join(dir.out.bench, "create_host_index.txt")
    log:
        os.path.join(dir.out.stderr, 'create_host_index.log')
    resources:
        mem_mb = config.resources.med.mem
    threads:
        config.resources.med.cpu
    conda:
        os.path.join(dir.env, "minimap2.yaml")
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
        r1 = os.path.join(dir.out.temp, "p01", "{sample}_good_out.fastq"),
        host = dir.dbs.host.index
    output:
        r1 = temp(os.path.join(dir.out.temp, "p02", "{sample}_R1.unmapped.fastq")),
        s = temp(os.path.join(dir.out.temp, "p02", "{sample}_R1.unmapped.singletons.fastq")),
        o = temp(os.path.join(dir.out.temp, "p02", "{sample}_R1.other.singletons.fastq"))
    benchmark:
        os.path.join(dir.out.bench, "host_removal_mapping.{sample}.txt")
    log:
        mm = os.path.join(dir.out.stderr, "host_removal_mapping.{sample}.minimap.log"),
        sv = os.path.join(dir.out.stderr, "host_removal_mapping.{sample}.samtoolsView.log"),
        fq = os.path.join(dir.out.stderr, "host_removal_mapping.{sample}.samtoolsFastq.log")
    resources:
        mem_mb = config.resources.med.mem
    threads:
        config.resources.med.cpu
    conda:
        os.path.join(dir.env, "minimap2.yaml")
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
        s = os.path.join(dir.out.temp, "p02", "{sample}_R1.unmapped.singletons.fastq"),
        o = os.path.join(dir.out.temp, "p02", "{sample}_R1.other.singletons.fastq")
    output:
        sr1 = temp(os.path.join(dir.out.temp, "p03", "{sample}_R1.u.singletons.fastq")),
        or1 = temp(os.path.join(dir.out.temp, "p03", "{sample}_R1.o.singletons.fastq")),
    benchmark:
        os.path.join(dir.out.bench, "nonhost_read_repair.{sample}.txt")
    log:
        os.path.join(dir.out.stderr, "nonhost_read_repair.{sample}.log")
    resources:
        mem_mb = config.resources.med.mem,
        javaAlloc = int(0.9 * config.resources.med.mem)
    threads:
        config.resources.med.cpu
    conda:
        os.path.join(dir.env, "bbmap.yaml")
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
        r1 = os.path.join(dir.out.temp, "p02", "{sample}_R1.unmapped.fastq"),
        sr1 = os.path.join(dir.out.temp, "p03", "{sample}_R1.u.singletons.fastq"),
        or1 = os.path.join(dir.out.temp, "p03", "{sample}_R1.o.singletons.fastq")
    output:
        t1 = temp(os.path.join(dir.out.temp, "p04", "{sample}_R1.singletons.fastq")),
        r1 = temp(os.path.join(dir.out.temp, "p04", "{sample}_R1.all.fastq"))
    benchmark:
        os.path.join(dir.out.bench, "nonhost_read_combine.{sample}.txt")
    log:
        os.path.join(dir.out.stderr, "nonhost_read_combine.{sample}.log")
    shell:
        """
        {{ cat {input.sr1} {input.or1} > {output.t1};
        cat {input.r1} {output.t1} > {output.r1}; }} 2> {log}
        rm {log}
        """


rule archive_for_assembly:
    """Copy the files that will be required in the assembly steps; fastq.gz files will be generated from these"""
    input:
        os.path.join(dir.out.temp,"p02","{sample}_R1.unmapped.fastq"),
        os.path.join(dir.out.temp,"p04","{sample}_R1.singletons.fastq"),
        os.path.join(dir.out.temp,"p04","{sample}_R1.all.fastq"),
    output:
        temp(os.path.join(dir.out.assembly,"{sample}_R1.unmapped.fastq")),
        temp(os.path.join(dir.out.assembly,"{sample}_R1.singletons.fastq")),
        temp(os.path.join(dir.out.assembly,"{sample}_R1.all.fastq")),
    params:
        dir.out.assembly
    shell:
        """cp {input} {params}"""
