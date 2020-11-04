"""

Snakefile based on [contaminant_removal.sh](../base/contaminant_removal.sh)

Rob Edwards, Jan 2020
Updated: Scott Handley, Nov 2020
"""

import os

# NOTE: bbtools uses "threads=auto" by default that typically uses all threads, so no need to specify. 
# -Xmx is used to specify the memory allocation for bbtools operations

rule remove_leftmost_primerB:
    """
    Step 01: Remove leftmost primerB. Not the reverse complements
    """
    input:
        r1 = os.path.join(READDIR, PATTERN_R1 + file_extension),
        r2 = os.path.join(READDIR, PATTERN_R2 + file_extension),
        primers = os.path.join(CONPATH, "primerB.fa")
    output:
        r1 = temp(os.path.join(TMPDIR, "step_01", PATTERN_R1 + ".s1.out.fastq")),
        r2 = temp(os.path.join(TMPDIR, "step_01", PATTERN_R2 + ".s1.out.fastq")),
        stats = os.path.join(STATS, "step_01", "{sample}.s1.stats.tsv")
    benchmark:
        "benchmarks/removeprimerB_{sample}.txt"
    log:
        "LOGS/step_01/{sample}.s1.log"
    resources:
        mem_mb=20000,
        cpus=8
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        bbduk.sh in={input.r1} in2={input.r2} \
            ref={input.primers} \
            out={output.r1} out2={output.r2} \
            stats={output.stats} \
            k=16 hdist=1 mink=11 ktrim=l restrictleft=20 \
            removeifeitherbad=f trimpolya=10 ordered=t rcomp=f ow=t \
            -Xmx{resources.mem_mb}m 2> {log}
        """

rule remove_3prime_contaminant:
    """
    Step 02: Remove 3' read through contaminant
    """
    input:
        r1 = os.path.join(TMPDIR, "step_01", PATTERN_R1 + ".s1.out.fastq"),
        r2 = os.path.join(TMPDIR, "step_01", PATTERN_R2 + ".s1.out.fastq"),
        primers = os.path.join(CONPATH, "rc_primerB_ad6.fa")
    output:
        r1 = temp(os.path.join(TMPDIR, "step_02", PATTERN_R1 + ".s2.out.fastq")),
        r2 = temp(os.path.join(TMPDIR, "step_02", PATTERN_R2 + ".s2.out.fastq")),
        stats = os.path.join(STATS, "step_02", "{sample}.s2.stats.tsv")
    benchmark:
        "benchmarks/remove_3prime_contaminant_{sample}.txt"
    log:
        "LOGS/step_02/{sample}.s2.log"
    resources:
        mem_mb=20000,
        cpus=8
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        bbduk.sh in={input.r1} in2={input.r2} \
            ref={input.primers} \
            out={output.r1} out2={output.r2} \
            stats={output.stats} \
            k=16 hdist=1 mink=11 ktrim=r removeifeitherbad=f ordered=t rcomp=f ow=t \
            -Xmx{resources.mem_mb}m 2> {log}
        """

rule remove_primer_free_adapter:
    """
    Step 03: Remove primer free adapter (both orientations)
    """
    input:
        r1 = os.path.join(TMPDIR, "step_02", PATTERN_R1 + ".s2.out.fastq"),
        r2 = os.path.join(TMPDIR, "step_02", PATTERN_R2 + ".s2.out.fastq"),
        primers = os.path.join(CONPATH, "nebnext_adapters.fa")
    output:
        r1 = temp(os.path.join(TMPDIR, "step_03", PATTERN_R1 + ".s3.out.fastq")),
        r2 = temp(os.path.join(TMPDIR, "step_03", PATTERN_R2 + ".s3.out.fastq")),
        stats = os.path.join(STATS, "step_03", "{sample}.s3.stats.tsv")
    benchmark:
        "benchmarks/remove_primer_free_adapter_{sample}.txt"
    log:
        "LOGS/step_03/{sample}.s3.log"
    resources:
        mem_mb=20000,
        cpus=8
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        bbduk.sh in={input.r1} in2={input.r2} \
            ref={input.primers} \
            out={output.r1} out2={output.r2} \
            stats={output.stats} \
            k=16 hdist=1 mink=10 ktrim=r removeifeitherbad=f ordered=t rcomp=t ow=t \
            -Xmx{resources.mem_mb}m 2> {log}
        """

rule remove_adapter_free_primer:
    """
    Step 04: Remove adapter free primer (both orientations)
    """
    input:
        r1 = os.path.join(TMPDIR, "step_03", PATTERN_R1 + ".s3.out.fastq"),
        r2 = os.path.join(TMPDIR, "step_03", PATTERN_R2 + ".s3.out.fastq"),
        primers = os.path.join(CONPATH, "rc_primerB_ad6.fa")
    output:
        r1 = temp(os.path.join(TMPDIR, "step_04", PATTERN_R1 + ".s4.out.fastq")),
        r2 = temp(os.path.join(TMPDIR, "step_04", PATTERN_R2 + ".s4.out.fastq")),
        stats = os.path.join(STATS, "step_04", "{sample}.s4.stats.tsv")
    benchmark:
        "benchmarks/remove_adapter_free_primer_{sample}.txt"
    log:
        "LOGS/step_04/{sample}.s4.log"
    resources:
        mem_mb=20000,
        cpus=8
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        bbduk.sh in={input.r1} in2={input.r2} \
            ref={input.primers} \
            out={output.r1} out2={output.r2} \
            stats={output.stats} \
            k=16 hdist=0 removeifeitherbad=f ordered=t rcomp=t ow=t \
            -Xmx{resources.mem_mb}m 2> {log}
        """

rule remove_vector_contamination:
    """
    Step 05: Vector contamination removal (PhiX + NCBI UniVecDB)
    """
    input:
        r1 = os.path.join(TMPDIR, "step_04", PATTERN_R1 + ".s4.out.fastq"),
        r2 = os.path.join(TMPDIR, "step_04", PATTERN_R2 + ".s4.out.fastq"),
        primers = os.path.join(CONPATH, "vector_contaminants.fa")
    output:
        r1 = temp(os.path.join(TMPDIR, "step_05", PATTERN_R1 + ".s5.out.fastq")),
        r2 = temp(os.path.join(TMPDIR, "step_05", PATTERN_R2 + ".s5.out.fastq")),
        stats = os.path.join(STATS, "step_05", "{sample}.s5.stats.tsv")
    benchmark:
        "benchmarks/remove_vector_contamination_{sample}.txt"
    log:
        "LOGS/step_05/{sample}.s5.log"
    resources:
        mem_mb=20000,
        cpus=8
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        bbduk.sh in={input.r1} in2={input.r2} \
            ref={input.primers} \
            out={output.r1} out2={output.r2} \
            stats={output.stats} \
            k=31 hammingdistance=1 ordered=t ow=t \
            -Xmx{resources.mem_mb}m 2> {log}
        """
        
rule remove_low_quality:
    """
    Step 06: Remove remaining low-quality bases and short reads
    """
    input:
        r1 = os.path.join(TMPDIR, "step_05", PATTERN_R1 + ".s5.out.fastq"),
        r2 = os.path.join(TMPDIR, "step_05", PATTERN_R2 + ".s5.out.fastq")
    output:
        r1 = os.path.join(QC, PATTERN_R1 + ".clean.out.fastq"),
        r2 = os.path.join(QC, PATTERN_R2 + ".clean.out.fastq"),
        stats = os.path.join(STATS, "step_06", "{sample}.s6.stats.tsv")
    benchmark:
        "benchmarks/remove_low_quality_{sample}.txt"
    log:
        "LOGS/step_06/{sample}.s6.log"
    resources:
        mem_mb=20000,
        cpus=8
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        bbduk.sh in={input.r1} in2={input.r2} \
            out={output.r1} out2={output.r2} \
            stats={output.stats} \
            ordered=t \
            qtrim=r maxns=2 \
            trimq={config[QSCORE]} \
            minlength={config[MINLENGTH]} 2> {log} 
        """
