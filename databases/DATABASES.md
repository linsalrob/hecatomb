# Databases used in hecatomb

We use a variety of databases that we have catologed here:

## Bacterial databases

*Database*: bac_uniquespecies_giant.masked_Ns_removed.fasta
*md5sum*: f7671c5507974903979295c586c729f0
*dbtype*: nucl
*Sequences*: 57,161
*bp*: 2,811,399,982
*Genomes*: 930
*Source*: Bacterial genomes from the gut masked to remove prophages
*Comments*:

Here are the most abundant genera in this dataset:

Genus | Genomes
--- | ---
Ruminococcus | 127
Clostridium | 87
Bacteroides | 84
Collinsella | 78
Coprobacillus | 39
Blautia | 37
Roseburia | 27
Eubacterium | 26
Parabacteroides | 23
Lachnospiraceae | 23
Firmicutes | 19
Coprococcus | 18
Butyricicoccus | 17
Streptococcus | 12
Clostridiaceae | 10


## Contaminant databases

*Database*:  nebnext_adapters.fa
*md5sum*: 96967a1372600346cf93126cb197a206
*dbtype*: nucl
*Sequences*: 49
*bp*: 3,153
*Source*:  NEBNext Multiplex Oligos for Illumina. 
*Comments*: The NEBNext Ultra Universal Primer plus 48 indexed primer sequences

*Database*:   primerB.fa
*md5sum*: e83366d892407f6f3a7244595f2e8a28
*dbtype*:  nucl
*Sequences*: 24
*bp*: 384
*Source*: The roundAB protocol
*Comments*: These are the 24 primers used in the roundAB protocol for amplifying viral DNA prior to sequencing.

*Database*:  rc_primerB_ad6.fa
*md5sum*: fb7a711e118bd8139ee652910c5ebcf5
*dbtype*: nucl
*Sequences*: 24
*bp*: 528
*Source*: The roundAB protocol
*Comments*: This is the reverse complement of each primer in primerB.fa with the AGATCG tag appended at the 3' end

*Database*:  line_sine.fasta
*md5sum*: 99c5e57d1f2040a7f392598726cf20bf
*dbtype*: nucl
*Sequences*: 238
*bp*: 60,526 
*Source*: A combination of http://sines.eimb.ru/banks/SINEs.bnk and http://sines.eimb.ru/banks/LINEs.bnk
*Comments*: A LINE/SINE database to remove those from the sequences. A few euk viruses hit strongly to LINES/SINES


## Human genomes

*Database*:  human_virus_masked.fasta
*md5sum*: 86dce9ed9cc52b0f9176b78fdfbf031e
*dbtype*: nucl
*Sequences*: 25
*bp*: 2,818,003,269
*Source*: Started with GRCh38 human reference (masked) we remasked against all known viruses.
*Comments*:


# Phage

*Database*:  phage_taxonomic_lineages.txt
*md5sum*: 26973c874c8718990a40cd6a759b58d7
*dbtype*: text
*Sequences*: -
*bp*: -
*Source*: A list of taxonomic levels that are purely phage. 
*Comments*:

# Proteins


uniref50_virus.fasta: concatentation of uniref50.fasta.gz and uniprot_virus_c99.faa


*Database*: uniref50.fasta.gz 
*md5sum*: 
*dbtype*: prot
*Sequences*: 39,992,926
*aa*: 11,307,670,958
*Source*: Downloaded from uniref: ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz
*Comments*: All uniref proteins clustered at 50% identity

*Database*: uniprot_virus_c99.faa
*md5sum*:
*dbtype*: prot
*Sequences*: 1,745,139
*aa*: 518,718,189
*Source*: Uniprot Viruses Clustered at 99%
*Comments*: 

All the uniprot viruses are downloaded using the [https://www.uniprot.org/uniprot/?query=taxonomy:%22Viruses%20[10239]%22&format=fasta&&sort=score&fil=reviewed:no](https://www.uniprot.org/uniprot/?query=taxonomy:%22Viruses%20[10239]%22&format=fasta&&sort=score&fil=reviewed:no) URL, and then clustered at 99% identity using `cd-hit`

The unclustered file has 4,383,357 sequences and 1,377,032,930 amino acids, so the clustered file is about 1/3 the size of the original (file sizes 688M vs. 1.6 G)

*Database*: uniref50_virus.fasta
*md5sum*:
*dbtype*: prot
*Sequences*:
*aa*:
*Source*: Concatentation of uniref50.fasta.gz and uniprot_virus_c99.faa
*Comments*: This is a temporary file that is currently deleted upon compilation into a mmseqs database

