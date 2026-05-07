# SSL Analysis for Staphylococcus aureus Vaccine Target Selection

Bioinformatic analysis pipeline for evaluating sequence diversity, phylogenetic coverage, and predicted immunogenicity of selected Staphylococcus aureus Superantigen-Like proteins (SSLs) as candidate vaccine targets.

This repository focuses on the analysis of SSL3, SSL7, and SSL11 sequence variation across global S. aureus isolates using data sourced from the PubMLST S. aureus database.

## Project Overview

### The goal of this project is to:

- Analyze the global sequence diversity of selected S. aureus SSL proteins
- Estimate sequence coverage provided by selected representative strains
- Identify the smallest set of SSL variants capable of maximizing global coverage
- Evaluate predicted immunogenicity using multiple epitope prediction tools
- Support rational vaccine antigen selection for S. aureus

The analysis was performed using representative SSL sequences from major clonal complexes (CCs):

- CC1
- CC3
- CC5
- CC8
- CC22
- CC30

### Repository Structure

SSL_analysis/
├── bin/
├── data/
└── README.md

bin/ — Analysis Scripts

Script | Description
-------|-------------
align_SSL.py |	Align SSL protein/nucleotide sequences
annotate_alignment_peptides.py | Annotate aligned peptide regions
batch_jackhmmer.py | Batch JackHMMER searches
bath_search_SSL_contig.py |	Search SSL-containing contigs
compile_6frame_genomes.py	| Generate six-frame translated genomes
compile_antigens.py	| Compile antigen candidate datasets
compile_pubMLST_loci.py	| Compile loci from PubMLST data
convert_sto_to_fasta.py	| Convert Stockholm alignments to FASTA
esl_alifilter.py	| Alignment filtering wrapper
extract_best_hit_seq.py	| Extract best sequence hits
extract_frame_n.py	| Extract specific translated reading frames
fetch_NCBI_strain_genomes.py	| Download genomes from NCBI
fetch_pubMLST_contigs.py	| Retrieve contigs from PubMLST
fetch_sequences_Uniprot.py	| Download sequences from UniProt
get_random_proteins.py	| Generate random background protein dataset
mafft_align_tree.py	| MAFFT alignment + phylogenetic tree generation
nail_locus_SSLs.py	| Identify SSL loci
plot_immunogenicity_tree.py	| Plot immunogenicity-aware phylogenetic trees
predict_immunogenicity.py	| Run immunogenicity prediction workflows
separate_and_join_fasta.py	| FASTA manipulation utilities
summarize_predictions.py	| Summarize prediction outputs
translate_6_frame.py	| Six-frame translation utility
trim_pubMLST_nail.py	| Trim PubMLST SSL regions
visualize_tree.py	| Tree visualization utilities

data/ — Analysis Data

Directory/File	| Description
----------------|-------------
6_strains/	| Representative SSL sequences from CC1, CC3, CC5, CC8, CC22, and CC30
pubMLST/	| SSL sequence variation sourced from the PubMLST S. aureus database
jackhmmer_results/	| Homology search outputs
ssl/	| SSL-specific processed datasets
SSL_seq.fasta	| SSL protein sequences
SSL_nuc.fasta	| SSL nucleotide sequences
aligned_SSLs.fasta	| Multiple sequence alignment
aligned_SSLs.aln	| Alignment output
aligned_SSLs.dnd	| Guide tree
aligned_SSLs.html	| Alignment visualization

### Analysis Workflow
1. Sequence Collection

Genome contigs from all S. aureus isolates deposited in PubMLST between 2016–2026 were retrieved.

Representative SSL3, SSL7, and SSL11 sequences were mapped against each genome to identify homologous variants.

2. Sequence Clustering

Sequences sharing greater than 90% identity were clustered to reduce redundancy.

One representative sequence per cluster was retained for downstream analysis.

3. Phylogenetic Analysis

Multiple sequence alignments and phylogenetic trees were generated to:

- visualize evolutionary relationships
- estimate global sequence coverage
- identify redundant representatives
- evaluate antigenic diversity

Tree labels contain:

- accession code
- country of origin
- year added
- sequence type
- mapped protein
- coverage

Missing metadata are labeled as "unknown".

Coverage was calculated as:

(cluster sequence count / total retrieved sequences) × 100


## Results
**PubMLST 2016–2026 S. aureus Coverage Trees**

1. The selected clonal complex SSL3 representatives capture:

89.8% of PubMLST SSL3 sequence variation (2016–2026)

*Key observations*

SSL3 CC8 and SSL3 CC1 cluster within the same branch
Strong sequence similarity suggests potential redundancy
A single representative may adequately capture both groups

2. The selected clonal complex SSL7 representatives capture:

93.39% of PubMLST SSL7 sequence variation (2016–2026)

*Key observations*

SSL7 CC5 and SSL7 CC93 cluster together
These representatives may be functionally redundant

3. The selected clonal complex SSL11 representatives capture:

86.17% of PubMLST SSL11 sequence variation (2016–2026)

**Immunogenicity Analysis**

Five bioinformatic prediction tools were used:

- AlgPred2
- BepiPred3
- IFNEpitope2
- NetMHCpan 4.1
- NetMHCIIpan 4.2

Predicted immunogenicity features were normalized using Z-scores generated from a background dataset of:

250 randomly selected S. aureus proteins

Proteins containing an OB-fold domain were included as positive controls.

**Outputs include:**

- immunogenicity heatmaps
- normalized prediction summaries
- comparative antigenicity analyses

*Key Conclusions*

- Selected SSL representatives provide approximately 90% global sequence coverage
- Several representatives display minimal antigenic distance
- Some clonal complex representatives may be removable without major loss of coverage
- Reducing antigen redundancy may simplify vaccine formulation

Examples of potentially redundant representatives include:

- SSL3 CC1 ↔ SSL3 CC8
- SSL3 CC22 ↔ SSL3 CC30

### Future Directions

Potential future improvements include:

- Expanding coverage toward 95% global sequence representation
- Further minimizing representative sequence count
- Integrating structural epitope prediction
- Evaluating experimentally validated immune responses
- Incorporating newer PubMLST genome releases

## Reference:
Gupta et al., 2006

Dependencies

The pipeline relies on several external bioinformatics tools, including:

- MAFFT
- HMMER / JackHMMER
- EMBOSS
- NetMHCpan
- NetMHCIIpan

Additional Python package requirements may include:

- biopython
- pandas
- numpy
- matplotlib
- seaborn
- scipy
- ete3

*Example Usage*

python bin/mafft_align_tree.py \
    --input data/SSL_seq.fasta \
    --output aligned_SSLs
    
python bin/predict_immunogenicity.py \
    --input data/SSL_seq.fasta
 
## Contact

For questions, collaboration, or issues related to the analysis pipeline, please open a GitHub issue or contact nadia.prasetyo@otago.ac.nz
