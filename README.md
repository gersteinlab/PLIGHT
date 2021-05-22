# PLIGHT
**P**rivacy **L**eakage through **I**nference across **G**enotypic **H**MM **T**rajectories

<img src="https://github.com/gersteinlab/PLIGHT/blob/main/images/PLIGHT_github_figure.png" width="600">

**PLIGHT** is a computational tool that employs a population-genetics-based Hidden Markov Model of recombination and mutation to find piecewise matches of a sparse query set of SNPs to a reference genotype panel. The premise of the tool is that even limited, noisy and sparsely distributed genotypic information carries with it a certain risk of identification and downstream inference.

Inspired by imputation methods such as *IMPUTE2* [[1]](#1) and *Eagle* [[2]](#2), the inference procedure in **PLIGHT** is based on the Li-Stephens model [[3]](#3), where an HMM is used to explore the space of underlying pairs of haplotypes in a diploid genome with the possibility of de novo mutations and recombination between haplotypes. A solution to the inference problem consists of a set of best-fit haplotype pairs at each observed locus, each pair being linked to another pair at the next locus, to form a set of piecewise matches to reference haplotypes. If multiple equally likely solutions exist, the method identifies all of them.  Collectively, these form a set of genotypic trajectories through reference haplotype space, where a trajectory is defined as a sequence of reference haplotype pairs (for a diploid genome) at each locus that best fit the observations.

For further details about the method and application cases, please refer to:
> Emani, P.S.; Gürsoy, G.; Miranker, A.; Gerstein, M.B. PLIGHT: A tool to assess privacy risk by inferring identifying characteristics from sparse, noisy genotypes. **2021** *biorxiv* **TBD**

## Code Usage
The code is written in Python 3, and consists of a set of three algorithms with special use cases:
1. **PLIGHT_Exact** performs the exact HMM inference process using the Viterbi algorithm [[4]](#4);
2. **PLIGHT_Truncated** phases in a process of truncating the set of all calculated trajectories to only those within a certain probability distance from the maximally optimal ones, resulting in a smaller memory footprint;
3. **PLIGHT_Iterative** iteratively partitions the reference search space into more manageable blocks of haplotypes and runs **PLIGHT_Exact** on each block, followed by pooling and repetition of the scheme on the resulting, smaller cohort of haplotypes.


 
## References
<a id="1">[1]</a>
Howie, B. N.; Donnelly, P.; Marchini, J. A Flexible and Accurate Genotype Imputation Method for the next Generation of Genome-Wide Association Studies. PLoS Genet. 2009, 5 (6). https://doi.org/10.1371/journal.pgen.1000529.

<a id="2">[2]</a>
Loh, P. R.; Danecek, P.; Palamara, P. F.; Fuchsberger, C.; Reshef, Y. A.; Finucane, H. K.; Schoenherr, S.; Forer, L.; McCarthy, S.; Abecasis, G. R.; et al. Reference-Based Phasing Using the Haplotype Reference Consortium Panel. *Nat. Genet.* **2016**, 48 (11), 1443–1448. https://doi.org/10.1038/ng.3679.

<a id="3">[3]</a>
Li, N.; Stephens, M. Modeling Linkage Disequilibrium and Identifying Recombination Hotspots Using Single-Nucleotide Polymorphism Data. *Genetics* **2003**, 165, 2213–2233.

<a id="4">[4]</a>
Viterbi, A. Error Bounds for Convolutional Codes and an Asymptotically Optimum Decoding Algorithm. *IEEE Trans. Inf. Theory* **1967**, 13 (2), 260–269. https://doi.org/10.1109/TIT.1967.1054010
