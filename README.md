# PLIGHT
**P**rivacy **L**eakage through **I**nference across **G**enotypic **H**MM **T**rajectories

<img src="https://github.com/gersteinlab/PLIGHT/blob/main/images/PLIGHT_github_figure.png" width="600">

**PLIGHT** is a computational tool that employs a population-genetics-based Hidden Markov Model of recombination and mutation to find piecewise matches of a sparse query set of SNPs to a reference diploid genotype panel. The premise of the tool is that even limited, noisy and sparsely distributed genotypic information carries with it a certain risk of identification and downstream inference.

Inspired by imputation methods such as *IMPUTE2* [[1]](#1) and *Eagle* [[2]](#2), the inference procedure in **PLIGHT** is based on the Li-Stephens model [[3]](#3), where an HMM is used to explore the space of underlying pairs of haplotypes in a diploid genome with the possibility of de novo mutations and recombination between haplotypes. A solution to the inference problem consists of a set of best-fit haplotype pairs at each observed locus, each pair being linked to another pair at the next locus, to form a set of piecewise matches to reference haplotypes. If multiple equally likely solutions exist, the method identifies all of them.  Collectively, these form a set of genotypic trajectories through reference haplotype space, where a trajectory is defined as a sequence of reference haplotype pairs (for a diploid genome) at each locus that best fit the observations.

For further details about the method and application cases, please refer to:
> Emani, P.S.; Gürsoy, G.; Miranker, A.; Gerstein, M.B. PLIGHT: A tool to assess privacy risk by inferring identifying characteristics from sparse, noisy genotypes. **2021** *biorxiv* **doi:** https://doi.org/10.1101/2021.07.18.452853

## Software Requirements
### OS requirements
We tested the code in the Red Hat Enterprise Linux Server 7.9 (Maipo) environment. We have tested the batch job submission using the Slurm job handler. An example of a Slurm script for submission is included. We highly recommend the deployment of multiple CPUs in a job, as the code is designed to parallelize the HMM optimization.

### External Dependencies
The following external programs employed by the algorithms and need to be installed before running the code:
```
bcftools (version 1.10.2)
tabix (version 1.11)
```

### Python Requirements
The libraries/modules required in the corresponding Python scripts are:
```
numpy (several versions work, we have used 1.18 and 1.19 at different stages of development)
subprocess
gzip
mmap
multiprocessing
```
## Installation guide
```
git clone https://github.com/gersteinlab/PLIGHT.git
```

## External data used in manuscript 
1. The genotype vcf files used as the reference database can be downloaded from the 1000 Genomes Consortium: http://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/
2. The identification of related individuals is based on the 1000 Genomes Phase 3 Related Samples Panel: https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/related_samples_vcf/related_samples_panel.20140910.ALL.panel.
3. 1000 Genomes Phase 3 related individuals pedigree file https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped.

## Code Description
The code is written in Python 3, and consists of a set of three algorithms with special use cases:
1. **PLIGHT_Exact** performs the exact HMM inference process using the Viterbi algorithm [[4]](#4);
    - **PLIGHT_InRef** is a specific individual-in-the-reference-database instance of **PLIGHT_Exact**, where the recombination rate is set to 0; that is, all SNPs are assumed to belong to one individual, and the most likely individual(s) in the reference database is(are) found.
2. **PLIGHT_Truncated** phases in a process of truncating the set of all calculated trajectories to only those within a certain probability distance from the maximally optimal ones, resulting in a smaller memory footprint;
3. **PLIGHT_Iterative** iteratively partitions the reference search space into more manageable blocks of haplotypes and runs **PLIGHT_Exact** on each block, followed by pooling and repetition of the scheme on the resulting, smaller cohort of haplotypes.


Below is an example of all the arguments that can be passed to the PLIGHT_Exact module:

```
usage: PLIGHT_Exact.py [-h] -c CHROMOSOMEFILE [-O OBSERVEDSAMPLE]
                       [-I CHROMOSOMEID] [-m METADATA] [-F GENFOLDER]
                       [-M MUTATIONRATE] [-e EFFPOP] [-r REFPOP]
                       [-b RECOMBRATE] [-s {distance,custom}] [-f RECOMBFILE]
                       [-t TOLERANCE] [-C CURRDIR] [--affilter AFFILTER]
                       [--numproc NUMPROC]

Identify closest related reference haplotypes

optional arguments:
  -h, --help            show this help message and exit
  -c CHROMOSOMEFILE, --chromosomefile CHROMOSOMEFILE
                        Chromosome file name
  -O OBSERVEDSAMPLE, --observedsample OBSERVEDSAMPLE
                        Observed Sample Genotype File
  -I CHROMOSOMEID, --chromosomeID CHROMOSOMEID
                        Chromosome ID
  -m METADATA, --metadata METADATA
                        Metadata file with ancestry information
  -F GENFOLDER, --genfolder GENFOLDER
                        Genotype folder
  -M THETAMUTATIONRATE, --thetamutationrate THETAMUTATIONRATE
                        Theta (Coalescent) Mutation rate
  -L LAMBDAMUTATIONRATE, --lambdamutationrate LAMBDAMUTATIONRATE
                        Lambda (direct error) Mutation rate
  -e EFFPOP, --effpop EFFPOP
                        Effective population
  -r REFPOP, --refpop REFPOP
                        Reference population
  -b RECOMBRATE, --recombrate RECOMBRATE
                        Recombination rate in cM/Mb (if /distance/ option is
                        chosen)
  -s {distance,custom}, --recombswitch {distance,custom}
                        Recombination Model Switch (distance = simple linear
                        increase of recombination rate with distance, custom =
                        alternative model)
  -f RECOMBFILE, --recombfile RECOMBFILE
                        Custom Recombination Model File
  -t TOLERANCE, --tolerance TOLERANCE
                        Fraction of maximum value used as allowance for
                        inclusion of a path
  -C CURRDIR, --currdir CURRDIR
                        Current working directory
  --numproc NUMPROC     Number of processes for parallelization
  --posspecific POSSPECIFIC
                        Position-specific mutation rates included in
                        observation file? (True/False)
  --prefix PREFIX       String prefix to append to output Best_trajectories
                        file, in addition to chromosome number
```
Additional options for the other modules include

```
PLIGHT_Truncated.py
--truncation TRUNCATION
                        Fraction of trajectories carried forward from each
                        step

PLIGHT_Iterative.py
--subgroup SUBGROUP   Number of haplotypes in each iterative scheme subgroup
--niter NITER         Number of iterations of bootstrapping process
```

In the case of ```PLIGHT_InRef.py```, the recombination rate parameters are all ignored.

## Explanation of the parameters
1. -c CHROMOSOMEFILE, --chromosomefile CHROMOSOMEFILE : This is the name of the reference database vcf file, either a composite of all chromosomes or separated out by chromosome; note that the code only runs the model for a single chromosome at a time. 
2. -O OBSERVEDSAMPLE, --observedsample OBSERVEDSAMPLE : This file contains the set of observed SNPs in tab-delimited format (no header should be included) with columns of 
   >```Chromosome_Number Genome_Position Genome_Position Alternate_Allele Observed_Genotype (if POSSPECIFIC=True)Alternative_Genotype_1:Probability of Observing this Genotype (if POSSPECIFIC=True)Alternative_Genotype_2:Probability of Observing this Genotype```
   
    For example, let us assume we have a query SNP from chromosome 5 at position 33951693, with an observed alternate allele of ```G``` and called genotype of ```1```. If the user chooses to pass a position-specific error rate (as indicated by the ```POSSPECIFIC``` flag described below), with ```98%``` probability of the called genotype and ```2%``` probability of a genotype ```2```, they would enter the line as 
   >```5       33951693        33951693        G       1       0:0.0   2:0.02```
  
   Note that in this example, the genotype of ```0``` is thought to have ```0``` probability, and all the probabilities sum to 1.
  
3.  -I CHROMOSOMEID, --chromosomeID CHROMOSOMEID : This is the chromosome being analyzed, with the format 'chr\<Chromosome Number\>'. This is appended as a prefix to the output files.

4. -F GENFOLDER, --genfolder GENFOLDER : Folder location of the reference genotype vcf file, with a trailing '/' character
5. -M THETAMUTATIONRATE, --thetamutationrate MUTATIONRATE : This is the mutation rate per haplotype as used in the coalescent model. The resulting mutation rate ```lambda``` at a particular site is given by ```lambda = theta/(2(N + theta))``` where ```N``` is the number of reference haplotypes. The option exists to omit this and directly pass the value for ```lambda```.
6. -L LAMBDAMUTATIONRATE, --lambdamutationrate LAMBDAMUTATIONRATE : This is the ```lambda``` (direct error) Mutation rate. If the user intends to associate a general error rate due to genotyping or the possibility of *de novo* mutation at any site without reference to the coalescent model, this is the parameter to set.
7. Note about mutation rates: In the absence of both the mutation rate parameters and when the ```posspecific``` parameter (see below) is set to ```False```, the mutation rates are set as follows: ```thetamutrate = (sum(1.0/i for i in range(1,2*refpop-1)))**(-1); lambdamutrate = 0.5*mutrate/(mutrate + 2*refpop)```. If ```posspecific=True```, the above mutation rate parameters are ignored, and the position-specific mutation rates in the observed sample file are used.
9. -e EFFPOP, --effpop EFFPOP : Effective size of the human population. Default value = 11,418 [[3]](#3).
10. -r REFPOP, --refpop REFPOP : Size of the reference population, i.e. the number of genotypes in the reference database (not the number of haplotypes).
11. -s {distance,custom}, --recombswitch {distance,custom} : Choose whether to use a linear growth in recombination rate with genomic distance, or a custom file of recombination values. If ```custom``` is chosen, provide a file for the recombination rate between adjacent SNPs (for L SNPs, there will be L-1 such recombination values).
12. -b RECOMBRATE, --recombrate RECOMBRATE : If the ```distance``` option is chosen for the ```-s``` flag, then provide the recombination rate in ```cM/Mb```, i.e. centimorgans/Megabase. In the paper, we mainly used a value of 0.5 cM/Mb. This is set as the default value.
13. -f RECOMBFILE, --recombfile RECOMBFILE : If the ```custom``` option is chosen for the ```-s``` flag, provide the name of the file here. In the ```distance``` option, each value is set to ```4 * Effective Population Size * distance in Mb * Recombination rate in cM/Mb```. Thus, in the ```custom``` case, the user should set each recombination values between loci L and L+1 to ```4 * Effective Population Size * distance in cM between loci L and L+1```.
14. -t TOLERANCE, --tolerance TOLERANCE : This is the tolerance factor that determines the cutoff for the sub-optimal paths to be included. That is, if the score of the best-fit trajectory is ```S```, all paths with a score ```S - TOLERANCE*S``` will be included. For the analyses in the paper, we chose a tolerance of 0.01. This is set as the default value.
15. -C CURRDIR, --currdir CURRDIR : The current working directory to run the code. The default is set to './'.
16. --numproc NUMPROC : Number of processes for parallelization. The more the merrier. The default is set to 1.
17. --posspecific POSSPECIFIC : Are position-specific mutation rates included in observation file? (True/False). The default is 'False'.
18. --prefix PREFIX : String prefix to append to output Best_trajectories file, in addition to chromosome number. The default is '' (empty string).
19. --truncation TRUNCATION (**Only in PLIGHT_Truncated**): Fraction of trajectories carried forward from each step, after the phase-in period (see Supplementary Methods of the **PLIGHT** paper.
20. --subgroup SUBGROUP (**Only in PLIGHT_Iterative**): Number of haplotypes in each iterative scheme subgroup, where the subgroups are defined as the partitions into which the overall reference haplotype set is divided for the **PLIGHT_Iterative** run.
21. --niter NITER (**Only in PLIGHT_Iterative**): Number of iterations of bootstrapping process. The default is 20, though values of 30 were also considered in the paper.

## Example run
An example of one run of the code for the 1000 Genomes Phase 3 reference database is:
```
python3 PLIGHT_Iterative.py --metadata integrated_call_samples_v3.20130502.ALL.panel --genfolder Genotypes/ --lambdamutationrate 0.1 --effpop 11418 --refpop 2504 --recombrate 0.5 --recombswitch distance --tolerance 0.01 --currdir ./ --numproc 20 --subgroup 300 --niter 1 --posspecific False -c ALL.chr3.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz -O chr3_Observed_SNPs.txt --chromosomeID chr3 --prefix Mutrate0.1_Subgroup300
```
## PLIGHT_Vis
The **PLIGHT_Vis** module carries out downstream analyses on the output of the HMM algorithms. This includes plotting the trajectories, and calculating the optimal haplotypes for each chromosome, as well as across all chromosomes under consideration. The help menu is as follows:

```
usage: PLIGHT_Vis.py [-h] [-C CHROMOSOMEIDS [CHROMOSOMEIDS ...]]
                     [-t TRAJECTORYPATTERN] [-T TRAJECTORYFOLDER]
                     [-P PLOTFOLDER]

Unravel the trajectories, find the identities of the best-fit reference
haplotypes, and plot the results

optional arguments:
  -h, --help            show this help message and exit
  -C CHROMOSOMEIDS [CHROMOSOMEIDS ...], --chromosomeIDs CHROMOSOMEIDS [CHROMOSOMEIDS ...]
                        List of Chromosome IDs to be considered: format chr
                        followed by number, listed singly, separated by spaces
                        (for eg. -C chr1 chr2 chr19)
  -t TRAJECTORYPATTERN, --trajectorypattern TRAJECTORYPATTERN
                        Pattern for trajectory files, where the chromosome ID
                        would fit into the curly braces
  -T TRAJECTORYFOLDER, --trajectoryfolder TRAJECTORYFOLDER
                        Folder for storing the trajectory files
  -P PLOTFOLDER, --plotfolder PLOTFOLDER
                        Folder for storing the plots
```

An example of running **PLIGHT_Vis** is provided below:
```
python3 PLIGHT_Vis.py -C chr3 chr6 -t InRef_{}_Best_trajectories.tsv -T Trajectories -P HMM_Plots
```
## Example data included
### Inputs
We include examples of input SNP lists for both cases of (a) non-position-specific and (b) position-specific mutation rates. See the parameter descriptions above for explanations on how these two input files differ.

### Outputs
We provide examples of the best-fit trajectories files produced by the algorithms, where each line represents the pointer from the best-fit pair of haplotypes at the previous query SNP to the corresponding best-fit pair(s) at the current query SNP. Multiple possible best-fit pairs of haplotypes at the current query SNP are separated by semi-colons. The query SNPs are shown in reverse order as this is the manner in which the Viterbi algorithm identifies the optimal trajectories.

## References
<a id="1">[1]</a>
Howie, B. N.; Donnelly, P.; Marchini, J. A Flexible and Accurate Genotype Imputation Method for the next Generation of Genome-Wide Association Studies. PLoS Genet. 2009, 5 (6). https://doi.org/10.1371/journal.pgen.1000529.

<a id="2">[2]</a>
Loh, P. R.; Danecek, P.; Palamara, P. F.; Fuchsberger, C.; Reshef, Y. A.; Finucane, H. K.; Schoenherr, S.; Forer, L.; McCarthy, S.; Abecasis, G. R.; et al. Reference-Based Phasing Using the Haplotype Reference Consortium Panel. *Nat. Genet.* **2016**, 48 (11), 1443–1448. https://doi.org/10.1038/ng.3679.

<a id="3">[3]</a>
Li, N.; Stephens, M. Modeling Linkage Disequilibrium and Identifying Recombination Hotspots Using Single-Nucleotide Polymorphism Data. *Genetics* **2003**, 165, 2213–2233.

<a id="4">[4]</a>
Viterbi, A. Error Bounds for Convolutional Codes and an Asymptotically Optimum Decoding Algorithm. *IEEE Trans. Inf. Theory* **1967**, 13 (2), 260–269. https://doi.org/10.1109/TIT.1967.1054010
