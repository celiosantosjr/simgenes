# Simgenes 

Pipeline to simulate and create random genes with selected features.

[![Simgenes - Python package](https://github.com/celiosantosjr/simgenes/actions/workflows/python-package.yml/badge.svg)](https://github.com/celiosantosjr/simgenes/actions/workflows/python-package.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

If you use this software in a publication please cite:

>   Santos-Júnior C.D., 2023.
>   Simgenes: Random gene sequences simulator.
>   GitHub Repository. Available in: https://github.com/celiosantosjr/simgenes/

---
## Table of contents

1. [Introduction](#introduction)
2. [Installing](#install)
3. [Usage](#usage)
4. [Contact](#contact)

---
## Introduction

Simulate random gene and protein sequences controlling for real parameters, such as GC content.
The package simgenes rely on very common libraries:

   - Bio
   - tqdm
   - numpy

In a test with 1 million sequences generated using simgenes:

```
from simgenes import batch_random

batch_random(ofile='test_check', L=(30,1002), GC=None, transl=False, n=1_000_000)
```

The time to generate 1 million gene sequences was 2765 s, equivalent to say that to simgenes can
generate on average 361.66 sequences per second.

We obtained a file with gene sequences in a length ranging from 30 to 1002 bp,
with a very even distribution:

<p align="center">
  <img src="https://github.com/celiosantosjr/simgenes/blob/main/.assets/random_gene_length_distribution.svg" alt="drawing" width="450"/>
</p>
<p align="center">
  <em><b>Fig. 1</b> - Length distribution of generated genes in a simulation.</em>
</p>

Also, the contents of the nucleotides in the genes ranged largely, but yet inside the limits established 
internally in the pipeline according to [other authors](https://pubmed.ncbi.nlm.nih.gov/23028785/).

<p align="center">
  <img src="https://github.com/celiosantosjr/simgenes/blob/main/.assets/GC_distribution.svg" alt="drawing" width="450"/>
</p>
<p align="center">
  <em><b>Fig. 2</b> - GC contents distribution in genes from the simulation.</em>
</p>

The 1 million genes generated did not match to any family of Uniprot, and did not cluster under
cdhit at 95% of identity and 90% of coverage of the shorter sequence. This shows that random genes
from simgenes definitely differ from the real genes, being interesting to work as a filter of
artifacts in gene prediction. This system can easily be used to generate decoy sequences for
multiple purposes.

---
## Install

To install the package simgenes, you can easily create a conda environment:

```
conda env create -f simgenes_environment.yml
conda activate simgenes
```

Then, just install the simgenes:

```
python3 -m pip install simgenes
```

---
# Usage

Simgenes has several functions to help the user to quantify randomness in the sequences or
generate new sequences.

**1. Entropy measurements**

This function can be used by simply entering the protein or nucleotide sequence
and returns the entropy measured, the maximum possible entropy, and the relative
entropy measure - Note that for all entropy calculations, we rely on the Log2 
approximating the measures to the Shannon Entropy. Example of usage:

```
from simgenes import entropy

# creating example sequences
seq = 'ATACAGATGCGAGTAAACGAGCAAAAA'
seq_aa = 'MKLMKLMLKLMLLLLLLAHHHHHHH'

H, maxH, Hrel = entropy(seq)

print('''For nucleotides
We measured an entropy of: {H},
a possible maximum entropy of: {maxH},
and as relative entropy of {Hrel}''')

H, maxH, Hrel = entropy(seq_aa)

print('''For amino acids
We measured an entropy of: {H},
a possible maximum entropy of: {maxH},
and as relative entropy of {Hrel}''')
```

**2. Translate nucleotides into amino acid sequence**

This function can be used by simply entering the nucleotide sequence
and returns the translated sequence by using Standard Genetic Code, 
translation table 11. Example of usage:

```
from simgenes import get_prot

# creating example sequence
seq = 'ATACAGATGCGAGTAAACGAGCAAAAA'

protein = get_prot(seq)

print(protein)
# IQMRVNEQK
```

**3. Check internal stop codons**

This function returns the indexes of nucleotides from which
an internal stop codon is formed in the first frame. For example
in the sequence ATGTAAAAATTTTGA, this function would return a list
with the indexes of stop codons in the sequence [1, 4]. With this
information, one can easily (by multiplying by 3) check the position
of such codons in the sequence: [3 - TAA, 12 - TGA]. Example of
usage:

```
from simgenes import check_internal_stops

# get example sequence
seq = 'ATGAGGGAGGGGGTGGGGATGAGGGATGAAAGTGTAGGATAAGTGACGGGGTGTAAGGCTGGCTTGGTGGACC
TGCAATTGAATGGCAAGGGCGGAGTAGATGTGGTGAGTGGAGAGAAGGTGGGAAAGTGTCAAATTCATGTGGGACGGTGGT
GGTTGGAGGGTTGCGGGAAGATCGGTGTGTGGAAGTTAGTGGGTTAAAGGGAACGTGCAGTAGGGGGAGCAGCTGGACAGC
TGGAAGGCGGGTCTGAGCGGCGGGGTGGGGTGGGTGGGGCGGGTTGGGGGGTTGGTGGTGTGCGGGAAAATCGGTGGGGGG
TATGGGTTAGGGCGAGGGGGGGGTTCAGAAGTGTGGGTGGGCGAGGCGAATGAGGAAAGGAGGGGGAGGTGGGAATGTGTG
GGGGACATGCGGGTAGGGAGAGCGAGCGGGAGTTGGTGGGGGGGAGTGGGGTACGGGGTTGGATGAAGGAGTGATGGGAGC
GGTAGAGAAAGGCGCTGGAGGGAAGCGAAGGGTTAGTAGCGAAGCAAGGTGGTTTGGGCGACGAAATTCGGGTGGTAAGTT
ATAGTGGTATAGGCAGCAGGTAA'

stop_codons = check_internal_stops(seq)

# getting stops in the sequence:
for x in stop_codons:
   print(seq[x*3: (3+(x*3))])

# 'TAA'
# 'TAA'
# 'TGA'
# 'TGA'
# 'TAG'
# 'TAA'

# getting indexes of the stop codons:
for x in stop_codons:
    print((x*3, 3+(x*3)))

# (39, 42)
# (198, 201)
# (366, 369)
# (468, 471)
# (480, 483)
# (579, 582)
```

**4. Determine random nucleotides distribution**

Although users can easily attribute a GC content for the generated random sequence,
the exact distribution of A, T, C, G is not allowed to be set. This is because we are
interested on generating the most random gene sequences possibly, therefore, the
distribution is mostly attributed by randomly selecting the proportions that add up to
100%, also respecting the GC proportion previously established. When user gives
GC = -1, then the proportions are set to 25% for each nucleotide, alternatively user
can attribute 0 or no value (None) to get a completely random distribution that respects
the rule from [Mir et al. (2012)](https://pubmed.ncbi.nlm.nih.gov/23028785/) that says
bacterial ORFs present a GC content between 21.4% and 74.9% of GC. For example:

```
from simgenes import det_props

# using fixed proportions
# [a = 25.00%, t = 25.00%, c = 25.00%, g = 25.00%]
props = det_props(GC=-1)

# using completely random assignment
props = det_props(GC=None)
# alternatively:
props = det_props()
# as result, for example:
# [23.2, 22.5, 18.8, 35.5]
# in this case GC was randomly set to 54.3%

# using an established GC content
props = det_props(GC=66.7)
# as a result, we can have, for example:
# [a = 22.60%, t = 10.70%, c = 27.60%, g = 39.10%]
```

**5. Getting a random length**

Unless user specifies clearly one single length for the generated
sequence (multiple of 3 because of codons counting), simgenes will attribute
a random length to the gene. It is also possible to create a range for the
length of this gene in a form of a tuple, e.g. (30, 300) in this case meaning
that the gene should encode proteins of 10 to 100 codons in total (counting start 
and stop as well). If a length that is not divisible by 3 is selected by the 
user, then the program overules the input and select a random length. Usage
example follows:

```
from simgenes import get_len

# getting a random length between 30 and 300
L = get_len(30, 300)
```

In simgenes, the length of genes by default varies between 300 and 900 bp,
with an average of 819 bp as suggested by
[Rajic et al. (2005)](https://pubmed.ncbi.nlm.nih.gov/15883855/).
It is also important to mention that it is possible to generate small
genes as well by setting low multiple of 3 numbers, such as 30, 120, and so on.

**6. Random gene**

This function returns a gene (and) protein sequence. User can provide a range to the
length of random genes, by attributing to L a tuple with min and max gene lengths in
base pairs, respectively. Users also can provide a fixed gene length by using an
integer, or yet conform with the randomly selected length from the range between
300 to 900 base pairs. Another parameter to be selected is the GC gc content (0-100)
of the desired genes, in case it is not given, the program will select randomly the
proportion of nucleotides between 21.4 <= GC <= 74.9 (typical GC content of ORFs). 
There is also the option GC=-1, which means that all nucleotides are expected at same
frequency. If user wants the corresponding protein sequence for the generated
gene, then just set the "trans" parameter as True. Examples of usage:

```
from simgenes import random_gene

L = [(30, 300), 333, 40, None]
gc = [-1, 75.2, None]
for m in L:
    for g in gc:
        s, p = random_gene(L=m, gc=g, trans=True)
        print(f'L={m}bp, GC={g}%, nt={s}, prot={p}')
```
    
**7. Batch random - To generate a set of genes**

The option "batch_random" is interesting when the users need to generate a set of random
genes. It takes the same parameters of the random_gene function, such as length in bp (L),
GC content (GC=0-100), translation (transl=True/False), and additionally accepts the
outputfile address tag (as a string) and the number of desired sequences (n=int). If no
parameter is given, this function should generate 1 million nucleotide sequences from 30 to
300 bp with GC=50.8%, not translating them, and saving them into a file named
'output_testseq.fna.xz'. Example of usage:

```
from simgenes import batch_random

batch_random(ofile='output_testseq',
             L=1500,
             GC=65.2,
             transl=False,
             n=1_000)
```

---
## Contact

If you have any issues or want to contribute to the project, please contact the author [Celio Dias Santos Junior](mailto:celio.diasjunior@gmail.com).
