# FindEnzymes

## Description of scripts

- extract_SNP.jl: This script extracts the SNP/Indel/Ins/Del/Dup from a cosmic mutation file.
- modify_rawDNA.jl: This script modifies the raw DNA sequence to a new DNA sequence using mutation information.
- process_cond*.jl: This script find the mutations that are under different conditions.
- filter_step*.jl: This script filters the mutations that are under different conditions.
- check_cond4.jl: This script checks the results of condition 4.

## How to use the scripts

1. Install Julia from https://julialang.org/downloads/
2. Install the required packages by running the following command in the terminal:
```bash
add BioSequences
add FASTX
```
3. Run the scripts by running the following command in the terminal:
3.1 Extract the SNP/Indel/Ins/Del/Dup from a cosmic mutation file:
```bash
julia extract_SNP.jl
```
3.2 Made a bedfile from the mutation file and then extract the DNA sequence from the bedfile using bedtools
3.3 Made a new DNA sequence using the mutation information:
```bash
julia modify_rawDNA.jl
```
3.4 Find the mutations that are under different conditions:
```bash
julia process_cond1.jl
julia process_cond2.jl
julia process_cond3.jl
julia process_cond4_1.jl
julia process_cond4_2.jl
julia process_cond5.jl
```
3.5 Filter the mutations that are under different conditions:
```bash
julia filter_step1.jl
julia filter_step2.jl
julia filter_step3.jl
julia filter_step4.jl
julia filter_step5.jl
```
3.6 Check the results of condition 4:
for example:

```bash
python3 check_cond4.py chr22:41140144 REs.tsv AspS9I "p41140140|CGNCC|0" flIII
```

```
usage: check_cond4.py [-h] coordinate enzyme name meta matched

Fetch a sequence around a given position with the position highlighted.

positional arguments:
  coordinate  Chromosome coordinate in a format like: chr2:123456
  enzyme      The file of the enzyme
  name        The name of the enzyme
  meta        meta info
  matched     The name of the enzyme

options:
  -h, --help  show this help message and exit
```


