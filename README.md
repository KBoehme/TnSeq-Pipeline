
#TnSeq Pipeline


##Overview

This script is designed to take Tn-Seq data and map (using Bowtie2) and tally hop sites on a user-provided reference genome.

###Input
1. Fasta/Fastq read files.
2. Bowtie2 reference genome.
3. Genbank Ptt files for reference genome.

###Output

1. Gene-by-gene (-GENE). Hop counts are tallied within genes.
```
Num	GeneID	Condition1	Condition2	Condition3	Start	Stop	Strand	Length	PID	Gene	Function
...
12	SM_b20013	1843	282	2289	19079	20323	-	996	16263764	-	hypothetical protein
13	SM_b20014	1316	258	1176	20377	21378	-	801	16263765	-	transcriptional regulator
14	SM_b20015	1924	174	1770	21392	22348	-	766	16263766	-	sugar ABC transporter permease
15	SM_b20016	4263	6649	4907	22345	23337	-	794	16263767	-	sugar ABC transporter permease
16	SM_b20017	4580	5036	4044	23337	24842	-	1205	16263768	-	sugar ABC transporter ATP-binding protein
17	SM_b20018	1485	1632	2146	25092	26114	-	818	16263769	-	sugar ABC transporter substrate-binding protein
18	SM_b20019	2498	635	4363	26206	27342	-	910	16263770	sucB	dihydrolipoamide succinyltransferase
19	SM_b20020	3636	976	5119	27345	29423	-	1664	16263771	pdh	pyruvate dehydrogenase E1 component, subunits alpha and beta
...
```

2. Hop-by-hop (-HOPS). Breaks down gene-by-gene ouput, showing each individual hop site observed within that gene.
```
Num	GeneID	Condition1	Condition2	Condition3	Start	Stop	Strand	Length	PID	Gene	Function
...
12	SM_b20013	2660	3455	317	19079	20323	-	996	16263764	-	hypothetical protein
	SM_b20013	39	5	0	19206		-
	SM_b20013	94	249	0	19207		-
	SM_b20013	110	141	0	19222		+
	SM_b20013	12	0	0	19225		-
	...
	...
	SM_b20013	0	106	0	20113		-
	SM_b20013	0	0	16	20120		-
	SM_b20013	193	84	8	20146		+
	SM_b20013	135	52	0	20152		+
13	SM_b20014	1878	1770	303	20377	21378	-	801	16263765	-	transcriptional regulator
	SM_b20014	61	0	0	20480		-
	SM_b20014	51	5	6	20503		-
	SM_b20014	0	143	0	20521		-
...
```


##Getting Started

### 1. Download Repo

Option 1 - Run the following command in your terminal

```
git clone https://github.com/KBoehme/TnSeq-Pipeline.git
```

Option 2 - Download this repo as a zip file (Find button above and on the right hand side).


### 2. Get Bowtie2

See the Bowtie2 [homepage](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) and [download link](http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.4/) for information and instructions on how to install Bowtie2.

This script relies on Bowtie2 in order to map reads. Make sure it is installed and the executable is in the PATH. You should be able to run it simply by typing:
```
bowtie2
```

### 3. Prepare Reference Genome

##### Download organism's reference genome from Genbank

Download GFF version of the genome [here](http://www.ncbi.nlm.nih.gov/guide/howto/dwn-genome/).

##### Indexing a Reference Genome
 Bowtie2 requires an indexed reference genome to run and as such this pipeline requires that this indexed genome be created prior to running. Using the genbank files you just downloaded see [here](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#getting-started-with-bowtie-2-lambda-phage-example)
 under the section header "Indexing a reference genome" for instructions on how to create an indexed reference genome.


### 4. Create config file

The user must create a configuration file which tells the script where to find input files as well as specifying to the program other essential paramenters. Below is the example.config file that can be found in the example_data directory within this repository.

```
[input]

Reads           = ./example_data/test_data/*.fastq             ; Absolute or relative path (from TnSeq-Pipeline.py script) to TnSeq reads.
BowtieReference = ./example_data/bowtie2_reference/smeliloti   ; Path to bowtie2 reference genome. Use prefix.
Ptt             = ./example_data/bowtie2_reference/ptt/*.ptt   ; Path to PTT files of reference genome.
Out             = My_example_run                               ; Name for the output. (For example: My_example_run-GENE.txt and My_example_run-HOPS.txt)


################################

[parameters]

Transposon      = TCGAGATGTGTATAAGAGACAG   ; Transposon sequence
Mismatches      = 3                        ; Mismatches allowed when finding transposon
GeneTrim        = 10                       ; Percent of gene length truncated on both sides of gene where hops wont be counted.
ReadLength      = 25                       ; If read is shorter than this length after removing transposon it will be removed, otherwise it will be trimmed to this length and mapped.
MinimumHopCount = 1                        ; Any hop site with less than this number of occurrences (totaled across all conditions) will be removed from analysis and will not be outputed in results.

################################

[options]

Debug                   = False           ; Shows detailed running parameters for debugging purposes.
Normalize               = Intergenic      ; Options: [Intergenic, Total] Intergenic: Normalize based on intergenic hops only. Total: Normalize based on all hops.
DeleteIntermediateFiles = True            ; Delete intermediate fasta and sam files at the end of the run.
ReverseComplementReads  = True            ; If True this will take the reverse complement of all reads in the input fasta/fastq file before searching for the transposon sequence.
```




