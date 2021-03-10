# AptaDesign
AptaDesign is a python-based software tool for in silico aptamer design based on motif extraction and directed evolution. This software tool helps with the creation of aptamers targeting a specific molecule (based on previous data available) with the possibility of also hybridizing a specific oligonucleotide sequence for their use in genetic circuits.
This algorithm requires a fasta or fastq file containing previously experimentally-demonstrated aptamer sequences targeting the desired molecule or a similar molecule. From this file a given number of conserved motifs is extracted and used for evaluation of the directed evolution procedure.
In the directed evolution part of the algorithm, for each generation random mutations are generated, which are then evaluated using the Edit Distance with the conserved motifs, the MFE of the secondary structure and the hybridation MFE (if required). The evolution continues until any of the stop conditions are reached.

If you use this software for a publication, please cite:

Rueda-Silva, Juan Carlos. (2021). AptaDesign

## Tutorial
### Installaton
#### 1.1. Contents of the download bundle
To begin please download the compiled distribution. The distribution contains the following files, which will be required during the installation and for the succeful run of the program. The files are:
- AptaDesign.exe and build components: The compiled software tool and external utility files required for te execution of the program.
- Examples folder: a series of premade files that will be used for testing the installation
- Fc_Aptamers.fasta: A fasta file containing the compilation of aptamers previously reported in literature specific for the Fc stem of mammalian IgG. This file contains the aptamers reported by Yang et al. (2020), Miyakawa et al. (2008), Bognár & Gyurcsányi (2020), Yoshida et al. (2019) and Ma et al. (2013).
- HRP_DNAzyne.txt: The seuqnece of a HRP-mimicking DNAzyme, previously reported by Alizadeh et al. (2020).
- ViennaRNA_Package: A partial installation of the ViennaRNA Package conatining the minumum files required for the use of the compiled version of RNAfold, which it is used by this software for MFE calculation and secondary structure predeiction.
- RNAhybrid-2.1.2: A partial installation of the compiled RNAhybrid software tool. It is used by this software tool for hybridation MFE calculation.
- LICENSE: License information for this sofwtare tool. This software tool is licensed under the MIT license.
- README.md: This document.

#### 1.2. Installing the software
Extract the entire contents of the downloaded zip folder, it should take around 15 minutes.

After this the installation should be ready.


#### 1.3. Testing the installation
In orther to test that the installation was succefull, in a command line within the folder conatining the downloaded files, please run the following command:

```
> AptaDesign --help
```

Which should generate an output starting with:

```
Usage: AptaDesign [OPTIONS]

AptaDesign v.1.0.4 Aptamer design in silico based on directed evolution and previous data. With possibility for creating aptamers capable of hybridizing specific sequences.

USAGE: AptaDesign.exe [OPTIONS]* *Note: if no option are provided, the programs starts in its UI mode
...
```

Also test that the utilities are succefully installed by typing the following command:

```
> AptaDesign -s .\examples\Fc_Aptamers.fasta -g 2 -t GTGGGGCATTGTGGGTGGGTGTGG -ph <Your desired output path> -n test.1
Repleace <Your desired output path> with the path to the directory where you want to save the output.
```


Which should output after a few seconds something similar to:

```
Welcome to AptaDesign v1.0.4
------------------------------------------
```

And after around 30 minutes:

```
Output located in: <Your path>\\test.1
------------------------------------------
Thank you for using AptaDesign!
```

And in your output directory you should find the following 3 files with contents similar as the following:

- test.1.log
- test.1_results.xlsx
- test.1_fig.png

If these tests shows any errors please go through the installation again.


### Usage
#### 2.1. Using the interface mode
In order to start the software in the User Interface mode, type the following command:

```
> AptaDesign
```

And follow the instructions as they are promted the terminal. It will ask for the input fasta or fastq file, number of motifs to use for the evaluation, pool size (number of sequences generated per cycle) to use, number of generations (cycles) to simulate, number of candidates (number of sequences selected per cycle), number of generations without score increase to stop the program, hybridation target (if desired), the path to save the output, title of the project and wether to generate graphs. If you enter to Advanced Options you can enter the minumum length of the aptamers, the number of generations per hyperdiverse period, the maximum number of consecutive hyperdiverse periods before stopping the program and the minimum score that will cause the program to stop.


#### 2.2. Using the command line tool mode
In order to start the software in the function mode, type the following command:

```
> AptaDesign [OPTIONS]
```

With the options shown bellow.

Note: For the options not marked as required, if no value is provided, the default value will be used.

Options for the function mode
- -s, --source :(required) The path to the fasta or fastq file containing the sequences from which the algorithm will extract the motifs.
- -m, --motifs : The number of motifs to be used for the evaluation. (default: 8)
- -t, --target :(optional) The sequence of a hybridation oligonucleotide target if desired.
- -p, --pool : Number of sequences to generate per cycle. (default: 100)
- -g, --gens : Number of cycles to simulate. (default: 1000)
- -c, --candidates : Number of sequences to select after each cycle. (default: 10)
- -l, --length : Minimum length of the apatmers. (default: 10)
- -hg, --hyperd_gens: Number of subgenerations per hyperdiverse period. (default: 3)
- -ch, --cons_hyperd: Maximum number of consicutive hyperdiverse periods. (default: 10)
- -cs, --cons_score: Maximum number of consecutive generations without score improvement before stopping the program. (default: 10)
- -gp, --graphs : ('T'/'F') Generate graphs for the project. (default: 'True')
- -bc, --break_score : Minimum score that will cause the program to stop. (default: 0.99)
- -ph, --path :(required) Path to save the output of the project.
- -n, --name :(required) Name of the project.
- --help : Open the help page.

### 2.3. Output
After the succeful run of the algorithm, a folder with the name of the project is created on the given output directory. This folder contains the following files

- A log file: Containing messages generated during the execution with timestamps.
- A png image: A graph of each generation vs. maximum score acgieved.
- A xlsx file: Containing the final selected candidates, their score and the values used for the calculation of the score

Please read this software's documentation for more information: https://juancrueda.github.io/AptaDesign/documentation/index.html
