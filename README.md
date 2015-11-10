MetaFinisherSC
=============

## Introduction ##
This tool is for users to upgrade their metagenomics assemblies using long reads. This includes fixing mis-assemblies and scaffolding/gap-filling. If you encounter any issues, please contact me at kklam@eecs.berkeley.edu. My name is Ka-Kit Lam. 


## Commands to run the tool ##

Assume you have your long reads (LR.fasta) and contigs formed from short/long reads (LC.fasta). 

The final output is **abun.fasta** after running through the M-Fixer(mis-assembly fixing) and A-Splitter(join contigs based on abundance and overlap information). 

In order to do that, here are the steps. 

1. Run M-Fixer to fix misassemblies with information from reads (input is LR.fasta, LC.fasta ; output is mFixed.fasta).

        python -m srcRefactor.misassemblyFixerLib.mFixer destinedFolder mummerPath 

2. Run A-Splitter  to join contigs based on abundance and overlap information (input is mFixed.fasta, raw\_reads.fasta; output is abun.fasta). Note that raw\_reads.fasta is copied from LR.fasta after running M-Fixer, so there is no need to rename if you run M-Fixer first. But if you want to use A-Splitter alone, then the input filenames are as stated. 

        python -m srcRefactor.repeatPhaserLib.aSplitter destinedFolder mummerPath

	
Sometimes, if the names of raw reads and contigs consists of special characters/formats, MetaFinisherSC/MUMmer may not parse them correctly. In that case, you want to have a quick renaming of the names of contigs/reads in LC.fasta or LR.fasta using the following command. 

        perl -pe 's/>[^\$]*$/">Seg" . ++$n ."\n"/ge' LR.fasta > newLR.fasta
        cp newLR.fasta LR.fasta
        perl -pe 's/>[^\$]*$/">Seg" . ++$n ."\n"/ge' LC.fasta > newLC.fasta
        cp newLC.fasta LC.fasta

## Dependency ##
MetaFinisherSC requires the following libraries and softwares. Recommended way to download the packages is as follows.

1. Download the python libraries. 

        sudo pip install scipy numpy biopython python-Levenshtein

  Remark: For Linux users, you may also replace pip install step by "sudo apt-get install -qq python-numpy python-scipy python-biopython python-Levenshtein"
  
2. Download MUMmer3.23 

        wget http://sourceforge.net/projects/mummer/files/mummer/3.23/MUMmer3.23.tar.gz/download -O /tmp/mummer.tar.gz
        tar -zxvf /tmp/mummer.tar.gz -C /tmp/
        make install -C /tmp/MUMmer3.23/
       
3. Download ClustalW-2.1 (optional)

        wget http://www.clustal.org/download/current/clustalw-2.1.tar.gz -O /tmp/clustalw.tar.gz
        tar -zxvf /tmp/clustalw.tar.gz -C /tmp/
        cd /tmp/clustalw-2.1/
        ./configure
        make
        export PATH=$PATH:/tmp/clustalw-2.1/src/

## Example ##
Below is a step by step example on running MetaFinisherSC on the testset provided. In this example, there are two misassembled contigs and we will fix the misassemblies and join them back correctly. The reads are synthetic reads extracted from two synthetic species of different abundances(20X and 50X respectively). Both species share a common segment of length 12000 bp and that the readlength is 6000bp.  

1. Clone MetaFinisherSC
        
        git clone https://github.com/kakitone/MetaFinisherSC.git
        
2. Get into the MetaFinisherSC folder, create directory to hold the data, download data  and check to see if you get everything(LC.fasta, LR.fasta, reference.fasta). 
        
        cd MetaFinisherSC/        

        mkdir dataFolder
        
        cd dataFolder
        
        wget -O testData.zip http://www.eecs.berkeley.edu/~kakitone/testData.zip 
        
        unzip testData.zip
        
        ls -lt

3. Go back to the main directory containing MetaFinisherSC and run the following two commands. Here we assume that /tmp/MUMmer3.23/ is the path to your MUMmer. This step can take a few minutes. 
        
        cd ../

        python -m srcRefactor.misassemblyFixerLib.mFixer dataFolder/ /tmp/MUMmer3.23/
        
        python -m srcRefactor.repeatPhaserLib.aSplitter  dataFolder/ /tmp/MUMmer3.23/

4. Now verify that the original LC.fasta contain 2 contigs whereas the final abun.fasta contains 2 contigs

        fgrep -o ">" dataFolder/LC.fasta  | wc -lc
        
        fgrep -o ">" dataFolder/abun.fasta  | wc -lc


5. Align the abun.fasta with reference.fasta. 

        /tmp/MUMmer3.23/nucmer  -maxmatch dataFolder/abun.fasta dataFolder/reference.fasta         

        /tmp/MUMmer3.23/show-coords out.delta

        |   [S1]     [E1]  |     [S2]     [E2]  |  [LEN 1]  [LEN 2]  |  [% IDY]  | [TAGS] 
        |=====================================================================================
        | 2500006  2512005  |  2500001  2512000  |    12000    12000  |   100.00  | Segkk1	Segkk0
        |       1  5000002  |  5000000        1  |  5000002  5000000  |    99.99  | Segkk0	Segkk0
        |       1  5000005  |        1  5000000  |  5000005  5000000  |    99.99  | Segkk1	Segkk1
        | 2488001  2500002  |  2512000  2500001  |    12002    12000  |    99.83  | Segkk0	Segkk1

6. As a check, you may also want to see that the original contigs are misassembled. This can be seen by 

        /tmp/MUMmer3.23/nucmer  -maxmatch dataFolder/LC.fasta dataFolder/reference.fasta         

        /tmp/MUMmer3.23/show-coords out.delta

        |    [S1]     [E1]  |     [S2]     [E2]  |  [LEN 1]  [LEN 2]  |  [% IDY]  | [TAGS]
        |=====================================================================================
        |       1  2512000  |        1  2512000  |  2512000  2512000  |   100.00  | Segkk0	Segkk0
        | 2500001  5000000  |  2500001  5000000  |  2500000  2500000  |   100.00  | Segkk1	Segkk0
        |       1  2512000  |        1  2512000  |  2512000  2512000  |   100.00  | Segkk1	Segkk1
        | 2500001  5000000  |  2500001  5000000  |  2500000  2500000  |   100.00  | Segkk0	Segkk1

7. Also, you may also want to check that there is really a common segment of length 12000 across species

        /tmp/MUMmer3.23/nucmer  -maxmatch dataFolder/reference.fasta dataFolder/reference.fasta         

        /tmp/MUMmer3.23/show-coords out.delta
        
        |    [S1]     [E1]  |     [S2]     [E2]  |  [LEN 1]  [LEN 2]  |  [% IDY]  | [TAGS]
        |=====================================================================================
        |       1  5000000  |        1  5000000  |  5000000  5000000  |   100.00  | Segkk0	Segkk0
        | 2500001  2512000  |  2500001  2512000  |    12000    12000  |   100.00  | Segkk1	Segkk0
        |       1  5000000  |        1  5000000  |  5000000  5000000  |   100.00  | Segkk1	Segkk1
        | 2500001  2512000  |  2500001  2512000  |    12000    12000  |   100.00  | Segkk0	Segkk1
     
     
##Current build status##
It is an indicator on the current built powered by Travis-CI. If you issue a pull request, Travis-CI will evaluate your suggestion by automatically running the code on the default test case. 

![alt text](https://travis-ci.org/kakitone/MetaFinisherSC.svg?branch=master "Current build status")

Also, if you want to have a look at the commands and stdout of the sample run, you can visit https://travis-ci.org/kakitone/MetaFinisherSC

