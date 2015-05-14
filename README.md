# MetaFinisherSC
### Introduction ###
This tool is for users to upgrade their metagenomics assemblies using long reads. This includes fixing mis-assemblies and scaffolding/gap-filling. If you encounter any issues, please contact me at kklam@eecs.berkeley.edu. My name is Ka-Kit Lam. 


### Only long reads and contigs formed from long reads ###

Assume you have your long reads (LR.fasta) and contigs formed from long reads (LC.fasta). 

The final output is **abun.fasta** after running through the M-Fixer(mis-assembly fixing) and A-Splitter(join contigs based on abundance and overlap information). 

In order to do that, here are the steps. 

1. Run M-Fixer which fix misassemblies with information from reads(key intermediate output here is mFixed.fasta).

        python -m srcRefactor.misassemblyFixerLib.mFixer destinedFolder mummerPath 

2. Run A-Splitter to join contigs based on abundance information (key intermediate output here is abun.fasta)

        python -m srcRefactor.repeatPhaserLib.aSplitter destinedFolder mummerPath

### Example ###
Below is a step by step example on running MetaFinisherSC on the testset provided. In this example, there are two misassembled contigs and we will fix the misassemblies and join them back correctly. The reads are synthetic reads extracted from two synthetic species of different abundances(20X and 50X respectively). Both species share a common segment of length 12000 bp and that the readlength is 6000bp.  

1. Clone MetaFinisherSC
        
        git clone https://github.com/kakitone/MetaFinisherSC.git
        
2. Get into the MetaFinisherSC folder, create directory to hold the data, download data from Dropbox and check to see if you get everything(LC.fasta, LR.fasta, reference.fasta). 
        
        cd MetaFinisherSC/        

        mkdir dataFolder
        
        cd dataFolder
        
        wget -O testData.zip http://www.eecs.berkeley.edu/~kakitone/testData.zip 
        
        unzip testData.zip
        
        ls -lt

3. Go back to the main directory containing MetaFinisherSC and run the following two commands. Here we assume that /usr/bin/ is the path to your MUMmer. This step can take a few minutes. 
        
        cd ../

        python -m srcRefactor.misassemblyFixerLib.mFixer dataFolder/ /usr/bin/
        
        python -m srcRefactor.repeatPhaserLib.aSplitter  dataFolder/ /usr/bin/

4. Now verify that the original LC.fasta contain 2 contigs whereas the final abun.fasta contains 2 contigs

        fgrep -o ">" dataFolder/LC.fasta  | wc -lc
        
        fgrep -o ">" dataFolder/abun.fasta  | wc -lc


5. Align the abun.fasta with reference.fasta. 

        nucmer  -maxmatch dataFolder/abun.fasta dataFolder/reference.fasta         

        show-coords out.delta

        |   [S1]     [E1]  |     [S2]     [E2]  |  [LEN 1]  [LEN 2]  |  [% IDY]  | [TAGS] 
        |=====================================================================================
        | 2500006  2512005  |  2500001  2512000  |    12000    12000  |   100.00  | Segkk1	Segkk0
        |       1  5000002  |  5000000        1  |  5000002  5000000  |    99.99  | Segkk0	Segkk0
        |       1  5000005  |        1  5000000  |  5000005  5000000  |    99.99  | Segkk1	Segkk1
        | 2488001  2500002  |  2512000  2500001  |    12002    12000  |    99.83  | Segkk0	Segkk1

6. As a check, you may also want to see that the original contigs are misassembled. This can be seen by 

        nucmer  -maxmatch dataFolder/LC.fasta dataFolder/reference.fasta         

        show-coords out.delta

        |    [S1]     [E1]  |     [S2]     [E2]  |  [LEN 1]  [LEN 2]  |  [% IDY]  | [TAGS]
        |=====================================================================================
        |       1  2512000  |        1  2512000  |  2512000  2512000  |   100.00  | Segkk0	Segkk0
        | 2500001  5000000  |  2500001  5000000  |  2500000  2500000  |   100.00  | Segkk1	Segkk0
        |       1  2512000  |        1  2512000  |  2512000  2512000  |   100.00  | Segkk1	Segkk1
        | 2500001  5000000  |  2500001  5000000  |  2500000  2500000  |   100.00  | Segkk0	Segkk1

7. Also, you may also want to check that there is really a common segment of length 12000 across species

        nucmer  -maxmatch dataFolder/reference.fasta dataFolder/reference.fasta         

        show-coords out.delta
        
        |    [S1]     [E1]  |     [S2]     [E2]  |  [LEN 1]  [LEN 2]  |  [% IDY]  | [TAGS]
        |=====================================================================================
        |       1  5000000  |        1  5000000  |  5000000  5000000  |   100.00  | Segkk0	Segkk0
        | 2500001  2512000  |  2500001  2512000  |    12000    12000  |   100.00  | Segkk1	Segkk0
        |       1  5000000  |        1  5000000  |  5000000  5000000  |   100.00  | Segkk1	Segkk1
        | 2500001  2512000  |  2500001  2512000  |    12000    12000  |   100.00  | Segkk0	Segkk1

