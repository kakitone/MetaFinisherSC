# MetaFinisherSC
### Introduction ###
To run the tool, you should have input files(contigs.fasta, raw_reads.fasta) in destinedFolder. The tool will first run through FinisherSC to generate improved2.fasta to destinedFolder. Then it will go through the newly added functions suitable for metagenomics study and the final output file is **abun.fasta**
The command to run MetaFinisherSC to use 20 threads(-par 20) is thus as follows. 

        python finisherSCCoreLib/finisherSC.py -par 20 destinedFolder/ mummerPath/

        python aSplitter.py -par 20 -rp improved2.fasta -ar True -rs 0 -rd True destinedFolder/ mummerPath/


A test set is located here https://www.dropbox.com/sh/slahcgv55037aiv/AAAzkaP6HILdBlH2G_xapKYHa?dl=0  

If you have enough computing resources, you can run 

        python finisherSCCoreLib/finisherSC.py -par 48 destinedFolder/ mummerPath/

        python aSplitter.py -par 48 -rp improved2.fasta destinedFolder/ mummerPath/

### Example ###
Below is a step by step example on running MetaFinisherSC on the testset provided. 
1. Clone MetaFinisherSC
        
        git clone https://github.com/kakitone/MetaFinisherSC.git
        
2. Get into the MetaFinisherSC folder, create directory to hold the data, download data from Dropbox and check to see if you get everything(contigs.fasta, raw_reads.fasta, reference.fasta). 
        
        cd MetaFinisherSC/        

        mkdir dataFolder
        
        wget -O testset.zip  https://www.dropbox.com/sh/slahcgv55037aiv/AAAzkaP6HILdBlH2G_xapKYHa?dl=1
        
        unzip testset.zip
        
        ls -lt

3. Go back to the main directory containing MetaFinisherSC and run the following two commands
        
        python finisherSCCoreLib/finisherSC.py -par 20 dataFolder/ /usr/bin/
        
        python aSplitter.py -par 20 -rp improved2.fasta -ar True -rs 0 -rd True dataFolder/ /usr/bin/

4. Now verify that the original contigs.fasta contain 4 contigs whereas the final abun.fasta contains 2 contigs

        fgrep -o ">" datafolder/contigs.fasta  | wc -l
        
        fgrep -o ">" dataFolder/abun.fasta  | wc -l


5. Align the abun.fasta with reference.fasta. 

        nucmer dataFolder/abun.fasta dataFolder/reference.fasta         

        show-coords out.delta

6. You should now see 

NUCMER

[S1]     [E1]  |     [S2]     [E2]  |  [LEN 1]  [LEN 2]  |  [% IDY]  | [TAGS]
=====================================================================================
1  5000000  |  5000000        1  |  5000000  5000000  |   100.00  | Segkk0	Segkk0
1  5000000  |  5000000        1  |  5000000  5000000  |   100.00  | Segkk1	Segkk1

