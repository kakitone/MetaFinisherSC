# MetaFinisherSC
### Mixed read type and mixed contig type###

Assume you have short reads (SR.fasta), long reads(LR.fasta) , contigs formed from short reads(SC.fasta), contigs formed from long reads(LC.fasta). 

The final output is **abun.fasta** after running through the M-Fixer(mis-assembly fixing), FinisherSC and A-Splitter(join contigs based on abundance information). 

In order to do that, here are the steps. 

1. Run M-Fixer which merge contigs with information from reads. It also fixes misassembly

        python -m srcRefactor.misassemblyFixerLib.mFixer -par 20 destinedFolder mummerPath 

2. Run FinisherSC to join the contigs together based on overlap information

        python -m srcRefactor.repeatPhaserLib.finisherSCCoreLib.finisherSC -par 20 destinedFolder mummerPath
            
3. Run A-Splitter to join contigs based on abundance information (to run faster, you can use the -ar True -rs 0 -rd True )

        python -m srcRefactor.repeatPhaserLib.aSplitter -par 20 -rp improved2.fasta  destinedFolder mummerPath



### Only long reads and contigs formed from long reads ###

Assume you have your long reads (LR.fasta) and contigs formed from long reads (LC.fasta). 

The final output is **abun.fasta** after running through the M-Fixer(mis-assembly fixing), FinisherSC and A-Splitter(join contigs based on abundance information). 

In order to do that, here are the steps. 

1. Run M-Fixer which merge contigs with information from reads(key intermediate output here is noEmbed.fasta). It also fixes misassembly

        python -m srcRefactor.misassemblyFixerLib.mFixer -par 20 -t LR destinedFolder mummerPath 


2. Run FinisherSC to join the contigs together based on overlap information(key intermediate output here is improved2.fasta)

        python -m srcRefactor.repeatPhaserLib.finisherSCCoreLib.finisherSC -par 20 destinedFolder mummerPath


3. Run A-Splitter to join contigs based on abundance information (to run faster, you can use the -ar True -rs 0 -rd True ) (key intermediate output here is abun.fasta)

        python -m srcRefactor.repeatPhaserLib.aSplitter -par 20 -rp improved2.fasta destinedFolder mummerPath



### Example ###
Below is a step by step example on running MetaFinisherSC on the testset provided(for simplicity we use the only long read mode here). 

1. Clone MetaFinisherSC
        
        git clone https://github.com/kakitone/MetaFinisherSC.git
        
2. Get into the MetaFinisherSC folder, create directory to hold the data, download data from Dropbox and check to see if you get everything(LC.fasta, LR.fasta, reference.fasta). 
        
        cd MetaFinisherSC/        

        mkdir dataFolder
        
        cd dataFolder
        
        wget -O testset.zip  https://www.dropbox.com/sh/slahcgv55037aiv/AAAzkaP6HILdBlH2G_xapKYHa?dl=1
        
        unzip testset.zip
        
        ls -lt

3. Go back to the main directory containing MetaFinisherSC and run the following two commands. Here we assume that /usr/bin/ is the path to your MUMmer. This step can take a few minutes. 
        
        cd ../

        python -m srcRefactor.misassemblyFixerLib.mFixer -par 20 -t LR dataFolder/ /usr/bin/
        
        python -m srcRefactor.repeatPhaserLib.finisherSCCoreLib.finisherSC -par 20  dataFolder/ /usr/bin/
        
        python -m srcRefactor.repeatPhaserLib.aSplitter -par 20 -rp improved2.fasta -ar True -rs 0 -rd True dataFolder/ /usr/bin/

4. Now verify that the original LC.fasta contain 4 contigs whereas the final abun.fasta contains 2 contigs

        fgrep -o ">" dataFolder/LC.fasta  | wc -lc
        
        fgrep -o ">" dataFolder/abun.fasta  | wc -lc


5. Align the abun.fasta with reference.fasta. 

        nucmer  -maxmatch dataFolder/abun.fasta dataFolder/reference.fasta         

        show-coords out.delta

6. You should now see 

        NUCMER
        
        |[S1]     [E1]  |     [S2]     [E2]  |  [LEN 1]  [LEN 2]  |  [% IDY]  | [TAGS] |
        |---------------|--------------------|--------------------|-----------|--------|
        |       1  4999997  |  4999999        1  |  4999997  4999999  |    99.99  | Segkk0	Segkk0|
        | 2490000  2500003  |  2510000  2500001  |    10004    10000  |    99.59  | Segkk1	Segkk0|
        |       1  5000004  |  4999999        1  |  5000004  4999999  |    99.99  | Segkk1	Segkk1|
        | 2490000  2499997  |  2510000  2500001  |     9998    10000  |    99.29  | Segkk0	Segkk1|


7. As a check, you may also want to see that there is really a 10K long repeat in the reference across the species. This can be seen by 

        nucmer  -maxmatch dataFolder/reference.fasta dataFolder/reference.fasta         

        show-coords out.delta



        NUCMER
        
        |[S1]     [E1]  |     [S2]     [E2]  |  [LEN 1]  [LEN 2]  |  [% IDY]  | [TAGS] |
        |---------------|--------------------|--------------------|-----------|--------|
        |       1  5000000  |        1  5000000  |  5000000  5000000  |   100.00  | Segkk0	Segkk0|
        | 2500001  2510000  |  2500001  2510000  |    10000    10000  |   100.00  | Segkk1	Segkk0|
        |       1  5000000  |        1  5000000  |  5000000  5000000  |   100.00  | Segkk1	Segkk1|
        | 2500001  2510000  |  2500001  2510000  |    10000    10000  |   100.00  | Segkk0	Segkk1|
