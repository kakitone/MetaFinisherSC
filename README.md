# MetaFinisherSC
### Mixed read type and mixed contig type###

Assume you have short reads (SR.fasta), long reads(LR.fasta) , contigs formed from short reads(SC.fasta), contigs formed from long reads(LC.fasta). 

The final output is **abun.fasta** after running through the mis-assembly fixing, FinisherSC and A-Splitter(join contigs based on abundance information). 

In order to do that, here are the steps. 

1. Run merger.py which merge contigs with information from reads. It also fixes misassembly

	python -m srcRefactor.misassemblyFixerLib.mFixer -par 20 destinedFolder mummerPath 

2. Run FinisherSC to join the contigs together based on overlap information

	python -m srcRefactor.repeatPhaserLib.finisherSCCoreLib.finisherSC -par 20 destinedFolder mummerPath
            
3. Run ASplitter to join contigs based on abundance information (to run faster, you can use the -ar True -rs 0 -rd True )

	python -m srcRefactor.repeatPhaserLib.aSplitter -par 20 -rp improved2.fasta  destinedFolder mummerPath



### Only long reads and contigs formed from long reads ###

Assume you have your long reads (LR.fasta) and contigs formed from long reads (LC.fasta). 

The final output is **abun.fasta** after running through the mis-assembly fixing, FinisherSC and A-Splitter(join contigs based on abundance information). 

In order to do that, here are the steps. 

1. Run merger.py which merge contigs with information from reads. It also fixes misassembly

	python -m srcRefactor.misassemblyFixerLib.mFixer -par 20 -t LR destinedFolder mummerPath 


2. Run FinisherSC to join the contigs together based on overlap information

	python -m srcRefactor.repeatPhaserLib.finisherSCCoreLib.finisherSC -par 20 destinedFolder mummerPath


3. Run ASplitter to join contigs based on abundance information (to run faster, you can use the -ar True -rs 0 -rd True )

	python -m srcRefactor.repeatPhaserLib.aSplitter -par 20 -rp improved2.fasta destinedFolder mummerPath





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
        
        cd dataFolder
        
        wget -O testset.zip  https://www.dropbox.com/sh/slahcgv55037aiv/AAAzkaP6HILdBlH2G_xapKYHa?dl=1
        
        unzip testset.zip
        
        ls -lt

3. Go back to the main directory containing MetaFinisherSC and run the following two commands. Here we assume that /usr/bin/ is the path to your MUMmer. This step can take a few minutes. 
        
        cd ../

        python finisherSCCoreLib/finisherSC.py -par 20 dataFolder/ /usr/bin/
        
        python aSplitter.py -par 20 -rp improved2.fasta -ar True -rs 0 -rd True dataFolder/ /usr/bin/

4. Now verify that the original contigs.fasta contain 4 contigs whereas the final abun.fasta contains 2 contigs

        fgrep -o ">" dataFolder/contigs.fasta  | wc -l
        
        fgrep -o ">" dataFolder/abun.fasta  | wc -l


5. Align the abun.fasta with reference.fasta. 

        nucmer dataFolder/abun.fasta dataFolder/reference.fasta         

        show-coords out.delta

6. You should now see 

NUCMER

|[S1]     [E1]  |     [S2]     [E2]  |  [LEN 1]  [LEN 2]  |  [% IDY]  | [TAGS] |
|---------------|--------------------|--------------------|-----------|--------|
|1  5000000     | 5000000        1  |  5000000  5000000  |   100.00  | Segkk0	Segkk0|
|1  5000000  |  5000000        1  |  5000000  5000000  |   100.00  | Segkk1	Segkk1|

