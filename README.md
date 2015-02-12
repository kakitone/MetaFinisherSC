# workingMetaFinisherSC
=======
# MetaFinisherSC
To run the tool, you should have input files(contigs.fasta, raw_reads.fasta) in destinedFolder. The tool will first run through FinisherSC to generate improved2.fasta to destinedFolder. Then it will go through the newly added functions suitable for metagenomics study and the final output file is **abun.fasta**
The command to run MetaFinisherSC to use 20 threads(-par 20) is thus as follows. 

        python finisherSCCoreLib/finisherSC.py -par 20 destinedFolder/ mummerPath/

        python aSplitter.py -par 20 -rp improved2.fasta -ar True -rs 0 -rd True destinedFolder/ mummerPath/


A test set is located here https://www.dropbox.com/sh/slahcgv55037aiv/AAAzkaP6HILdBlH2G_xapKYHa?dl=0  

If you have enough computing resources, you can run 

        python finisherSCCoreLib/finisherSC.py -par 48 destinedFolder/ mummerPath/

        python aSplitter.py -par 48 -rp improved2.fasta destinedFolder/ mummerPath/
