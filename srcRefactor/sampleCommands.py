1) Run MFixer with options 

python -m srcRefactor.misassemblyFixerLib.mFixer May13Test/ /usr/bin/ 


nohup python -m srcRefactor.misassemblyFixerLib.mFixer \
-t LR  -par 20  \
-op "toRunAggressive=False" \
May13Test/ /usr/bin/ &


2) Run ASplitter with options

python -m srcRefactor.repeatPhaserLib.aSplitter May11TestB/ /usr/bin/ 


nohup python -m srcRefactor.repeatPhaserLib.aSplitter \
-par 20 -rp mFixed.fasta  -ar True -rs 0 -rd False \
-pk map \
May13Test/ /usr/bin/ &

3) Run Evaluator to print headers or do evalulation 

python -m srcRefactor.evaluator --option header \
--quastPath /data/kakitone/download2/quast-2.3/quast.py \
--folderName /data/kakitone/May07-2015/workingMetaFinisherSC/May11TestB/ \
--paraMFixerFileName optionMFixer.json \
--paraASplitterFileName option.json \
--outputFilename /data/kakitone/May07-2015/workingMetaFinisherSC/results.csv &


nohup python -m srcRefactor.evaluator --option evaluate \
--quastPath /data/kakitone/download2/quast-2.3/quast.py \
--folderName /data/kakitone/May07-2015/workingMetaFinisherSC/May11TestB/ \
--paraMFixerFileName optionMFixer.json \
--paraASplitterFileName option.json \
--outputFilename /data/kakitone/May07-2015/workingMetaFinisherSC/results.csv &

