import os

commandList = []

### Download data
commandList.append("wget https://berkeley.box.com/shared/static/m3x4czbxbh87vq078li3llr3w0ba3i8n.tar  -O dataset0.tar")
commandList.append("tar -xvf dataset0.tar")

commandList.append("wget https://berkeley.box.com/shared/static/pq98sgh3lb6vqvxrsczfrphidf1ik1gj.tar -O dataset1.tar")
commandList.append("tar xvf dataset1.tar")

commandList.append("wget https://berkeley.box.com/shared/static/h05c4hl6fk8h8k95w0y17jtkeuzsz0kp.tar -O dataset2.tar")
commandList.append("tar xvf dataset2.tar")

commandList.append("wget https://berkeley.box.com/shared/static/7ub5gqnd61pfvay6knl0ipklmrbzcsjq.tar -O dataset3.tar")
commandList.append("tar xvf dataset3.tar")


### Install tools 
commandList.append("wget http://sourceforge.net/projects/mummer/files/mummer/3.23/MUMmer3.23.tar.gz/download -O ./mummer.tar.gz")
commandList.append("tar -zxvf ./mummer.tar.gz -C ./")
commandList.append("make install -C ./MUMmer3.23/")
commandList.append("wget https://sourceforge.net/projects/quast/files/quast-2.3.tar.gz/download -O quast-2.3.tar.gz")
commandList.append("tar -zxvf quast-2.3.tar.gz ")


### Run tools 
for i in range(4):
	commandList.append("python -m srcRefactor.misassemblyFixerLib.mFixer dataset"+str(i) +"/ ./MUMmer3.23/")
	commandList.append("python -m srcRefactor.repeatPhaserLib.aSplitter  dataset"+str(i) +"/ ./MUMmer3.23/")
	commandList.append("python ./quast-2.3/quast.py -R dataset"+str(i)+"/reference.fasta dataset"+str(i) +"/LC.fasta dataset"+str(i) +"/mFixed.fasta dataset"+str(i) +"/abun.fasta")



### Final print
commandList.append("find . -name 'report.txt'  | xargs cat > allinoneeport.txt ")

for eachCommand in commandList : 
    os.system(eachCommand)







