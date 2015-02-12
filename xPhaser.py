import time
import sys
import os
import polyPhaser

t0 = time.time()
print 'Number of arguments:', len(sys.argv), 'arguments.'
print 'Argument List:', str(sys.argv)
print os.path.abspath(os.path.dirname(sys.argv[0]))  
folderName = sys.argv[1]
mummerLink = sys.argv[2]


polyPhaser.mainFlow(folderName, mummerLink)
print  "Time", time.time() - t0
  
