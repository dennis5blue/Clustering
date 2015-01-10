#!/usr/bin/python2.7

from subprocess import Popen, PIPE
from subprocess import call
import re
import sys

def main(argv=None):
	pngSeqPath = "../sourceData/image_thesis/"
	yuvSeqPath = "../sourceData/image_yuv/"
	
	output = Popen(["ls", pngSeqPath],stdout=PIPE).communicate()[0] 
	fileList = output.split('\n')
	fileList = filter(None,fileList)
	fileList.pop() # the last item is the position information txt file 
 
	processList = len(fileList)*[None] 
	for i in range(len(fileList)):
		idx = int(re.sub(r'camera_|\.png','',fileList[i]))
		processList[idx] = re.sub(r'\.png','',fileList[i])
		#print (processList[idx])
		#print (len(processList))

	for i in range(len(fileList)):
        	print pngSeqPath+processList[i]+".png"
        	print yuvSeqPath+processList[i]+".yuv"
        	call(["../../ffmpeg","-i", pngSeqPath+processList[i]+".png","temp.mp4"])
        	call(["../../ffmpeg","-i", "temp.mp4",yuvSeqPath+processList[i]+".yuv"])
        	call(["rm","temp.mp4"])

if __name__ == "__main__":
	sys.exit(main())
