#!/usr/bin/python2.7

from subprocess import Popen, PIPE
from subprocess import call
import re
import sys

encoderPath = "/home/dennisyu/Documents/WMSN/JMVC/bin/H264AVCEncoderLibTestStaticd"
configurePath = "./cfgFile/"

# Independent coding
def idtCoding(processList, outFile):
    output = ""
    for i in range(len(processList)):
        call(["cp",testSeqPath+processList[i], "../../JMVCtempFiles/tmp_0.yuv"])
	print(testSeqPath+processList[i])
        #output = Popen([encoderPath,"-vf","01Ex4Idt.cfg","0"],stdout=PIPE).communicate()[0]
	output = Popen([encoderPath,"-vf",configurePath+"01Ex4Idt.cfg","0"],stdout=PIPE).communicate()[0]
	print (output)        
	print "encode: ", i
	lines = [line for line in output.split('\n') if "byte" in line]
	byteInfo = " ".join(lines[0].strip().replace("[","").replace("]","").split())
	outFile.write(" Code: "+str(i)+" "+byteInfo+"\n")
	
# Correlation discovering
def corrDiscover(camIndex, processList, outFile):
	encodeCamera = processList[int(camIndex)]
	print "Encoding Camera: ",encodeCamera
	
	for i in range(len(processList)):
		call(["rm","../../JMVCtempFiles/cam"+camIndex+"rec_0.yuv","../../JMVCtempFiles/cam"+camIndex+"rec_1.yuv"])
		call(["rm","../../JMVCtempFiles/cam"+camIndex+"stream_0.264","../../JMVCtempFiles/cam"+camIndex+"stream_1.264"])
		call(["rm","../../JMVCtempFiles/cam"+camIndex+"tmp_0.yuv","../../JMVCtempFiles/cam"+camIndex+"tmp_1.yuv"])
		call(["cp",testSeqPath+processList[i], "../../JMVCtempFiles/cam"+camIndex+"tmp_0.yuv"])
		output_0 =  Popen([encoderPath,"-vf",configurePath+"01cam"+camIndex+"Diff.cfg","0"],stdout=PIPE).communicate()[0]
		call(["cp",testSeqPath+encodeCamera, "../../JMVCtempFiles/cam"+camIndex+"tmp_1.yuv"])
		output_1 =  Popen([encoderPath,"-vf",configurePath+"01cam"+camIndex+"Diff.cfg","1"],stdout=PIPE).communicate()[0]
		
		print "ref: ", i,"code: ", camIndex
		lines = [line for line in output_1.split('\n') if "byte" in line]
		byteInfo = " ".join(lines[0].strip().replace("[","").replace("]","").split())
		outFile.write("Ref "+str(i)+" Code: "+str(camIndex)+" "+byteInfo+"\n")



def main(argv = None): 
    if sys.argv[1].startswith('--'):
		option = sys.argv[1][2:]
		if option == "phase":
			global testSeqPath
			global fileNamePrefix
			global saveFilePath
			testSeqPath = "../sourceData/image_yuv/"
			fileNamePrefix = "image"
			saveFilePath = "../sourceData/image_correlation/"
		if option == "help":
			print "Usage: --phase [0|1]: 0 for independent coding, 1 for correlation discovering"
			print "Usage: --camera k: verify which camera we aim to encode (no need to use this command when indep coding)"
			exit()
		
		output = Popen(["ls", testSeqPath],stdout=PIPE).communicate()[0] 
		fileList = output.split('\n')
		fileList = filter(None,fileList)
		#fileList.pop() # the last item is the position information txt file 
		processList = len(fileList)*[None]
		#print (fileList)
		for i in range(len(fileList)):
			idx = int(re.sub(r'camera_|\.yuv','',fileList[i]))
			processList[idx] = fileList[i]
			#print (processList)
		if sys.argv[2] == "0":
			call(["rm","../../JMVCtempFiles/rec_0.yuv","../../JMVCtempFiles/rec_1.yuv"])
			outFile = open(saveFilePath+fileNamePrefix+"_phase_"+str(0)+".out","w+")
			result = idtCoding(processList, outFile)
		if sys.argv[2] == "1":
			if sys.argv[3] == "--camera":
				camIndex = sys.argv[4]
				#print "Differential encode camera", camIndex
			else:
				print "Please verify which camera are encoding now"
				exit()
			call(["rm","../../JMVCtempFiles/rec_0.yuv","../../JMVCtempFiles/rec_1.yuv"])
			outFile = open(saveFilePath+fileNamePrefix+"_phase_"+str(1)+"_camera_"+str(camIndex)+".out","w+")
			#print (processList)
			corrDiscover(camIndex, processList, outFile)
        
if __name__ == "__main__":
    sys.exit(main())
