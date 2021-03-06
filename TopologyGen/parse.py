#!/usr/bin/python2.7

import re
import sys
def channelGain(x, y):
    g = 10**(-13.11-4.281)*(x**x+y**y)**(1/2)/1000
    return g

def main(argv=None):
	indepFilePath = "../sourceData/image_correlation/image_phase_0.out"
	corrFilePath = "../sourceData/image_correlation/image_phase_1_camera_"
	saveFilePath = "../sourceData/image_correlation/"
	posFilePath = "../sourceData/image_png/log.txt"

	W = 180 # kHz bandwidth
	P = 10  # dbm power
	indepCompression = open(indepFilePath,"r")
	ifPosition = open(posFilePath,'r')
	indepRate = open(saveFilePath+"indepRateAndByte.out","w+")
	ofPos = open(saveFilePath+"pos.out","w+")

	lines = indepCompression.readlines();
	for line in lines:
		tokens = re.split(' |\[',line)
		tokens = filter(None,tokens);
		indepRate.write(str(int(tokens[7]))+'\n') #tokens[5] is Rate, tokens[7] is number of Bytes when indep encoding

	lines = ifPosition.readlines()
	for line in lines:
		line.strip()
		tokens = line.split()
		if tokens[0] == "Pos":
			xx = float(tokens[1])*100
			yy = float(tokens[2])*100
			ofPos.write(str(xx)+' '+str(yy)+'\n')

	for camera in range(20):
		corrCompression = open(corrFilePath+str(camera)+".out","r")
		corrMatrix = open(saveFilePath+"corrMatrix","a+")
		lines = corrCompression.readlines()
		for line in lines:
			tokens = re.split(' |\[',line)
			tokens = filter(None,tokens)
			corrMatrix.write(str(int(tokens[9]))+" ")
		corrMatrix.write('\n') #c_{ij} is obtain by Ref j to code i

if __name__ == "__main__":
	sys.exit(main())
