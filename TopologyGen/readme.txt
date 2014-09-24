This is document for generating correlation matrix

 1. Execute the blender python script (cameraClusterScenario.py) and get the test sequence, and log.txt
	Saved in TestSequences/test_png
 2. Execute yuv_gen.py to convert the png file to the yuv raw data
 3. Execute encode.py to try all the combination of multi view encoding
	$python encode.py --phase 0 is to perform the independent coding
	
	Note that the encoder config file need to specified, 
    for indendent coding    : 00Ex4Idt.cfg, or 01Ex4Idt.cfg
    for differential coding : 00Ex4Diff.fig, or 01Ex4Diff.cfg
    the resoultion setting must match to the property of *.yuv file
	search range can be changed for higher performance

    $python encode.py --phase 1 --camera k is to perform the differential coding for camera k
	(Can use the encode_phase1.sh for phase 1 encoding)
 4. Execute parse.py to generate the rate.out and pos.out
	Give the path by $python parse.py
	Note that do not execute this script two times, since the written mode is append, it may cause duplicate data
 
