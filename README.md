anamul
======

MG fit, FFT fit, Detector response fit for SiPM spes and DCR

Programs and usage
==================

<<getfit>>: c program do multi-gaussian and FFT fit on spes-<i>.root files, get DCR and XT from spes-<i>-dcr.root files. Make sure the same voltage applied to same <i> file, otherwise DCR is not correct calculated. All output parameters are in <DIR>/spes.log
 
usage: getfit <DIR> <nfiles> (<nfiles> can be larger than n files in the folder)

<<getfitsingle>>: c program do multi-gaussian and FFT fit on <file>.root, leave the session in ROOT. Use this one to fine tune fitting parameters.

usage: getfitsingle <DIR> <filename> [ped] [gain] [npks]
(give the 0pe peak, gain,(in qdcch) and number of peaks, otherwise the program try to find by itself)

<<SiPMSpesFit.py>>: use predefined parameters do detector response fit on spes-<i>.root files. <spes.log> file must exist in the folder. All output parameters are written in <DIR>/pyspes.log, all fitted histograms are in <DIR>out.pdf 

usage: ./SiPMSpesFit.py <DIR>

<<SiPMPlot.py>>: plot the results, command out the functions in the file

usage: ./SiPMPlot.py <DIR>