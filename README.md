anamul
======

MG fit, FFT fit, Detector response fit for SiPM spes and DCR

Programs and usage
==================
. is anamul folder

<b>./bin/getfit</b>: c program do multi-gaussian and FFT fit on <em>spes-[i].root</em> files, get DCR and XT from <em>spes-[i]-dcr.root</em> files. Make sure the <b>same voltage</b> applied to <em>spes-[i].root</em> and <em>spes-[i]-dcr.root</em> file, otherwise DCR is not correctly calculated. All output parameters are in &lt;DIR>/spes.log
 
usage: getfit &lt;DIR> &lt;n-files> (&lt;n-files> can be larger than actual n files in the folder)

<<getfitsingle>>: c program do multi-gaussian and FFT fit on <file>.root, leave the session in ROOT. Use this one to fine tune fitting parameters.

usage: getfitsingle <DIR> <filename> [ped] [gain] [npks]
(give the 0pe peak, gain,(in qdcch) and number of peaks, otherwise the program try to find by itself)

<<SiPMSpesFit.py>>: use predefined parameters do detector response fit on spes-[i].root files. <spes.log> file must exist in the folder. All output parameters are written in <DIR>/pyspes.log, all fitted histograms are in <DIR>out.pdf 

usage: ./SiPMSpesFit.py <DIR>

<<SiPMPlot.py>>: plot the results, command out the functions in the file

usage: ./SiPMPlot.py <DIR>