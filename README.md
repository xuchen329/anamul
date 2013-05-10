anamul
======

MG fit, FFT fit, Detector response fit for SiPM spes and DCR

Programs and usage
==================
. is anamul folder

<b>./bin/getfit &lt;DIR> &lt;n-files></b>
    &lt;DIR> is the folder contains <em>spes-[i].root</em> and <em>spes-[i]-dcr.root</em> files, i in the range of  [0,&lt;n-files>], files not exist will be ignored. The program do multi-gaussian and FFT fit on <em>spes-[i].root</em> files, get DCR and XT from <em>spes-[i]-dcr.root</em> files. Make sure the <b>same voltage</b> applied to <em>spes-[i].root</em> and <em>spes-[i]-dcr.root</em> file, otherwise DCR is not correctly calculated. All output parameters are in <em>&lt;DIR>/spes.log</em>

<b>./bin/getfitsingle &lt;path/filename> [ped] [gain] [npks]</b>
    Run the same fit on single root file and leave the session in ROOT afterwards. 
    An entry is printed in the end, can be used to substitute the entry in <em>&lt;DIR>/spes.log</em> directly. 
    <em>[ped] [gain]</em> (must be given together) will be used instead of finding peaks. 
    <em>[npks]</em> tells the program how many peaks <b>after 0pe</b> are expected in the spes histogram. All in qdc channel number.

<b>./bin/getdcr &lt;DIR> &lt;n-file> &lt;fitpar0> &lt;fitpar1> &lt;effgate>
    Get DCR from <em>spes-[i]-dcr.root</em> files in <em>&lt;DIR></em>, use this in case the dcr are measured seperately.
    The program determine 0.5thr from the linear dependance of gain on Vbias.
    <em>&lt;fitpar0> &lt;fitpar1></em> are fit parameters in <em>Gain = &lt;fitpar0>+Vbias*&lt;fitpar1></em>, all in qdcch unit. Mind that Vbias may be negtive.
    <em>&lt;effgate></em> is the effective QDC gate in DCR measurement.

<b>./pythonfit/SiPMSpesFit.py &lt;DIR></b>
    <em>&lt;spes.log></em> generated by <b>getfit</b> must exist in <em>&lt;DIR></em>
    The program will then apply detector response fit on the root files. 
    Output parameters are written in <em>&lt;DIR>/pyspes.log</em>, fitted histograms are printed in <em>&lt;DIR>/out.pdf</em> 

<b>./pythonfit/SiPMPlot.py &lt;DIR></b>
    plot the results. results are taken from <em>pyspes.log</em> or <em>spes.log</em>