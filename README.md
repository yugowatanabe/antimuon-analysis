# antimuon-analysis
Python analysis script to demonstrate parity violation of weak decay of cosmic antimuons from data.

EXPERIMENT:
Read lab report (report.pdf) for detailed explanation of experiment. In essence, experiment setup captures antimuons in a 
copper target with a horizontal magnetic field applied through it--this causes the antimuons' magnetic moment to precess about it. 
The muons decay into a positron, electron neutrino, and a muon antineutrino. Scintillator panels with photomultiplier tubes (PMTs) 
above and below the target detect whether the antimuon decays upwards or downwards. If decay assymmetry is observed, parity 
violation, what was once considered a symmetry of the universe, is violated.

DETAILS:
- A separate C program was written previously by an unknown author to read data collected by the experiment probes and is not 
included in the repository (instead, I include example outputs of the program from real data collected).
- muplusAnalysis.py creates several graphs with results from the data.
- If some of the options for the script are confusing, refer to report.pdf to better understand the fundamentals of the analysis 
process.

FILES:
- 370hrsTinMon.txt is 370 hours of data WITH the copper target WITH a magnetic field applied.
- 268hrsTinMoff.txt is 268 hours of data WITH the copper target WITHOUT a magnetic field applied.
- 172hrsTout.txt is 172 hours of data WITHOUT the copper target (and therefore no magnetic field). Used for background subtraction.

EXAMPLE:
Run:
     python muplusAnalysis.py 370hrsTinMon.txt --bgsub=172hrsTout.txt 
     
as an example result of the script.
