#Blazar Variability / BlaVar
Scripts: python3 , bash


###Here are the specific scripts for data-analysis of the time-dependent one-zone leptonic analysis. 
CAUTION 1 ! Directories are structured with a specific format. Please adjust to your environment.
CAUTION 2 ! Many scripts are employed for analyzing data have dependencies. Specifically,
there are dependencies on "code_clean", "code_var" fortran codes for time-dependent radiative tran-
sfer and executables and on the Emmanoulopoulos code (Emmanoulopoulos et al. 2013), adopted in the 
Python environment as a "DELC" library.
-------------------------------------------------------------------------------------------------------


1) Obs4fit_v2.py :  Finds steady state parameters using input observables: slope of uncooled electrons, 
maximum emitting frequency, ratio of synchrotron/compton luminosities. Manual values can be used in 
combination of the suggested parameters. There is an option to run the steady state code (10 tcross) &
compare with the observed SED data points (chi^2). The results(.inp & .ps) are saved with a specific 
number (date) using  " bash ./run_copy.sh" in the /Steady_States folder.
See documentation in file.

2)Compare2timecurves.py : Used for comparing two LCs computing the DCF and FVs. Mainly used to compare
original and edited Fermi-LAT LCs.

3) Thesis_Plotter.py : Used to plot combined plots from different simulations. The code architecture
allows for adding as many simulation results as desired. Fractional Variabilities are inserted by-hand
from simulation results (Results_analysis.py). Dependencies: CI.txt, FERMI.txt, fakeLC time curves in the 
/Results/Name_obj/Type_of_sim folders of the simulations used.

4) /Lightcurve_Simulation/LC_analysis.py: Used automatically from BlaVar_single.sh and BlaVar_multi.sh, but
can be employed independently.
Uses real Fermi-LAT energy flux(.fits files) to produce  fake LCs. Also, edites the real LCs (corrects 
for extremely low fluxes) and creates parameter time-series which is then used as time-dependent input
for modelling.
Input file "fkTC.inp" in same folder (provided by other script)
See documentation in file.

5)/Var/run_var.sh Given the steady state parameters /object/code_new.inp,steady.inp + data points 
/object/name.ascii, it runs the full simulation for the requested time-interval and save results in 
/Results/object folder. Caution to the restart switch! Tries different tolerence
values if code stucks numerically at some point. 
Employs LC_analysis.py, code_clean, code_var, Results_analysis.py.
PROVIDE DIRECTORY IN FILE.


6)/Var/BlaVar_multi.sh If run_var.sh stucks at point. Simply use BlaVar_multi.sh to create multiple 
fort_index.81,85,89,55 files that you will eventually add together to a single timecurve. "multi.log" file 
contain the time (in tcross) the code crush and generated a new fort file. For the "gluing" of the results 
see Results_analysis.py.
Employs LC_analysis.py, code_clean, code_var, Results_analysis.py.
PROVIDE DIRECTORY IN FILE.

7) /Var/Results/Object/Name_simulation/Results_analysis.py, delta_analysis.py: Are python scripts that do all the data-
analysis (lightcurves, contour plots, timelapse, DCFs , SFs , FVs , stingray, index analysis ++)
See documentation in file.
Send to all directories updated version of /Var/Results_analysis.py or delta_analysis.py by using the ./send.sh script.

------------------------------------------------------------------------------------------------------------
