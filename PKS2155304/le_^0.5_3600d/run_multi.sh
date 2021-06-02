#!/bin/bash
#OPTIONS
restart=1 #1 to activate,  #have to generate LC first, try run_var.sh
steady=1 #1: new steady state
route='/home/markos/Desktop/BLAVAR/Var'
obj='PKS2155304' #no '-' or '_' or '/', raw letters and numbers
lendays=3600 #Desired time interval (days) for reproducing lightcurves
namevar='le' #Variable to be varying
POWER=0.5 #power p where: param propto (Flux_gamma)^p / 1.0 for simple use 		   #example:  for delta POWER= '0.25'
TOLinit=-4 #TOL= 1eTOLinit steady state tolerance
#
#
#
#
#
#
#SCRIPT
#Steady State (uncomment if changed object)
if [ $steady -eq 1 ]
then
rm steady*
cp ./code.inp ./temp.inp
cp ./$obj/steady.inp ./code.inp
./code_clean
cp ./temp.inp ./code.inp
rm temp.inp
cp ./fort.81 ./steady.81
cp ./fort.85 ./steady.85
cp ./fort.89 ./steady.89
echo "Steady State ready/ Var Code next"
sleep 3
fi

#used variables in script
R=$(nawk 'NR==3{print $1}' $route/$obj/code_new4.inp)
delta=$(nawk 'NR==8{print $1}' $route/$obj/code_new4.inp)
tcross=$(awk "BEGIN {print $R/3e10/24./3600.}")

if [ $restart -eq 1 ]
then
#Run code from where it has crashed or from the beginning
rm fort*
rm fakeTC*
rm SEDs.txt
rm SMARTS*
rm X_*
rm FERMI.txt
rm FRAMES.txt
rm BB*
#Create/Display_current Timecurve
remT=$(awk "BEGIN {print $lendays/$tcross*$delta}")
remT=$((${remT%.*}+1))
runnum=0
rm ~/Desktop/BLAVAR/Lightcurve_Simulation/$obj/multi.log
fi

#Continue using log
if [ $restart -eq 0 ]
then 
runnum=$(awk 'END{print NR}' multi.log)
runnum=$(($runnum-1))
echo "WARNING: Restart at the run number $runnum"
echo "For different treatment (not overwriting the crashed run) use run_var.sh"
sleep 1
lenprev=$(awk 'END{print $1}' multi.log)
remT=$(awk "BEGIN {print $lendays/$tcross*$delta-$lenprev}")
remT=$((${remT%.*}+1))
fi

while [ $remT -gt 0 ] #outer loop
do
rm fort.*
rm fakeTC.txt

newD="multi$remT"
awk -v var=$namevar 'NR==10{$1=var};{print}' $route/$obj/code_new4.inp > $route/$obj/temp.inp
awk -v var=$lendays 'NR==10{$2=var};{print}' $route/$obj/temp.inp > $route/$obj/code_new4.inp
awk -v var=$newD 'NR==10{$3=var};{print}' $route/$obj/code_new4.inp > $route/$obj/temp.inp
awk -v var=$POWER 'NR==10{$4=var};{print}' $route/$obj/temp.inp > $route/$obj/code_new4.inp
cd ~/Desktop/BLAVAR/Lightcurve_Simulation/$obj
cp $route/$obj/steady.inp ./code.inp
cp $route/$obj/code_new4.inp ./fkTC.inp
source ~/anaconda3/etc/profile.d/conda.sh
conda activate astro
python ./my_analysis_v3.py
conda deactivate
cp ./code.inp $route/$obj/code.inp
fktc=$(echo fakeTC"$lendays"_"$namevar"_"$obj".txt)
cp ./$fktc $route/fakeTC.txt
cp ./$fktc $route/fakeTC_"$runnum".txt
cp ./code_new4.inp $route/code_new4.inp
cp ./multi.log $route/multi.log
cd $route
echo "New Run $runnum initialized"
echo "Check Steepness Diagram in Lightcurve folder to find potential NAG-failure points"


TIME=$(nawk 'NR==1{print $4}' $route/$obj/code.inp)
STEPS=$(nawk 'NR==1{print $3}' $route/$obj/code.inp)
awk -v var=$TIME 'NR==1{$4=var};{print}' $route/$obj/code_new4.inp > $route/$obj/temp.inp
awk -v var=$STEPS 'NR==1{$3=var};{print}' $route/$obj/temp.inp > $route/$obj/code_new4.inp
echo "Time interval for the run created $TIME tcross"
sleep 3
crashes=0
stuck=0
SECONDS=0
tolchange=0
TOL="$TOLinit"
NlnsTC=$(awk 'END{print NR}' fakeTC.txt)
tab='      '
NTCtext="$tab ifake=$NlnsTC"
sed -i -e "205s/.*/$NTCtext/" code_var_$namevar.f
g77-3.4 -w -L ~/Downloads/nag17/ -o code_var code_var_$namevar.f -lnag -O3 -funroll-loops
cp ./$obj/code.inp ./code.inp
time=0.
mod=0

if [ "$namevar" == "gmax" ]
then
#sed -i -e "2s/7./9./1" code.inp
sed -i -e "1s/5/10/1" code.inp
fi

while [ $mod -eq 0 ]
do
oldtime=$(nawk 'NR==1{print $2}' $route/dum_r4.dat)
./code_var
awk 'NR==1{$1=1};{print}' $route/code.inp > $route/temp.inp
cp $route/temp.inp $route/code.inp
time=$(nawk 'NR==1{print $2}' $route/dum_r4.dat)
modstep=$(echo "($time-$oldtime)/($TIME/$STEPS)" |bc)
if [ $modstep -eq 0 ]
then
TOL=-2
nplus=10.
while [ $modstep -eq 0 ]
do
tolchange=$(($tolchange+1))
TOL=$(($TOL-1))
TOLs="3e$TOL"

if [ $TOL -eq -5 ]
then
nplus=1.
fi
shT=$(echo "scale=1;$time+$nplus*$TIME/$STEPS" |bc)
shstp=$(echo "$shT/($TIME/$STEPS)" |bc)
awk -v var=$TOLs 'NR==2{$5=var};{print}' $route/code.inp > $route/temp.inp
awk -v var=$shT 'NR==1{$4=var};{print}' $route/temp.inp > $route/code.inp
awk -v var=$shstp 'NR==1{$3=var};{print}' $route/code.inp > $route/temp.inp
cp ./temp.inp ./code.inp
rshT=$(echo "($shT-$time)/($TIME/$STEPS)" |bc)
awk -v var=$rshT 'NR==1{$1=var};{print}' dum_r4.dat> temp.dat
cp ./temp.dat ./dum_r4.dat 
oldtime=$(nawk 'NR==1{print $2}' $route/dum_r4.dat)
./code_var
time=$(nawk 'NR==1{print $2}' dum_r4.dat)
modstep=$(echo "($time-$oldtime)*$STEPS/$TIME" |bc)

if [ $TOL -lt -5 ]
then
modstep=1
else
TOLs="2e$TOL"
shT=$(echo "scale=1;$time+$nplus*$TIME/$STEPS" |bc)
shstp=$(echo "$shT/($TIME/$STEPS)" |bc)
awk -v var=$TOLs 'NR==2{$5=var};{print}' $route/code.inp > $route/temp.inp
awk -v var=$shT 'NR==1{$4=var};{print}' $route/temp.inp > $route/code.inp
awk -v var=$shstp 'NR==1{$3=var};{print}' $route/code.inp > $route/temp.inp
cp ./temp.inp ./code.inp
rshT=$(echo "($shT-$time)/($TIME/$STEPS)" |bc)
awk -v var=$rshT 'NR==1{$1=var};{print}' dum_r4.dat> temp.dat
cp ./temp.dat ./dum_r4.dat 
oldtime=$(nawk 'NR==1{print $2}' $route/dum_r4.dat)
./code_var
time=$(nawk 'NR==1{print $2}' dum_r4.dat)
modstep=$(echo "($time-$oldtime)*$STEPS/$TIME" |bc)
fi

if [ $TOL -lt -5 ]
then
modstep=1
fi

done
fi

if [ $TOL -gt -6 ]
then
TOL="$TOLinit"
TOLs="1e$TOL"
awk -v var=$TOLs 'NR==2{$5=var};{print}' $route/code.inp > $route/temp.inp
awk -v var=$TIME 'NR==1{$4=var};{print}' $route/temp.inp > $route/code.inp
awk -v var=$STEPS 'NR==1{$3=var};{print}' $route/code.inp > $route/temp.inp
cp ./temp.inp ./code.inp
rshT=$(echo "($TIME-$time)/($TIME/$STEPS)" |bc)
awk -v var=$rshT 'NR==1{$1=var};{print}' dum_r4.dat> temp.dat
cp ./temp.dat ./dum_r4.dat
fi

mod=$(echo "$time/$TIME" |bc)
if [ $TOL -gt -6 ]
then 
crashes=$(($crashes+1))
else
echo "Code STUCK at time $time tcross"
mod=1
fi
done

cp fort.81 "fort_$runnum.81"
cp fort.85 "fort_$runnum.85"
cp fort.89 "fort_$runnum.89"
cp fort.55 "fort_$runnum.55"
if [ "$namevar" == "lext" ]
then
cp fort.86 "fort_$runnum.86"
fi

runnum=$(echo "$runnum+1" |bc)
remT=$(echo "$TIME/1.0-$time/1.0" |bc)
echo "Simulations to conjoin so far: $runnum"
stuck=$(($stuck+1))
done #outer loop end


if [ $remT -lt 2 ]
then
echo "Code run to the endtime"
echo "Code crashed $crashes times due to NAG FAILURE"
echo "Code stucked $stuck times due to NAG FAILURE"
echo "Tolerance changed $tolchange times"
echo "Total elapsed time: $SECONDS secs"

cp $route/$obj/code.inp ./code.inp
mkdir $route/Results/$obj/"$namevar"_^"$POWER"_"$lendays"d
cd $route/Results/$obj/"$namevar"_^"$POWER"_"$lendays"d
cp $route/fort* ./
cp $route/steady* ./
cp $route/code.inp ./
cp $route/multi.log ./
cp $route/run_multi.sh ./
cp $route/fakeTC* ./
cp $route/Var4.py ./Var4.py
cp $route/code_var ./
cp $route/code_var_"$namevar".f ./
rm temp.inp
rm fort.*


echo "Results ready in directory Results/Object/Varname  // Analysis (Var4.py) next , Ctrl+Z to stop here"
sleep 3
python3 Var4.py
fi


