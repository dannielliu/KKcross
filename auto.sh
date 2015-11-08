#!/bin/bash

#cd /moose/Bes3User/liud/work/KKCross/KKcross
rm cutflow2
echo "############################" >> cutflow.dat
echo "## uncertainty tof cut  ####" >> cutflow.dat
echo "############################" >> cutflow.dat
 energys="20000 20500 21000 21500 21750 \
 22000 22324 23094 23864 23960 \
 25000 26444 26464 27000 28000 \
 29000 29500 29810 30000 30200 \
 30800"
#energys=" 23960 \
#25000 26444 26446 27000 28000 \
#29000 29500 29810 30000 30200 \
#30800"
#energys="30800"


#filelist="KK_20000.root KK_20500.root KK_21000.root KK_21500.root KK_21750.root KK_22000.root KK_22324.root KK_23094.root KK_23864.root KK_23960.root KK_25000.root KK_26444.root KK_26446.root KK_27000.root KK_28000.root KK_29000.root KK_29500.root KK_29810.root KK_30000.root KK_30800.root"

 for ene in $energys;do
   file="KK_$ene.root"
   echo "current file is $file"
   echo "current file is $file" >> log
   #mkdir -p outputcut
   ./anaKK /Volumes/data/data2/ee2KK/ee2KKv4/$file >> log
 done


#mcdir="/Volumes/data2/ee2KK_mc/root"
#mcdir="/Volumes/data2/ee2KK_mc/ee2KK_mcConExc_665/root"

#mcdir="/Volumes/data2/ee2KK_mcuserxs/root"
#mcdir="/Volumes/data2/ee2KK_mc2_3/root"
#mcdir="/Volumes/data2/ee2KK/ee2KK_Phokhara/root"
 mcdir="/Volumes/data/data2/ee2KK/ee2KK_ConExcmyxs/root"

 for file in `ls $mcdir/*.root`;do
   echo "current file is $file"
   echo "current file is $file" >> log
   ./anaKK $file >> log
 done
 
 mcdir="/Volumes/data2/ee2KK_ConExcmyxs2/root"

#for file in `ls $mcdir/*.root`;do
#  echo "current file is $file"
#  echo "current file is $file" >> log
#  ./anaKK $file >> log
#done

mcdir="/Volumes/data2/ee2KK/ee2KK_ConExc_0rad/root"

#for file in `ls $mcdir/*.root`;do
#  echo "current file is $file"
#  echo "current file is $file" >> log
#  ./anaKK $file >> log
#done



#mcdir="/moose/Bes3User/liud/ee2KK_Phokhara_born/root"
#mcdir="/Volumes/data2/mc2015/ee2KK_BB/root"

#for file in `ls $mcdir/*.root`;do
#  echo "current file is $file"
#  echo "current file is $file" >> log
#  ./AnaCuts $file >> log
#done

#mcdir="/Volumes/data2/mc2015/ee2KK_dimu/root"

#for file in `ls $mcdir/*.root`;do
#  echo "current file is $file"
#  echo "current file is $file" >> log
#  ./AnaCuts $file >> log
#done

#mcdir="/Volumes/data2/mc2015/ee2KK_had/root"

#for file in `ls $mcdir/*.root`;do
#  echo "current file is $file"
#  echo "current file is $file" >> log
#  ./AnaCuts $file >> log
#done

#mcdir="/Volumes/data2/ee2KK_ConExcmyxs/root"

#for file in `ls $mcdir/*.root`;do
#  echo "current file is $file"
#  echo "current file is $file" >> log
#  ./AnaCuts $file >> log
#done

#mcdir="/Volumes/data2/ee2KK_mcuserxs/root"

#for file in `ls $mcdir/*.root`;do
#  echo "current file is $file"
#  echo "current file is $file" >> log
#  ./AnaCuts $file >> log
#done


cat cutflow2 >> cutflow.dat
