#!/bin/bash



# Open  mobaexterm shell from local laptop
#Run the transfer files GBM_transfer_codes.sh


figDir=Fig6

for f in /projects/compsci/jgeorge/USERS/alizae/GBM/matching_patients/R_Code/sc_pipeline/$figDir/*/*/Configs/*GEX_preprocessing*.R;do 
echo renaming $(basename $f):
mv $f $(dirname $f)/GEX_preprocessing.R 
done



for f in /projects/compsci/jgeorge/USERS/alizae/GBM/matching_patients/R_Code/sc_pipeline/$figDir/*/*/Configs/*ADT_preprocessing*.R;do
echo renaming $(basename $f):
mv $f $(dirname $f)/ADT_preprocessing.R
done


for f in /projects/compsci/jgeorge/USERS/alizae/GBM/matching_patients/R_Code/sc_pipeline/$figDir/*/*/Configs/*standard_viz*.R;do
echo renaming $(basename $f):
mv $f $(dirname $f)/standard_viz.R
done


for f in /projects/compsci/jgeorge/USERS/alizae/GBM/matching_patients/R_Code/sc_pipeline/$figDir/*/*/Configs/*utils*.R;do
echo renaming $(basename $f):
mv $f $(dirname $f)/utils.R
done


for f in /projects/compsci/jgeorge/USERS/alizae/GBM/matching_patients/R_Code/sc_pipeline/$figDir/*/*/Configs/*analysis*.json;do
echo renaming $(basename $f):
mv $f $(dirname $f)/analysis.json
done


for f in /projects/compsci/jgeorge/USERS/alizae/GBM/matching_patients/R_Code/sc_pipeline/$figDir/*/*/Configs/*config*.json;do
echo renaming $(basename $f):
mv $f $(dirname $f)/config.json
done


for f in /projects/compsci/jgeorge/USERS/alizae/GBM/matching_patients/R_Code/sc_pipeline/$figDir/*/*/Configs/*sessionInfo.txt;do
echo removing $(basename $f):
rm $f 
done


for f in /projects/compsci/jgeorge/USERS/alizae/GBM/matching_patients/R_Code/sc_pipeline/*/*/ADTResAll/Configs/GEX_preprocessing.R;do
echo renaming $(basename $f):
mv $f $(dirname $f)/ADT_preprocessing.R
done
