#!/bin/bash

####################################################
#                                                  #
# This is a compilation script of the processing   #
# scripts required to produce the Hurst exponent   #
# values on the Glasser parcellation of MRC_AIMS   #
# data.                                            #
#																									 #
# The preprocessing steps include:                 #
#  1) Manually removing 8 volumes from the data.   #
#																									 #
#  2) Running speedy for wavelet denoising,        #
# calculating DVARS and Mean Framewise             #
# displacement along with standard preprocessing   #
# steps.   																				 #
#																									 #
#  3) Resampling the Glasser 2016 template to	   	 #
# native space to produce subject specific maps    #
# that are used in the calculation of the Hurst	   #
# exponent. 																			 #
#																									 #
#  4) Calculating the Hurst exponent on the masked #
# resting data. 																	 #
#												   												 #
# Dependencies																		 #
# 	All steps depend on afni, fsl and matlab       #
# 	and python 			    													 #
#																									 #
# 	Depends on the Brain Wavelet toolbox and the   #
#   fMRI signal Preprocessing Toolbox              #
#   from www.brainwavelet.org                      #
####################################################

## Constant variables
# Define paths
#Path for directories of AIMS data
rootpath=/home/stavros/Documents/ABIDE_1/MRC_AIMS/rsfMRI/raw_data/TEST_compilation_script
#path to matlab Hurst file
path_matlab="/home/stavros/Documents/MatlabCode/DWT_Matlab_fractal-master/C_1a_parcEst.m"

#RESAMPLING
#Atlas path and filename to Glasser temp
atlasroot=/home/stavros/rsfmri/Templates/GlasserHCP
atlasfile=$atlasroot/MMP_in_MNI_symmetrical_1.nii.gz

#PREPROCESSING
#specify paths
#path to speedyppX in fmri spt
speedy_path=/home/stavros/rsfmri-master.OLD/speedyppX.py
#path to fd.py in fmri spt
fd_path=/home/stavros/rsfmri-master.OLD/fd.py
#path to dvars_se.py in fmri spt
dvars_path=/home/stavros/rsfmri-master.OLD/dvars_se.py
#path to plot_fd_dvars.py custom MVL script
plot_fd_dvars_path=/home/stavros/fmri_spt/plot_fd_dvars.py


#enter num of CPU's parallel processing env viarable for afni functions
OMP_NUM_THREADS=8
#enter TR
AIMS_TR=1.302000

#Format speedyppX.py's options for preprocessing
FWHM=6
basetime=10
speedyoptions="-f ${FWHM}mm -o --coreg_cfun=lpc+ $basetime --betmask --nobandpass --ss MNI_caez --align_ss --qwarp --rmot --rmotd --keep_means --wds --threshold=10 --SP --OVERWRITE"


##Take subjects
#cd into path
cd $rootpath
#Make a listof the subjects
sublist=$(printf "%s\n" *0* | paste -sd " ")
#%%

##Remove first 8 volumes
for i in $sublist
do
	# define and cd into individual directory
	subpath=$rootpath/$i

	preprocpath2=$subpath/preproc_2
	mkdir $preprocpath2
	cd $preprocpath2
	# Copy mprage and rest in preprocpath2
	cp $subpath/rest.nii.gz $preprocpath2 >> $preprocpath2/${subid}_preproc.sh
	cp $subpath/rest.nii.gz $preprocpath2
	cp $subpath/mprage.nii.gz $preprocpath2 >> $preprocpath2/${subid}_preproc.sh
	cp $subpath/mprage.nii.gz $preprocpath2


	#Specify fslroi variables
	tmin=8
	tsize=612

	#Cut volumes with fslroi
	#output name is Erest.nii.gz
	fslroi rest.nii.gz Erest.nii.gz $tmin $tsize

	#Check output
	NVOL=$(3dinfo -nv Erest.nii.gz)
	ENVOL=$(echo "${NVOL}")

	if [ $ENVOL == 612 ]
		then
		echo "${i} has 612 in Erest.nii.gz"
		echo "${i} has 612 in Erest.nii.gz"
	else
		echo "ERROR in ${i}"
		echo "CUT first vol ERROR in ${i}" >> CUT_first_vol_error_report.txt
	fi
	cd $rootpath
done #for loop ends
#%%

##Preprocessing
cd $rootpath
# Loop over subjects
for subid in $sublist
do
	# Path for specific subject to process
	subpath=$rootpath/$subid
	preprocpath2=$subpath/preproc_2

	# cd into directory
	cd $preprocpath2



	# Make sure -space field in the header is set to ORIG
	echo 3drefit -space ORIG $preprocpath2/mprage.nii.gz >> $preprocpath2/${subid}_preproc.sh
	3drefit -space ORIG $preprocpath2/mprage.nii.gz
	echo 3drefit -space ORIG $preprocpath2/Erest.nii.gz >> $preprocpath2/${subid}_preproc.sh
	3drefit -space ORIG $preprocpath2/Erest.nii.gz
	#Change TR to be sure
	echo 3drefit -TR $AIMS_TR $preprocpath2/Erest.nii.gz >> $preprocpath2/${subid}_preproc.sh
	3drefit -TR $AIMS_TR $preprocpath2/Erest.nii.gz
	# Make sure orientation is LPI
	echo "Checking orientation of rest.nii.gz and reorienting"
	fslorient -getorient Erest.nii.gz
	fslorient -forceneurological Erest.nii.gz
	echo fslorient -forceneurological Erest.nii.gz >> $preprocpath2/${subid}_preproc.sh

	echo "Checking orientation of mprage.nii.gz and reorienting"
	fslorient -getorient mprage.nii.gz
	fslorient -forceneurological mprage.nii.gz
	echo fslorient -forceneurological mprage.nii.gz >> $preprocpath2/${subid}_preproc.sh

	echo "**+ DONE +**"




	# Call speedyppX.py
	echo python $speedy_path -d Erest.nii.gz -a mprage.nii.gz $speedyoptions >> $preprocpath1/${subid}_preproc.sh
	python2  $speedy_path -d Erest.nii.gz -a mprage.nii.gz $speedyoptions
	# Compute framewise displacement with summary statistics
	echo python $fd_path -d Erest_motion.1D >> $preprocpath2/${subid}_preproc.sh
	python $fd_path -d Erest_motion.1D
	# cd into spp.rest
	echo cd spp.Erest >> $preprocpath2/${subid}_preproc.sh
	cd spp.Erest
	# Compute DVARS with summary statistics
	echo python $dvars_path -d Erest_sm.nii.gz >> $preprocpath2/${subid}_preproc.sh
	python $dvars_path -d Erest_sm.nii.gz
	echo python $dvars_path -d Erest_noise.nii.gz >> $preprocpat21/${subid}_preproc.sh
	python $dvars_path -d Erest_noise.nii.gz
	echo python $dvars_path -d Erest_wds.nii.gz >> $preprocpath2/${subid}_preproc.sh
	python $dvars_path -d Erest_wds.nii.gz

	# Run complete preprocessing batch script
	#bash $preprocpath2/${subid}_preproc.sh
	cd $preprocpath2

	# Making PLOTS
	# path to PLOTS

    plotspath=$preprocpath2/PLOTS

	# Make $plotspath
    mkdir $plotspath

    # Change directory to $preprocpath2

    cd $preprocpath2

    # Write commands into bash script

    # Copy paste fd and dvars txt files into $plotspath
    mv Erest_motion_fd.txt $plotspath
    mv spp.Erest/Erest_sm_dvars.txt $plotspath
    mv spp.Erest/Erest_noise_dvars.txt $plotspath
    mv spp.Erest/Erest_wds_dvars.txt $plotspath
    # cd into $plotspath
    cd $plotspath
	# Format arguments for plot_fd_dvars.py
    FD=Erest_motion_fd.txt
    DVARSSM=Erest_sm_dvars.txt
    DVARSNOISE=Erest_noise_dvars.txt
    DVARSWDS=Erest_wds_dvars.txt
    plotpyoptions="--fd $FD --dvars_sm $DVARSSM --dvars_noise $DVARSNOISE --dvars_wds $DVARSWDS"

    # Run plot_fd_dvars.py
    python $plot_fd_dvars_path $plotpyoptions --pdf2save motion_plot.pdf




# cd back to $rootpath
    cd $rootpath
done # for loop ends

#%%
cd $rootpath

##Resampling atlas to native space
# loop over subjects
for subid in $sublist
do
	# Path for each id
	 subpath=$rootpath/$subid
	 preprocpath2=$subpath/preproc_2



	 cd $preprocpath2

	 echo "** In ${subid} Ready to start resampling **"

	 #3dresample is called
	 #output file name is Erest_TEMP.nii.gz
	 3dresample -master Erest_pp.nii.gz -prefix Erest_TEMP.nii.gz -input $atlasfile


	mkdir FINALHURSTS && cp Erest_TEMP.nii.gz FINALHURSTS && ls FINALHURSTS
	cd $rootpath
done # for loop ends
#%%

# Run MATLAB script to parcellate and estimate Hurst exponent
matlab -nodisplay -nosplash -nodesktop -r "run('$path_matlab');"
