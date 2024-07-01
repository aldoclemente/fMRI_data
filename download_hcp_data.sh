#!/bin/bash

user="$(whoami)"
m=80
head -$m subject_idxs.txt > input_idxs.txt

if ! [ -d "/media/${user}/EXTERNAL_USB/data/thickness/" ]; then
    mkdir "/media/${user}/EXTERNAL_USB/data/thickness/"
fi

while read subject; do
    file=/media/${user}/EXTERNAL_USB/data/thickness/${subject}.thickness.32k_fs_LR.dscalar.nii
    if ! [ -f $file ]; then
        aws s3 cp s3://hcp-openaccess/HCP/${subject}/MNINonLinear/fsaverage_LR32k/${subject}.thickness.32k_fs_LR.dscalar.nii $file
    fi
done < input_idxs.txt 
  
if ! [ -d "/media/${user}/EXTERNAL_USB/data/rfMRI_REST/" ]; then
    mkdir "/media/${user}/EXTERNAL_USB/data/rfMRI_REST/"
fi  

while read subject; do
    file=/media/${user}/EXTERNAL_USB/data/rfMRI_REST/${subject}.rfMRI_REST1_LR_Atlas.dtseries.nii
    if ! [ -f $file ]; then
    aws s3 cp s3://hcp-openaccess/HCP/${subject}/MNINonLinear/Results/rfMRI_REST1_LR/rfMRI_REST1_LR_Atlas.dtseries.nii $file
    fi
done < input_idxs.txt

if ! [ -d /media/${user}/EXTERNAL_USB/data/DTI/ ]; then
  mkdir /media/${user}/EXTERNAL_USB/data/DTI/
fi

while read subject; do
    if ! [ -d /media/${user}/EXTERNAL_USB/data/DTI/${subject} ]; then
        mkdir /media/${user}/EXTERNAL_USB/data/DTI/${subject}
        aws s3 sync s3://hcp-openaccess/HCP/$subject/T1w/Diffusion /media/${user}/EXTERNAL_USB/data/DTI/${subject}
    fi
done < input_idxs.txt



