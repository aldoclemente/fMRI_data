#!/bin/bash

user="$(whoami)"

if ! [ -d "/media/${user}/EXTERNAL_USB/data/thickness_preprocessed/" ]; then
  mkdir "/media/${user}/EXTERNAL_USB/data/thickness_preprocessed/"
fi

if ! [ -d "/media/${user}/EXTERNAL_USB/data/thickness_processed/" ]; then
  mkdir "/media/${user}/EXTERNAL_USB/data/thickness_preprocessed/"
fi

if ! [ -d "/media/${user}/EXTERNAL_USB/data/rfMRI_REST_preprocessed/" ]; then
  mkdir "/media/${user}/EXTERNAL_USB/data/rfMRI_REST_preprocessed/"
fi

if ! [ -d "/media/${user}/EXTERNAL_USB/data/rfMRI_REST_processed/" ]; then
  mkdir "/media/${user}/EXTERNAL_USB/data/rfMRI_REST_processed/"
fi

for filename in /media/${user}/EXTERNAL_USB/data/thickness/*.nii; do
    wb_command -cifti-separate "$filename" COLUMN -metric CORTEX_LEFT "/media/${user}/EXTERNAL_USB/data/thickness_preprocessed/$(basename "$filename" .nii).gii"
done

for filename in /media/${user}/EXTERNAL_USB/data/thickness_preprocessed/*.gii; do
    gifti_tool -infile "$filename" -write_1D "/media/${user}/EXTERNAL_USB/data/thickness_processed/$(basename "$filename" .gii).1D"
done

for filename in /media/${user}/EXTERNAL_USB/data/rfMRI_REST/*.nii; do
    wb_command -cifti-separate "$filename" COLUMN -metric CORTEX_LEFT "/media/${user}/EXTERNAL_USB/data/rfMRI_REST_preprocessed/$(basename "$filename" .nii).gii"
done

for filename in /media/${user}/EXTERNAL_USB/data/rfMRI_REST_preprocessed/*.gii; do
    gifti_tool -infile "$filename" -write_1D "/media/${user}/EXTERNAL_USB/data/rfMRI_REST_processed/$(basename "$filename" .gii).1D"
done
