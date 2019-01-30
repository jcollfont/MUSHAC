#! /bin/bash
# $1 is the path

# source python
source ~/anaconda2_startup.sh

crlConvertImageToVectorImage $1/connectom/st/dwi_modb_x.nii.gz $1//connectom/st/dwi_modb_x_vec.nii
crlResampler2 --voxelsize 1.0,1.0,1.0 -i $1/connectom/st//dwi_modb_x_vec.nii -o $1/connectom/st//dwi_modb_x_1mm.nii --interp sinc -p 40

crlConvertImageToVectorImage $1/connectom/st//dwi_modb_y.nii.gz $1//connectom/st/dwi_modb_y_vec.nii
crlResampler2 --voxelsize 1.0,1.0,1.0 -i $1/connectom/st//dwi_modb_y_vec.nii -o $1//connectom/st/dwi_modb_y_1mm.nii --interp sinc -p 40

crlConvertImageToVectorImage $1//connectom/st/dwi_modb_z.nii.gz $1//connectom/st/dwi_modb_z_vec.nii
crlResampler2 --voxelsize 1.0,1.0,1.0 -i $1//connectom/st/dwi_modb_z_vec.nii -o $1//connectom/st/dwi_modb_z_1mm.nii --interp sinc -p 40



crlConvertImageToVectorImage $1//connectom/sa/dwi_modb_x.nii.gz $1//connectom/sa/dwi_modb_x_vec.nii
crlResampler2 --voxelsize 1.0,1.0,1.0 -i $1//connectom/sa/dwi_modb_x_vec.nii -o $1//connectom/sa/dwi_modb_x_1mm.nii --interp sinc -p 40

crlConvertImageToVectorImage $1//connectom/sa/dwi_modb_y.nii.gz $1//connectom/sa/dwi_modb_y_vec.nii
crlResampler2 --voxelsize 1.0,1.0,1.0 -i $1//connectom/sa/dwi_modb_y_vec.nii -o $1//connectom/sa/dwi_modb_y_1mm.nii --interp sinc -p 40

crlConvertImageToVectorImage $1//connectom/sa/dwi_modb_z.nii.gz $1//connectom/sa/dwi_modb_z_vec.nii
crlResampler2 --voxelsize 1.0,1.0,1.0 -i $1//connectom/sa/dwi_modb_z_vec.nii -o $1//connectom/sa/dwi_modb_z_1mm.nii --interp sinc -p 40