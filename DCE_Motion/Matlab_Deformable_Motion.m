addpath /home/abeers/Projects/Matlab_DCE/Matlab_QTIM_Script/NIfTI_20140122/

input_file='/home/abeers/Projects/DCE_Motion_Phantom/DCE_MRI_Phantom_Regenerated_Signal.nii.gz'

input_data=load_untouch_nii(input_file);

deformation_matrix = deformation_matrix / 2;
final_deformation_matrix = deformation_matrix;

for t = 1:size(input_data.img,4)
    input_data.img(:,:,:,t)=imwarp(input_data.img(:,:,:,t), final_deformation_matrix(:,:,:,:,t));
end

save_untouch_nii(input_data, '/home/abeers/Projects/DCE_Motion_Phantom/DCE_MRI_Phantom_Regenerated_Signal_Weak_Motion.nii.gz')

% input_data.hdr.dime.dim(1)=4;
% input_data.hdr.dime.dim(5)=10;
% input_data.hdr.dime.pixdim(1)=1;
input_data.hdr.dime.datatype=16;
input_data.hdr.dime.bitpix=32;
input_data.hdr.dime.cal_max=0;
input_data.hdr.dime.glmax=0;     

for vector = 1:3
	input_data.img = final_deformation_matrix(:,:,:,vector,1);
	save_untouch_nii(input_data, strcat('/home/abeers/Projects/DCE_Motion_Phantom/Deformation_Vectors_Weak', num2str(vector),'.nii.gz'));