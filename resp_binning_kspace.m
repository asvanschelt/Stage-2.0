clear all; close all; clc; 

%PARAMETERS TO CHANGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% addpath(genpath('/home/jmdenharder/lood_storage/divi/Ima/parrec/Chiel/DCE_Knie/Reconstruction_code/DCE_code/'))
addpath(genpath('/home/pdeheer/lood_storage/divi/Users/pdeheer/Code/Github/dhrpaul/R4D/'))

load('/home/pdeheer/lood_storage/divi/Ima/parrec/Pauldh/2017-10-03_Liver_radial/Raw_data/Recon/matlab_variables_binning.mat');
[kdatau,ku] = ksp2frames(raw_data_kspace,trajectory_kspace,P.binparams);
