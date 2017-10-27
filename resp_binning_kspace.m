clear all; close all; clc; 

%PARAMETERS TO CHANGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%addpath(genpath('/home/jmdenharder/lood_storage/divi/Ima/parrec/Chiel/DCE_Knie/Reconstruction_code/DCE_code/'))
%addpath(genpath('/home/pdeheer/lood_storage/divi/Users/pdeheer/Code/Github/dhrpaul/R4D/'))

load('matlab_variables_binning.mat');
signal1 = 'SCANPHYSLOG_li_03102017_1804326_14_1_dce_spair_70ms_radial_navV4.log';
signal2 = 'devlogcurrent_1804-1807.log';

% figure(93); montage(abs(permute(raw_data_kspace(:,:,:,1),[1 2 4 3])),[0 200])
% figure(93); montage(abs(permute(ifft(raw_data_kspace(:,:,:,1),[],3),[1 2 4 3])),[0 200])
% image_space = ifft(raw_data_kspace(:,:,:,1),[],3);

P.binparams.nBins = 6;
P.binparams.sortingmethod = 'v1';
P.binparams.time_blok = 0.186;
P.binparams.cutoff = 1.2;

[kdatau,ku, gating_signal, ~, gating_signal_s] = ksp2frames(raw_data_kspace,trajectory_kspace,P.binparams);

[store_MRI, store_Fys, store_Nav, store_Fys_d, store_Nav_d] = ReadAll(signal1, signal2, gating_signal_s, P.binparams, gating_signal, trajectory_kspace, raw_data_kspace);

for i=1:P.binparams.nBins
    MRI(i) = size(store_MRI{i},2);
    Fys(i) = size(store_Fys{i},2);
    Nav(i) = size(store_Nav{i},2);
    Fys_d(i) = size(store_Fys_d{i},2);
    Nav_d(i) = size(store_Nav_d{i},2);
end

