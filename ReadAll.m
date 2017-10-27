function [store_MRI, store_Fys, store_Nav, store_Fys_d, store_Nav_d] = ReadAll(signal1, signal2, gating_signal_s, params, gating_signal, k, raw_data_kspace)

%% PARAMETERS
% signal1 = 'SCANPHYSLOG_li_03102017_1804326_14_1_dce_spair_70ms_radial_navV4.log';
% signal2 = 'devlogcurrent_1804-1807.log';
% binning = 'v2';
% nBins = 6;

%% START CODE
nBins = params.nBins;
binning = params.sortingmethod;
time_total = [0:2*params.time_blok:2*params.time_blok*260-params.time_blok];

[dcol_fys, t_fys] = ReadFysLog(signal1);
[dcol_nav, t_nav] = ReadNav(signal2);

% Equalize the signals
dcol_nav(:) = dcol_nav(:)-mean(dcol_nav); %remove drift
dcol_fys(:) = dcol_fys(:)/(max(dcol_fys)/max(dcol_nav));
gating_signal_s(:) = gating_signal_s(:)/(max(gating_signal_s)/max(dcol_nav));

% Align all signals
if t_nav(length(t_nav))>t_fys(length(t_fys))
    diff = t_nav(length(t_nav))-t_fys(length(t_fys));
    t_fys(:) = t_fys(:)+diff;
elseif t_nav(length(t_nav))<=t_fys(length(t_fys))
    diff = t_fys(length(t_fys))-t_nav(length(t_nav));
    t_nav(:) = t_nav(:)+diff;
end

time_total(:) = time_total(:)-(time_total(length(time_total))-t_nav(length(t_nav)));

switch binning  
    case 'p1'; p=3;
    case 'p2'; p=4;
    case 'v1'; p=5;
    case 'v2'; p=6;
    case 'v3'; p=7;
    case 'b1'; p=8;
    case 'b2'; p=9;
end

% Sort each signal in Bins
[F_nav, peaks_nav, mins_nav] = SortInBins(dcol_nav,t_nav,params,1,p);
[F_fys, peaks_fys, mins_fys] = SortInBins(dcol_fys,t_fys,params,500,p);
[F_MRI, peaks_MRI, mns_MRI] = SortInBins(gating_signal_s,time_total,params,1,p);

% Translate both the fys and nav to MRI signal for reconstruction
[F_nav_new, F_fys_new, bins_260_fys, bins_260_nav] = ComparingBins(F_nav, F_fys, F_MRI, gating_signal, params, k);

% Sort the raw data in the bins defined by the different methods
[store_MRI] = TrySortRawK(F_MRI(:,3), raw_data_kspace, params);
[store_Nav] = TrySortRawK(F_nav_new(:,3), raw_data_kspace, params);
[store_Fys] = TrySortRawK(F_fys_new(:,3), raw_data_kspace, params);
[store_Fys_d] = TrySortRawK(bins_260_fys, raw_data_kspace, params);
[store_Nav_d] = TrySortRawK(bins_260_nav, raw_data_kspace, params);


%% Figures
% Figure with the different binning methods in subplots
figure(110); hold on;
E = zeros(nBins,1);
E1 = zeros(nBins,1);
E3 = zeros(nBins,1);
C = hsv(nBins);
subplot(511);
    plot(t_nav, dcol_nav); hold on;
    title('Navigator Binning')
subplot(513);
    plot(t_fys, dcol_fys); hold on;
    title('FysLog Binning')
subplot(515);hold on;
    plot(time_total, gating_signal_s); hold on;
    axis([0 max(time_total) -inf inf]);
    title('MRI Binning')
    for i=1:nBins
        subplot(511);
        B = F_nav(:,3)==i;
        X = F_nav(B,1);
        T = t_nav(B);
        scatter(T,X,[],C(i,:),'o');
        axis([0 max(t_nav) min(dcol_nav) max(dcol_nav)]);
        E(i,1)=length(find(F_nav(:,3)==i));
        
        subplot(513);
        B1 = F_fys(:,3)==i;
        X1 = F_fys(B1,1);
        T1 = t_fys(B1);
        scatter(T1,X1,[],C(i,:),'.');
        axis([0 max(t_fys) min(dcol_fys) max(dcol_fys)]);
        E1(i,1)=length(find(F_nav(:,3)==i));
        
        subplot(515);
        B3 = F_MRI(:,3)==i;
        X3 = F_MRI(B3,1);
        T3 = time_total(B3);
        scatter(T3,X3,[],C(i,:),'o');
        E3(i,1)=length(find(F_MRI(:,3)==i));
    end
    hold off


% Figure with all 3 different signals in one with the peaks and mins
figure(150); hold on;
plot(t_nav,dcol_nav,'r'); hold on;
scatter(peaks_nav(:,1),peaks_nav(:,2), 'r*'); hold on;
scatter(mins_nav(:,1),mins_nav(:,2), 'r*'); hold on;
plot(t_fys,dcol_fys,'b'); hold on;
scatter(peaks_fys(:,1),peaks_fys(:,2), 'b*'); hold on;
scatter(mins_fys(:,1),mins_fys(:,2), 'b*'); hold on;
plot(time_total,F_MRI(:,1),'g'); hold on;
scatter(peaks_MRI(:,1),peaks_MRI(:,2), 'g*'); hold on;
scatter(mns_MRI(:,1),mns_MRI(:,2), 'g*'); hold off;
    legend('Navigator','Peaks navigator','','Fyslog', 'Peaks fyslog','','MRI', 'Peaks MRI','');
end