function [F_nav, F_fys, F_MRI] = ReadAll(signal1, signal2, gating_signal_s, params)

% signal1 = 'SCANPHYSLOG_li_03102017_1804326_14_1_dce_spair_70ms_radial_navV4.log';
% signal2 = 'devlogcurrent_1804-1807.log';
% binning = 'v2';
% nBins = 6;

nBins = params.nBins;
binning = params.sortingmethod;
time_total = [0:2*params.time_blok:2*params.time_blok*260-params.time_blok];

[dcol_fys, t_fys] = ReadFysLog(signal1);
[dcol_nav, t_nav] = ReadNav(signal2);

dcol_nav(:) = dcol_nav(:)-mean(dcol_nav);
dcol_fys(:) = dcol_fys(:)/(max(dcol_fys)/max(dcol_nav));

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
end

[F_nav, peaks_nav, mins_nav] = SortInBins(dcol_nav,t_nav,nBins,1,p);
[F_fys, peaks_fys, mins_fys] = SortInBins(dcol_fys,t_fys,nBins,500,p);
[F_MRI, peaks, mns] = SortInBins(gating_signal_s,time_total,nBins,1,p);
        
figure(110); hold on;
E = zeros(nBins,1);
E1 = zeros(nBins,1);
E3 = zeros(nBins,1);
C = hsv(nBins);
subplot(511);
    plot(t_nav, dcol_nav); hold on;
    title('Navigator Binning')
    for i=1:nBins
        B = F_nav(:,p)==i;
        X = F_nav(B,1);
        T = t_nav(B);
        scatter(T,X,[],C(i,:),'o');
        axis([0 max(t_nav) min(dcol_nav) max(dcol_nav)]);
        E(i,1)=length(find(F_nav(:,p)==i));
    end
    hold off
subplot(513);
    plot(t_fys, dcol_fys); hold on;
    title('FysLog Binning')
    for ii=1:nBins
        B1 = F_fys(:,p)==ii;
        X1 = F_fys(B1,1);
        T1 = t_fys(B1);
        scatter(T1,X1,[],C(ii,:),'.');
        axis([0 max(t_fys) min(dcol_fys) max(dcol_fys)]);
        E1(ii,1)=length(find(F_nav(:,p)==ii));
    end
    hold off
subplot(515);hold on;
plot(time_total, gating_signal_s); hold on;
axis([0 max(time_total) -inf inf]);
title('MRI Binning')
for iii=1:nBins
    B3 = F_MRI(:,p)==iii;
    X3 = F_MRI(B3,1);
    T3 = time_total(B3);
    scatter(T3,X3,[],C(iii,:),'o');
    E3(iii,1)=length(find(F_MRI(:,p)==iii));
end
hold off


figure(150); hold on;
plot(t_nav,dcol_nav,'r'); hold on;
scatter(peaks_nav(:,1),peaks_nav(:,2), 'r*'); hold on;
scatter(mins_nav(:,1),mins_nav(:,2), 'r*'); hold on;
plot(t_fys,dcol_fys,'b'); hold on;
scatter(peaks_fys(:,1),peaks_fys(:,2), 'b*'); hold on;
scatter(mins_fys(:,1),mins_fys(:,2), 'b*'); hold on;
plot(time_total,F_MRI(:,1)*100,'g'); hold on;
scatter(peaks(:,1),peaks(:,2)*100, 'g*'); hold on;
scatter(mns(:,1),mns(:,2)*100, 'g*'); hold off;
    legend('Navigator','Peaks navigator','','Fyslog', 'Peaks fyslog','','MRI', 'Peaks MRI','');
end