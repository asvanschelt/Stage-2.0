function [F_nav_new, F_fys_new, bins_260_fys, bins_260_nav] = ComparingBins(F_nav, F_fys, F_MRI, gating_signal, params, k)

%% Binning
nBins = params.nBins;
time_total = F_MRI(:,2);
end_time = max(time_total);
binning = params.sortingmethod;
switch binning  
    case 'p1'; p=3;
    case 'p2'; p=4;
    case 'v1'; p=5;
    case 'v2'; p=6;
    case 'v3'; p=7;
    case 'b1'; p=8;
    case 'b2'; p=9;
end

% Find the closest time value in the Nav or Fys signal which corresponds to
% the MRI signal
edges_nav = [-Inf, mean([F_nav(2:end,2) F_nav(1:end-1,2)],2)', +Inf];
I_nav = discretize(time_total, edges_nav); 
edges_fys = [-Inf, mean([F_fys(2:end,2) F_fys(1:end-1,2)],2)', +Inf];
I_fys = discretize(time_total, edges_fys);

%% New FysLog and Nav signal with undersampling
new_nav = F_nav(I_nav,1);
new_nav_time = F_nav(I_nav,2);
new_fys = F_fys(I_fys,1);
new_fys_time = F_fys(I_fys,2);

figure(1);
subplot(211);
    plot(new_nav_time, new_nav);
    title('new navigator signal');
subplot(212);
    plot(new_fys_time,new_fys);
    title('new fyslog signal');

% Sort both signal in Bins 
[F_nav_new] = SortInBins(new_nav,new_nav_time,params,1,p);
[F_fys_new] = SortInBins(new_fys,new_fys_time,params,1,p);

%% Sort k-space in Bins by FysLog
bins_260_fys = nan(260,1);
bins_260_nav = nan(260,1);
bins_compare = nan(260,1);
for i=1:260
    bins_260_fys(i) = F_fys(I_fys(i),3);
    bins_260_nav(i) = F_nav(I_nav(i),3);
    if bins_260_fys(i) == bins_260_nav(i)
        bins_compare(i) = bins_260_fys(i);
    else
        bins_compare(i) = 0;
    end
end

% Trajectory plotting per Bin
figure(10); hold on;
for b=1:nBins
    subplot(ceil((nBins/3)*2),3,b); hold on;
    trajectories = find(bins_260_fys(:)==b);
    plot(k(:,trajectories));
    title('FysLog');
    subplot(ceil((nBins/3)*2),3,nBins+b); hold on;
    trajectories_nav = find(bins_260_nav(:)==b);
    plot(k(:,trajectories_nav));
    title('Navigator');
end
hold off;

%% Figures
% The Nav and Fys rebinned with undersampling according to the MRI signal
figure(1); hold on;
subplot(211); hold on;
E = zeros(nBins,1);
E1 = zeros(nBins,1);
C = hsv(nBins);
    for ii=1:nBins
        B = F_nav_new(:,3)==ii;
        X = F_nav_new(B,1);
        T = new_nav_time(B);
        subplot(211); hold on;
        scatter(T,X,[],C(ii,:),'o');
        axis([0 max(end_time) min(new_nav) max(new_nav)]);
        E(ii,1)=length(find(F_nav_new(:,3)==ii));
        
        B1 = F_fys_new(:,3)==ii;
        X1 = F_fys_new(B1,1);
        T1 = new_fys_time(B1);
        subplot(212);
        scatter(T1,X1,[],C(ii,:),'o');
        axis([0 max(end_time) min(new_fys) max(new_fys)]);
        E1(ii,1)=length(find(F_fys_new(:,3)==ii));
    end
    hold off

% FIGURE to plot the MRI signal binned to the different signals 
figure(110); hold on;
E = zeros(params.nBins,1);
E1 = zeros(params.nBins,1);
C = hsv(params.nBins);
subplot(512); hold on;
    plot(time_total,gating_signal); hold on;
    title('Binning direct from k-space according to navigator')
    axis([0 end_time -inf inf]);
subplot(514); hold on;
    plot(time_total,gating_signal); hold on;
    title('Binning direct from k-space according to FysLog')
    axis([0 end_time -inf inf]);
for iii=1:params.nBins
        B = bins_260_nav(:) == iii;
        X = gating_signal(B);
        T = time_total(B);
        subplot(512); hold on;
        scatter(T,X,[],C(iii,:),'o'); hold on;
        E(iii,1)=length(find(bins_260_nav==iii));
    subplot(514); hold on;
        B1 = bins_260_fys(:) == iii;
        X1 = gating_signal(B1);
        T1 = time_total(B1);
        scatter(T1,X1,[],C(iii,:),'o'); hold on;
        E1(iii,1)=length(find(bins_260_fys==iii));
end
hold off;
end