function [sort_nav, sort_fys, sort_nav_new, sort_fys_new, F_nav_new, bins_260_fys, bins_260_nav] = ComparingBins(F_nav, F_fys, time_total, gating_signal, params, k)

%% Binning
nBins = params.nBins;
end_time = max(time_total);
binning = params.sortingmethod;
switch binning  
    case 'p1'; p=3;
    case 'p2'; p=4;
    case 'v1'; p=5;
    case 'v2'; p=6;
    case 'v3'; p=7;
end

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

figure(1)
subplot(212);
plot(new_fys_time,new_fys);
title('new fyslog signal');

% Sort both signal in Bins 
[F_nav_new] = SortInBins(new_nav,new_nav_time,params.nBins,1,p);
[F_fys_new] = SortInBins(new_fys,new_fys_time,params.nBins,1,p);

%% sort raw data



%% Figures
figure(1); hold on;
subplot(211); hold on;
E = zeros(nBins,1);
C = hsv(nBins);
    for iii=1:nBins
        B = F_nav_new(:,6)==iii;
        X = F_nav_new(B,1);
        T = new_nav_time(B);
        scatter(T,X,[],C(iii,:),'o');
        axis([0 max(end_time) min(new_nav) max(new_nav)]);
        E(iii,1)=length(find(F_nav_new(:,6)==iii));
    end
    hold off
figure(1); hold on;
subplot(212);
E1 = zeros(nBins,1);
    for ii=1:nBins
        B1 = F_fys_new(:,6)==ii;
        X1 = F_fys_new(B1,1);
        T1 = new_fys_time(B1);
        scatter(T1,X1,[],C(ii,:),'o');
        axis([0 max(end_time) min(new_fys) max(new_fys)]);
        E1(ii,1)=length(find(F_fys_new(:,6)==ii));
    end
    hold off

% Set k_space in Bins
for t=1:nBins
    all = find(F_fys_new(:,6)==t);
    none = find(F_nav_new(:,6)==t);
    for a=1:length(all)
    sort_nav_new(:,a,t)= k(:,a);
    end
    for b=1:length(none)
    sort_fys_new(:,b,t)=k(:,b);
    end
    figure(2); title('FysLog new');
    subplot(ceil(nBins/2),3,t);
    plot(sort_nav_new(:,:,t));
    figure(3); title('Navigator new');
    subplot(ceil(nBins/2),3,t);
    plot(sort_fys_new(:,:,t));
end

%% Sort k-space in Bins by FysLog
bins = nan(1,260,71,8);
bins_260_fys = nan(260,1);
bins_260_nav = nan(260,1);
bins_compare = nan(260,1);
for i=1:260
    bins(:,i,:,:) = F_nav(I_nav(i),p);
    bins_260_fys(i) = F_fys(I_fys(i),p);
    bins_260_nav(i) = F_nav(I_nav(i),p);
    if bins_260_fys(i) == bins_260_nav(i)
        bins_compare(i) = bins_260_fys(i);
    else
        bins_compare(i) = 0;
    end
    
end

% Sort K trajectory
for t1=1:params.nBins
    all = find(bins_260_fys==t1);
    none = find(bins_260_nav==t1);
    for a1=1:length(all)
    sort_nav(:,a1,t1)= k(:,a1);
    end
    for b1=1:length(none)
    sort_fys(:,b1,t1)=k(:,b1);
    end
    figure(800); 
    title('FysLog');
    subplot(ceil(params.nBins/2),3,t1);
    plot(sort_nav(:,:,t1));
    figure(801);
    title('Navigator');
    subplot(ceil(params.nBins/2),3,t1);
    plot(sort_fys(:,:,t1));
end


% FIGURES
figure(110); hold on;
E = zeros(params.nBins,1);
C = hsv(params.nBins);
subplot(512); hold on;
plot(time_total,gating_signal); hold on;
title('Binning direct from k-space according to navigator')
axis([0 end_time -inf inf]);
for iv=1:params.nBins
    B = bins_260_nav(:) == iv;
    X = gating_signal(B);
    T = time_total(B);
    scatter(T,X,[],C(iv,:),'o'); hold on;
    E(iv,1)=length(find(bins_260_nav==iv));
end
hold off;

figure(110); hold on;
E1 = zeros(params.nBins,1);
C1 = hsv(params.nBins);
subplot(514); hold on;
plot(time_total,gating_signal); hold on;
title('Binning direct from k-space according to FysLog')
axis([0 end_time -inf inf]);
for v=1:params.nBins
    B1 = bins_260_fys(:) == v;
    X1 = gating_signal(B1);
    T1 = time_total(B1);
    scatter(T1,X1,[],C1(v,:),'o'); hold on;
    E1(v,1)=length(find(bins_260_fys==v));
end
hold off;

figure(900); hold on;
Ec = zeros(params.nBins,1);
Cc = hsv(params.nBins);
plot(time_total,gating_signal); hold on;
    title('Comparing the binning')
    axis([0 end_time -inf inf]);
for c=1:params.nBins
    Bc = bins_compare(:) == c;
    Xc = gating_signal(Bc);
    Tc = time_total(Bc);
    scatter(Tc,Xc,[],Cc(c,:),'o'); hold on;
    Ec(c,1)=length(find(bins_compare==c));
end
N = bins_compare(:) == 0;
scatter(time_total(N), gating_signal(N), [],'ko'); hold on;
end