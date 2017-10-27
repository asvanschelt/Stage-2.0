function [F, peak_array, min_array] = SortInBins(dcol,t,params,SR,p)

%% Find peaks and minimals in the signal
L = length(dcol);
Prom = max(dcol)*0.35;
minpeakdistance = SR*(60/21);
[peaks,maxlocs, ~, ~] = findpeaks(dcol, 'MinPeakProminence', Prom, 'MinPeakDistance', minpeakdistance);
[mins,minlocs, ~, ~] = findpeaks(-dcol, 'MinPeakProminence', Prom, 'MinPeakDistance', minpeakdistance);
npeaks = length(maxlocs);
nmins = length(minlocs);

peak_array(:,1) = t(maxlocs);
peak_array(:,2) = peaks;
min_array(:,1) = t(minlocs);
min_array(:,2) = -mins;

%% Start defining
F = nan(L,3);
F(:,1) = dcol;
F(:,2) = t;

nBins = params.nBins;

tops = nan((npeaks+nmins),1);
tops(1:npeaks,1) = maxlocs;
tops(npeaks+1:length(tops),1) = minlocs;
sort_tops = sortrows(tops,1); %an array of all peaks and mins
ntops = length(sort_tops);

up_down = nan(L,4);
up_down(:,1) = dcol;
up_down(:,2) = (1:L); %matrix separately defined for value binning

% Set a limit for the height of the peaks and mins to eliminate large
% deviations
median_tops = params.cutoff*median(peaks);
median_mins = params.cutoff*median(-mins);
binsize = ((median_tops-median_mins)/nBins); 

%% Binning methods
% Phases 1
if p==3 
for i=2:npeaks-1
   phase_length_1 = maxlocs(i-1):maxlocs(i);
   phase_signal_1 = dcol(phase_length_1);
   empty_1=nan(length(phase_length_1),3);
   empty_1(:,1)=phase_signal_1;
   empty_1(:,2)=phase_length_1;
   bindex_1 = round(linspace(0.5,nBins+0.5,length(phase_length_1)));
   empty_1_sort = sortrows(empty_1,1);
   empty_1_sort(:,3) = bindex_1;
   empty_1_sortback = sortrows(empty_1_sort,2);
   
   F(phase_length_1,3) = empty_1_sortback(:,3);
end

% Phases 2
elseif p==4
for i=2:ntops
    phase_length = sort_tops(i-1):sort_tops(i)-1;
    phase_signal = dcol(phase_length);
    
    empty=nan(length(phase_length),3);
    empty(:,1)=phase_signal;
    empty(:,2)=phase_length;
    
    if dcol(sort_tops(i))>dcol(sort_tops(i-1)) % phase up
        bindex = round(linspace(0.5,floor(nBins/2)+0.5,length(phase_length)));
        bindex(bindex>nBins/2) = nBins;
        up_down(phase_length,4)=1;
    else % phase down
        bindex = round(linspace(floor(nBins/2)+0.5,nBins+0.5,length(phase_length)));
        bindex(bindex>nBins) = nBins;
        up_down(phase_length,4)=0;
    end
    
    empty_sort = sortrows(empty,1);
    empty_sort(:,3) = bindex;
    empty_sortback = sortrows(empty_sort,2);
    F(phase_length,3)=empty_sortback(:,3);
end

%% Value binning
% Value 1
elseif p==5
V1 = nan(L,5);
V1(:,1) = dcol;
V1(:,2) = t;
bindex_v = round(linspace(0.5,nBins+0.5,L));
sort_A = sortrows(V1,1);
sort_A(:,3) = bindex_v;
V1 = sortrows(sort_A,2);
F(:,3) = V1(:,3);

% Value 2
elseif p==6
    
for i=2:ntops
    phase_length = sort_tops(i-1):sort_tops(i)-1;    
    if dcol(sort_tops(i))>dcol(sort_tops(i-1)) % phase up
        up_down(phase_length,4)=1;
    else up_down(phase_length,4)=0; % phase down
    end
end

up = up_down(up_down(:,4)==1,:);
down = up_down(up_down(:,4)==0,:);
sort_up = sortrows(up,1);
sort_down = sortrows(down,1);

bindex_up = round(linspace(0.5,floor(nBins/2)+0.5,length(up(:,1))));
bindex_down = round(linspace(floor(nBins/2)+0.5,nBins+0.5,length(down(:,1))));
    bindex_up(bindex_up>floor(nBins/2)) = floor(nBins/2);
    bindex_down(bindex_down>nBins) = nBins;
    
sort_up(:,3) = bindex_up;
sort_down(:,3) = bindex_down;

sortback_up = sortrows(sort_up,2);
sortback_down = sortrows(sort_down,2);

A2 = nan(length(sort_up)+length(sort_down),3);
    A2(1:length(sort_up),1) = sortback_up(:,1);
    A2(length(sort_up)+1:end,1) = sortback_down(:,1);
    A2(1:length(sort_up),2) = sortback_up(:,2);
    A2(length(sort_up)+1:end,2) = sortback_down(:,2);
    A2(1:length(sort_up),3) = sortback_up(:,3);
    A2(length(sort_up)+1:end,3) = sortback_down(:,3);
V2 = sortrows(A2,2);
F((V2(1,2):V2(length(V2),2)),3) = V2(:,3);

% Value 3
elseif p==7
Fdcol = gradient(dcol);
A3 = nan(L,4);
A3(:,1) = dcol;
A3(:,2) = (1:L);
A3(:,3) = Fdcol;

up3 = A3(A3(:,3)>=0,:);
down3 = A3(A3(:,3)<0,:);

sort_up3 = sortrows(up3,1);
sort_down3 = sortrows(down3,1);

bindex_up3 = round(linspace(0.5,floor(nBins/2)+0.5,length(up3(:,1))));
bindex_down3 = round(linspace(floor(nBins/2)+0.5,nBins+0.5,length(down3(:,1))));
    bindex_up3(bindex_up3>floor(nBins/2)) = floor(nBins/2);
    bindex_down3(bindex_down3>nBins) = nBins;
    
sort_up3(:,4) = bindex_up3;
sort_down3(:,4) = bindex_down3;

sortback_up3 = sortrows(sort_up3,2);
sortback_down3 = sortrows(sort_down3,2);

A4 = nan(length(sort_up3)+length(sort_down3),3);
    A4(1:length(sort_up3),1) = sortback_up3(:,1);
    A4(length(sort_up3)+1:end,1) = sortback_down3(:,1);
    A4(1:length(sort_up3),2) = sortback_up3(:,2);
    A4(length(sort_up3)+1:end,2) = sortback_down3(:,2);
    A4(1:length(sort_up3),3) = sortback_up3(:,4);
    A4(length(sort_up3)+1:end,3) = sortback_down3(:,4);
sort_A4 = sortrows(A4,2);
V1(:,5)=sort_A4(:,3);
F(:,3) = V1(:,5);


%% Phaselength binning
% Phaselength 1 
elseif p==8
F_select = F(F(:,1)>=median_mins&F(:,1)<=median_tops,1);
F(F(:,1)>=median_mins&F(:,1)<=median_tops,3)=ceil((F_select-median_mins)/binsize);
F(F(:,1)<median_mins|F(:,1)>median_tops,3)=0;

elseif p==9
for ii=2:ntops
    range = sort_tops(ii-1):sort_tops(ii);
    if F(sort_tops(ii),1)>F(sort_tops(ii-1),1) % Up
        binsize = (median_tops-median_mins)/(floor(nBins/2)); 
        F(range,3) = ceil((F(range,1)-median_mins)/binsize);
    else % Down
        binsize = (median_tops-median_mins)/(ceil(nBins/2)); 
        F(range,3) = ceil((F(range,1)-median_mins)/binsize)+floor(nBins/2);
    end
end
F(F(:,1)<median_mins|F(:,1)>median_tops,3)=0;

end

for l=1:length(F(:,2))
    if isnan(F(l,3))
        if l<(length(F(:,3))/2)
            [c index] = min(abs(F(sort_tops(2):sort_tops(3),1)-F(l,1)));
            F(l,3)=F(index+sort_tops(2),3);
        else
            [c index] = min(abs(F(sort_tops(ntops-2):sort_tops(ntops-1),1)-F(l,1)));
            F(l,3)=F(index+sort_tops(ntops-2),3);
        end
    end
end
end