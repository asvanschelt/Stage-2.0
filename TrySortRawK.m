function [store] = TrySortRawK(bins, raw_data_kspace, params)

nBins = params.nBins;
arrange(:,1) = [1:260];
arrange(:,2) = bins;
sort_arrange = sortrows(arrange,2);

% Find the size of the bins to later select all k-lines in that bin
size_bins = nan(nBins,1);
for i =1:nBins
    size_bins(i) = length(find(arrange(:,2)==i));
end

start_place = length(find(arrange(:,2)==0));

% Sort the data from bin 1 to nBins same as above
raw_data_sort = raw_data_kspace(:,sort_arrange(:,1),:,:);
start = 1 + start_place;
for ii = 1:nBins
    eind = start + size_bins(ii) - 1;
    store{ii} = raw_data_sort(:,start:eind,:,:);
    start = start+size_bins(ii);
end
if start_place > 0 
    store{nBins+1}= raw_data_sort(:,1:start_place,:,:);
end
end