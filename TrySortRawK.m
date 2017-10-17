arrange(:,1) = [1:260];
arrange(:,2) = bins_260_fys;

sort_arrange = sortrows(arrange,2);

size_bins = nan(6,1);
for i =1:6
    size_bins(i) = length(find(arrange(:,2)==i));
end

raw_data_sort = raw_data_kspace(:,sort_arrange(:,1),:,:);

start = 1;
eind = size_bins(1);
for ii = 1:6
    eind = start + size_bins(ii) -1 ;
    store{ii} = raw_data_sort(:,start:eind,:,:);
    start = start+size_bins(ii);
end