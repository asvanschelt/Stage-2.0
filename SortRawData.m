%%
clear a angles klines coils sort_raw_data angle select ans ii bin s iii x y
clc

sort_raw_data = nan(520,29,71,8,6);
length_angle = 0;
for i = 1:6
    angle = find(bins_260_fys(:,1)==i);
    select = raw_data_kspace(:,angle,:,:);
    length_angle = length_angle + length(angle);
    sort_raw_data(:,1:length(angle),:,:,i) = squeeze(raw_data_kspace(:,find(bins_260_fys(:,1)==i),:,:));
    
end

[a b c d e] = isnan(sort_raw_data);