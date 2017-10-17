function [dcol_fys, t_fys, t_mark1, t_mark2, t_start] = ReadFysLog(signal)

%% Read the scan of signal 1
D = ReadPhilipsScanPhysLog(signal);

% Plot the respiratory motion
L = length(D.C);
SR = 500;
t_fys = (1:L)*(1/SR);
col = D.C(1:L,6);
dcol_fys_1 = double(col);

fileID= fopen(signal); 
text = textscan(fileID,'%s');
fclose(fileID); 

search_markers = 'mark2';
search_start = '0010';
search_end = '0020';
search_mark1 = '0080';
search_mark2 = '0100';

ntitel = size(D.C,2);
nmark = 10;

loc_markers= find(strcmp(text{1},search_markers));
loc_start= find(strcmp(text{1},search_start));
loc_end_all= find(strcmp(text{1},search_end));
loc_mark1= find(strcmp(text{1},search_mark1));
loc_mark2= find(strcmp(text{1},search_mark2));

loc_start = ((loc_start-loc_markers)-nmark) / ntitel;
loc_end = ((loc_end_all(:)-loc_markers)-nmark) / ntitel;
loc_mark1 = ((loc_mark1(:)-loc_markers)-nmark) / ntitel;
loc_mark2 = ((loc_mark2(:)-loc_markers)-nmark) / ntitel;

t_mark1 = t_fys(loc_mark1);
t_mark2 = t_fys(loc_mark2);
t_start = t_fys(loc_start);

dcol_fys = dcol_fys_1(loc_start:loc_end(length(loc_end)));
dcol_fys = -dcol_fys;
t_fys = t_fys(1:length(dcol_fys)).';
end
