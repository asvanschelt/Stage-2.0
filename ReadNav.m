function [dcol_nav, t_nav] = ReadNav(signal)

%% Log File open
fileID= fopen(signal); % open the file in daataloc with name dataname. returns 3 when opened successfully
text = textscan(fileID,'%s'); % store content of currently opened file (fileID) as text ; '%s' as character vector
fclose(fileID); % close open file

search_nav= 'MMURNAV(debug_pos):'; % correspondance expression in txt file to find the time stamp for navigator 
search_end= 'RNAV:';
timestring= 7; % number of strings/cells for which you need to move to find the time stamp in the text vector
navstring= 7;

loc_nav= find(strcmp(text{1},search_nav));
loc_end= find(strcmp(text{1},search_end));

% remove initiation
timecell_end = text{1}(loc_end(length(loc_end))-timestring);
time_pos_end = char([timecell_end]');
time_pos_end= [time_pos_end(1,1:2), time_pos_end(1,4:5), time_pos_end(1,7:end)];
time_end(1) = str2double(time_pos_end(1,:));
h_end= fix(time_end/10000); % rounds to nearest integer 
tot_min_end= (time_end/10000-h_end)*100;
m_end= fix(tot_min_end);
s_end= (tot_min_end-m_end)*100;
time_end_s= h_end*60*60+ m_end*60+s_end;

for i = 1:length(loc_nav)
    timecell{i}= text{1}(loc_nav(i)-timestring); % go to "text" and search for 'search_nav', go 7 element back 
    navcell{i}=text{1}(loc_nav(i)+navstring);
end

time_pos = char([timecell{:}]');
nav_pos = char([navcell{:}]');

% Erase the double point ':' between time stamps
time_pos= [time_pos(:,1:2), time_pos(:,4:5), time_pos(:,7:end)]; 

for i = 1:length(loc_nav)
   navigator(i) = str2double(nav_pos(i,:)); % each number in one element
   time(i)= str2double(time_pos(i,:));
end

% Change dynamic to seconds and start from 0)
h= fix(time/10000); % rounds to nearest integer 
tot_min= (time/10000-h)*100;
m= fix(tot_min);
s= (tot_min-m)*100;

%Counter rounding issues
cor_m= s>=60; 
s(cor_m)=0;
m=m+cor_m;
time_s= h*60*60+ m*60+s;

% Calculate average frequency of navigator
f_nav= length(navigator)/max(time_s);

%Correction for drift in the set of, set set of around 0
nav_offset= (max(navigator)+min(navigator))/2;
navigator= navigator - nav_offset;
time_offset= time_s(1);
time_s= time_s - time_offset; % start at 0 sec 
time_s_end = time_end_s - time_offset;


%% Signal processing
dcol_nav = navigator.';
t_nav_old = time_s.';
t_nav_end = t_nav_old(length(t_nav_old));

% for tt=1:length(time_s)
%     if t_nav(tt,1) < time_s_end
%         dcol_nav(tt,1) = nan;
%     end
% end

dcol_nav(isnan(dcol_nav(:,1)))=[];
%t_nav = t_nav(1:length(dcol_nav));

steps = t_nav_end / length(t_nav_old);

t_nav = [steps:steps:(steps*length(dcol_nav))].';
[~,L] = size(t_nav);

end