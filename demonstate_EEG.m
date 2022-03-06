%clear all
%load (eyes open) EEG data from one of 10 subjects - (S02 - S11)_restingPre_EO.mat
data = load('S02_restingPre_EO.mat').dataRest;

%data = load('S_concat_restingPre_EO.mat').concat_data;
%sizeperfile = 38401;
%data(:,7*sizeperfile+1:8*sizeperfile) = [];
%data( = data(:,1:100000);

%channel names are copied from database documentation
channels = ["Fp1" "AF7" 'AF3' 'F1' 'F3' 'F5' 'F7' 'FT7' 'FC5' 'FC3' 'FC1' 'C1' 'C3' 'C5' 'T7' 'TP7' 'CP5' 'CP3' 'CP1' 'P1' 'P3' 'P5' 'P7' 'P9' 'PO7' 'PO3' 'O1' 'Iz' 'Oz' 'POz' 'Pz' 'CPz' 'Fpz' 'Fp2' 'AF8' 'AF4' 'Afz' 'Fz' 'F2' 'F4' 'F6' 'F8' 'FT8' 'FC6' 'FC4' 'FC2' 'FCz' 'Cz' 'C2' 'C4' 'C6' 'T8' 'TP8' 'CP6' 'CP4' 'CP2' 'P2' 'P4' 'P6' 'P8' 'P10' 'PO8' 'PO4' 'O2'];
left_ch = ["Fp1" "AF3" "AF7" 'F1' 'F3' 'F5' 'F7' 'F9' 'FC1' 'FC3' 'FC5' 'FT7' 'FT9' 'C1' 'C3' 'C5' 'T7' 'T9' 'CP1' 'CP3' 'CP5' 'TP7' 'TP9' 'P1' 'P3' 'P5' 'P7' 'P9' 'PO3' 'PO7' 'O1'];
[left_ch left_ia left_ib]= intersect(left_ch,channels,'stable');
right_ch = ["Fp2" 'AF4' 'AF8' 'F2' 'F4' 'F6' 'F8' 'F10' 'FC2' 'FC4' 'FC6' 'FT8' 'FT10' 'C2' 'C4' 'C6' 'T8' 'T10' 'CP2' 'CP4' 'CP6' 'TP8' 'TP10' 'P2' 'P4' 'P6' 'P8' 'P10' 'PO4' 'PO8' 'O2'];
[right_ch right_ia right_ib] = intersect(right_ch,channels,'stable');
central_ch = ["Nz" 'Fpz' 'Afz' 'Fz' 'FCz' 'Cz' 'CPz' 'Pz' 'POz' 'Oz' 'Iz'];
[central_ch central_ia central_ib] = intersect(central_ch,channels,'stable');

data(65:end,:) = [];
frontal = [find(channels == 'Fp1') find(channels == 'Fp2')];
central = [find(channels == 'C1') find(channels == 'C2')];
occipital = [find(channels == 'O1') find(channels == 'O2')];
figure('Name','Central electrode')
plot(data(central(1),:)')
figure('Name','Frontal electrode')
plot(data(frontal(1),:)')
figure('Name','Occipital electrode')
plot(data(occipital(2),:)')

data_frontal = data(frontal(1),:); %data_frontal is a frontal EEG electrode timeseries
%% find peaks of frontal electrode

%After finding the peaks (all local maxima) and computing their peak to peak distance,
%we implement a threshold so that only eyeblink spikes are left
%the threshold is implemented by the second argument of find_peak_ts
%(eyeblink_dev). find_peaks_ts will save all spikes with a peak to peak 
%distance value of at least eyeblink_dev multiplied by the maximum mean peak 
%to peak distance calculated across all local maxima.

eyeblink_dev = 0.5;
[distances indexes] = find_peaks_ts(data_frontal,eyeblink_dev);

i = 2;
%The following while loop will remove peaks that are too close to each
%other (less than 50 samples apart). These peaks likely correspond to the
%same eyeblink
while i <= size(indexes,2)
    if indexes(i)-indexes(i-1) < 50
        if distances(i) > distances(i-1)
            distances(i-1) = [];
            indexes(i-1) = [];
            i = i-1;
        else
            distances(i) = [];
            indexes(i) = [];
            i = i-1;
        end
    end
    i = i+1;
end
%% find peak to peak distance for all electrodes
%This time, the peak to peak distance is not calculated from the mean local
%maximum to the mean local minimum. Instead we compute the distance from the
%mean local maxima to the value of the timeseries 35 samples earlier
%(window variable). This variable is assumed to be the average number of
%samples between the beginning and maximum value of an eyeblink spike. 
%Since occipital and central channels exhibit no noticeable
%spikes during an eyeblink, using the local minimum to compute the peak
%to peak distance would not work. The window_mean variable is the number of
%samples around which we compute a mean (to lower the effect of noise).
%The distance calculated for each channel and each spike/eyeblink is stored
%in the matrix dist.
%The average distance along all spikes/eyeblinks for each of the 68
%electrodes is stored in the vector mean_ch_dist.

%%%%%%%%%%%%%%%%%%%%%%%
window_mean = 3;  %%%%% parameters
window = 35;      %%%%%
%%%%%%%%%%%%%%%%%%%%%%%

dist = zeros(size(data,1),size(indexes,2));
for ch=1:size(data,1)
    for i = 1:size(indexes,2)
        dist(ch,i) = mean(data(ch,indexes(i)-window_mean:indexes(i)+window_mean)) - mean(data(ch,indexes(i)-window_mean-window:indexes(i)-window+window_mean));
    end
end
mean_ch_dist = mean(dist');


%% Save projection matrix

%vectors containing the mean distance computed in the previous segment for
%left side electrodes, right side electrodes and central electrodes. The
%first elements in the vector correspond to the frontal electrodes and
%later elements correspond to occipital electrodes.
%The assymetry vector contains the difference in mean channel distance
%between each left side electrode and the corresponding right side
%electrode.
left_mean_dist = mean_ch_dist(left_ib);
right_mean_dist = mean_ch_dist(right_ib);
central_mean_dist = mean_ch_dist(central_ib);
assymetry = left_mean_dist - right_mean_dist;
figure('Name','Left brain side')
plot(left_mean_dist)
figure('Name','Right brain side')
plot(right_mean_dist)
figure('Name','Central brain side')
plot(central_mean_dist)
figure('Name','Left/Right assymetry')
plot(assymetry);

%occipital channels tend to have slightly negative values

%save('eyeblink_projection.mat','channels','mean_ch_dist');
%% Concatenate all data and repeat the code above

clear all

concat_data = load('S02_restingPre_EO.mat').dataRest;
next_data = load('S02_restingPre_EO.mat').dataRest;
difference = concat_data(:,end) - next_data(:,1);
next_data = next_data + difference;
next_data(:,1) = [];
i = 2;
while i < 12
    concat_data = [concat_data next_data];
    next_sub = strcat('S0',num2str(i),'_restingPre_EO.mat');
    if i > 9
        next_sub(2) = [];
    end
    next_data = load(next_sub).dataRest;
    difference = concat_data(:,end) - next_data(:,1);
    next_data = next_data + difference;
    next_data(:,1) = [];
    i = i+1;
end

data = concat_data;

channels = ["Fp1" "AF7" 'AF3' 'F1' 'F3' 'F5' 'F7' 'FT7' 'FC5' 'FC3' 'FC1' 'C1' 'C3' 'C5' 'T7' 'TP7' 'CP5' 'CP3' 'CP1' 'P1' 'P3' 'P5' 'P7' 'P9' 'PO7' 'PO3' 'O1' 'Iz' 'Oz' 'POz' 'Pz' 'CPz' 'Fpz' 'Fp2' 'AF8' 'AF4' 'Afz' 'Fz' 'F2' 'F4' 'F6' 'F8' 'FT8' 'FC6' 'FC4' 'FC2' 'FCz' 'Cz' 'C2' 'C4' 'C6' 'T8' 'TP8' 'CP6' 'CP4' 'CP2' 'P2' 'P4' 'P6' 'P8' 'P10' 'PO8' 'PO4' 'O2'];
left_ch = ["Fp1" "AF3" "AF7" 'F1' 'F3' 'F5' 'F7' 'F9' 'FC1' 'FC3' 'FC5' 'FT7' 'FT9' 'C1' 'C3' 'C5' 'T7' 'T9' 'CP1' 'CP3' 'CP5' 'TP7' 'TP9' 'P1' 'P3' 'P5' 'P7' 'P9' 'PO3' 'PO7' 'O1'];
[left_ch left_ia left_ib]= intersect(left_ch,channels,'stable');
right_ch = ["Fp2" 'AF4' 'AF8' 'F2' 'F4' 'F6' 'F8' 'F10' 'FC2' 'FC4' 'FC6' 'FT8' 'FT10' 'C2' 'C4' 'C6' 'T8' 'T10' 'CP2' 'CP4' 'CP6' 'TP8' 'TP10' 'P2' 'P4' 'P6' 'P8' 'P10' 'PO4' 'PO8' 'O2'];
[right_ch right_ia right_ib] = intersect(right_ch,channels,'stable');
central_ch = ["Nz" 'Fpz' 'Afz' 'Fz' 'FCz' 'Cz' 'CPz' 'Pz' 'POz' 'Oz' 'Iz'];
[central_ch central_ia central_ib] = intersect(central_ch,channels,'stable');
frontal = [find(channels == 'Fp1') find(channels == 'Fp2')];
central = [find(channels == 'C1') find(channels == 'C2')];
occipital = [find(channels == 'O1') find(channels == 'O2')];

data_frontal = data(frontal(1),:); %data_frontal is a frontal EEG electrode timeseries


eyeblink_dev = 0.5;
[distances indexes] = find_peaks_ts(data_frontal,eyeblink_dev);

i = 2;

while i <= size(indexes,2)
    if indexes(i)-indexes(i-1) < 50
        if distances(i) > distances(i-1)
            distances(i-1) = [];
            indexes(i-1) = [];
            i = i-1;
        else
            distances(i) = [];
            indexes(i) = [];
            i = i-1;
        end
    end
    i = i+1;
end



window_mean = 3;  
window = 35;      


dist = zeros(size(data,1),size(indexes,2));
for ch=1:size(data,1)
    for i = 1:size(indexes,2)
        dist(ch,i) = mean(data(ch,indexes(i)-window_mean:indexes(i)+window_mean)) - mean(data(ch,indexes(i)-window_mean-window:indexes(i)-window+window_mean));
    end
end
mean_ch_dist = mean(dist');


left_mean_dist = mean_ch_dist(left_ib);
right_mean_dist = mean_ch_dist(right_ib);
central_mean_dist = mean_ch_dist(central_ib);
assymetry = left_mean_dist - right_mean_dist;
figure('Name','Left brain side')
plot(left_mean_dist)
figure('Name','Right brain side')
plot(right_mean_dist)
figure('Name','Central brain side')
plot(central_mean_dist)
figure('Name','Left/Right assymetry')
plot(assymetry);


%save concatenated data
%save('S_concat_restingPre_EO.mat','concat_data');
