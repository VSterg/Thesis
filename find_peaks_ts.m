function [abs_dist locs] = find_peaks_ts(data,peak_dev)

%Find all peaks corresponding to an eyeblink in a frontal electrode channel
%timeseries. This is done by calculating peak to peak distance across all
%local maxima, then discarding the maxima with a relatively small (depending 
%on the peak_dev input) peak to peak difference. The window variable is the 
%number of samples before each local maximum, where we will search for a local
%minimum to compute the peak to peak distance. Since local maxima/minima 
%could be attributed to noise, I calculate the mean around the local 
%maximum/minimum and compute the peak to peak distance using these values.
%The window_mean variable is the number of samples around the local 
%maximum/minimum that will be used to compute the mean maximum/minimum.

%Inputs :
%data : the timeseries of a frontal electrode
%peak_dev : minimum height of spikes to be saved, relative to the mean of
%the 5 maximum spikes across all spikes.

%%%%%%%%%%%%%%%%%%%%%%%%% 
window = 60;       %%%%%% function 
window_mean = 3;   %%%%%% parameters
%%%%%%%%%%%%%%%%%%%%%%%%%


data = data - mean(data);
%data = abs(data);
data = detrend(data);
[peaks locs] = findpeaks(data);
locs = locs(find(locs>window+window_mean));
abs_dist = zeros(1,size(locs,2));
for i=1:size(locs,2)-window_mean
   [local_min idx]= min(data(locs(i)-window:locs(i)));
   local_mean_max = mean(data(locs(i)-window_mean:locs(i)+window_mean));
   local_mean_min = mean(data(locs(i)-window+idx-1-window_mean:locs(i)-window+idx-1+window_mean));
   abs_dist(i) = local_mean_max - local_mean_min;
end

av_max_peak = mean(maxk(abs_dist,5));
thres_dist = abs_dist > av_max_peak*peak_dev;
locs = thres_dist.*locs;
abs_dist = abs_dist(find(locs>0));
locs = locs(find(locs>0));


i = 2;
while i <= size(locs,2)
    if locs(i)-locs(i-1) < 150
        if abs(data(locs(i))) > abs(data(locs(i-1)))
            abs_dist(i-1) = [];
            locs(i-1) = [];
            i = i-1;
        else
            abs_dist(i) = [];
            locs(i) = [];
            i = i-1;
        end
    end
    i = i+1;
end

% while i <= size(locs,2)
%     if locs(i)-locs(i-1) < 150
%         if abs_dist(i) > abs_dist(i-1)
%             abs_dist(i-1) = [];
%             locs(i-1) = [];
%             i = i-1;
%         else
%             abs_dist(i) = [];
%             locs(i) = [];
%             i = i-1;
%         end
%     end
%     i = i+1;
% end