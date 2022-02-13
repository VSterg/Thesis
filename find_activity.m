function activity = find_activity(signal,threshold)
% The purpose of this function is to find the active periods in an eyeblink
% component timeseries.

% Inputs: signal is a 1xN vector representing an eyeblink component timeseries

%Outputs : activity is a 1xN logical vector indicating the time periods
%when the eyeblink component timeseries exhibits strong activity.

signal = abs(signal);

activity = smooth(signal,8,'rloess'); %smooth out curve to attenuate noise effect
activity = activity - median(activity); %center inactive periods at 0
activity = (activity - min(activity))/(max(activity) - min(activity)); %normalize to [0 1]

activity = activity > threshold;

%activity = gradient(cumsum(activity));
% curve_fit = fit((1:1:length(signal))',activity','cubicinterp');

% activity = curve_fit((1:1:10000)');
% activity = activity + circshift(activity,3) + circshift(activity,-3);
% activity = activity/3;
% activity = activity > threshold;
% activity = activity.*circshift(activity,window-10); % -10 because there's eyeblinks that are too close together





