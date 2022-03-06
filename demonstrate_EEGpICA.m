%The purpose of this script is to demonstrate that ICA's implicit
%localisation of components is inaccurate. For this purpose we load eyes open EEG data,
%perform ICA, remove the eyeblink component and then compare a frontal, an
%occipital and a central electrode to see how ICA cleaning the data has
%affected each of them.
%If the localisation is correct, we would expect minimal variance at the
%occipital electrode timeseries. However that is not always the case.
clear all
%load (eyes open) EEG dat
EEG_data = load('S03_restingPre_EO.mat').dataRest;

%channel names are copied from database documentation


channels = ["Fp1" "AF7" 'AF3' 'F1' 'F3' 'F5' 'F7' 'FT7' 'FC5' 'FC3' 'FC1' 'C1' 'C3' 'C5' 'T7' 'TP7' 'CP5' 'CP3' 'CP1' 'P1' 'P3' 'P5' 'P7' 'P9' 'PO7' 'PO3' 'O1' 'Iz' 'Oz' 'POz' 'Pz' 'CPz' 'Fpz' 'Fp2' 'AF8' 'AF4' 'Afz' 'Fz' 'F2' 'F4' 'F6' 'F8' 'FT8' 'FC6' 'FC4' 'FC2' 'FCz' 'Cz' 'C2' 'C4' 'C6' 'T8' 'TP8' 'CP6' 'CP4' 'CP2' 'P2' 'P4' 'P6' 'P8' 'P10' 'PO8' 'PO4' 'O2'];
EEG_data(65:end,:) = [];
channels(65:end) = [];
frontal = [find(channels == 'Fp1') find(channels == 'Fp2')];
central = [find(channels == 'C1') find(channels == 'C2')];
occipital = [find(channels == 'O1') find(channels == 'O2')];
figure('Name','Central electrode')
plot(EEG_data(central(1),:)')
figure('Name','Frontal electrode')
plot(EEG_data(frontal(1),:)')
figure('Name','Occipital electrode')
plot(EEG_data(occipital(2),:)')

figure
plotmts2(EEG_data(1:15,:)')

%% Find independent components
% Here for ICA I use fastICA algorithm (package included, for details see 
% the corresponding functions). Note: in the original paper runica from EEGLAB
% was used. You can also test other ICA algorithms at this step.
%
% Note, the use of long (in time) data sets REDUCES the quality of artifact
% suppression (for details see the abovementioned paper).
% Split long files into segments and clean them separately.
K = size(EEG_data,1);
Fs = 256;
[icaEEG, A, W] = fastica(EEG_data,'stabilization','on','verbose','on'); 
ICA_data = A*icaEEG;

%plot ICA components to manually find eyeblink (noise) component
figure
plotmts2(icaEEG(1:15,:)')
figure
plotmts2(icaEEG(16:30,:)')
figure
plotmts2(icaEEG(31:45,:)')
figure
plotmts2(icaEEG(45:end,:)')

%%
noiseindex = 5;
ICA_noise = icaEEG(noiseindex,:);
icaEEG2 = zeros(size(icaEEG));
icaEEG2(noiseindex,:) = ICA_noise;
elec_noise = A*icaEEG2;

figure
plotmts2(elec_noise(1:15,:)')
elec_sd_calm = std(elec_noise(:,13960:16700)');

%% Remove noise - normal ICA
%The purpose of this section is to demonstrate that cleaning eyeblink
%component via ICA will cause significant variance in occipital electrode
%timeseries (as well as frontal electrode timeseries)

noiseindex = 5; %change this manually to the noise component index
icaEEG3 = icaEEG;
icaEEG3(noiseindex,:) = 0; %remove eyeblink component
clean_data = A*icaEEG3; %calculate clean channel data

%compare cleaned data frontal and occipital electrode timeseries
%are there spikes in any of them?
figure('Name','Frontal electrode (ICA cleaned)')
plot(clean_data(frontal(1),:))
figure('Name','Occipital electrode (ICA cleaned)')
plot(clean_data(occipital(1),:))
figure('Name','Central electrode (ICA cleaned)')
plot(clean_data(central(1),:))
figure('Name','Central electrode (ICA cleaned)')
plot(EEG_data(central(1),:))
figure
plot(clean_data(25,:))
figure
plotmts2(EEG_data(25,:)')

%calculate variance caused by removing eyeblink component via ICA
diff = EEG_data - clean_data;
%Compare variance caused on frontal, occipital and central electrode timeseries
%Are there spikes in all of them?
figure('Name','Frontal variance (ICA cleaned)')
plot(diff(frontal(1),:))
figure('Name','Occipital variance (ICA cleaned)')
plot(diff(occipital(1),:))
figure('Name','Central variance (ICA cleaned)')
plot(diff(central(1),:))

%% Remove all but one component
load emptyEEG %contains channel location info

%The purpose of this section is to demonstrate that the ICA component
%contributes significantly to occipital electrodes as well 

keep = 7; %change this manually to the kept component index
icaEEG4 = zeros(size(icaEEG)); %zero out all ICA components
icaEEG4(keep,:) = icaEEG(keep,:); %restore eyeblink component
% dc_part = mean(icaEEG4(keep,1:50));
% icaEEG4(keep,:) = icaEEG4(keep,:) - dc_part; %center to 0
single_comp_data = A*icaEEG4; %calculate single-component channel data

%Compare frontal, central and occipital electrode timeseries
%Are spikes/eyeblink activations present in both?
%How do they compare?
figure('Name','Frontal electrode (eyeblink component projection)')
plot(single_comp_data(frontal(1),:))
figure('Name','Occipital electrode (eyeblink component projection)')
plot(single_comp_data(18,:))
figure('Name','Central electrode (eyeblink component projection)')
plot(single_comp_data(56,:))

channel_dev = std(single_comp_data');
%channel_dev = normalize(channel_dev);
channel_dev = channel_dev/max(channel_dev);

figure()
topoplotIndie(channel_dev, EEG.chanlocs,'numcontour',0,'electrodes','numbers','shading','interp');
set(gca,'clim',[min(channel_dev) max(channel_dev)])
%set(gca,'clim',[-1 1])
colorbar('southoutside')
title('ICA eyeblink component projection')
%%
kee = EEG;
chanlocs = EEG.chanlocs
eeglab
%%
EEG = kee;
%%
EEG.data = EEG_data;
EEG.trials = 1;
%%
ICA_comps = EEG.icawinv' * EEG.data;
figure
plotmts2(ICA_comps(1:15,:)')
figure
plotmts2(ICA_comps(16:30,:)')
figure
plotmts2(ICA_comps(31:45,:)')
figure
plotmts2(ICA_comps(46:end,:)')
%%
idx = 35;
A = EEG.icawinv;
icaEEG = ICA_comps;
%%
data_frontal = EEG_data(frontal(1),:);
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

aa = find_activity(data_frontal,0.5);

durs = [];
i = 1;
while i < size(aa,1)
    if aa(i) == 1
        count = 0;
        while aa(i) == 1
            count = count + 1;
            i = i + 1;
        end
        if count>20
            durs = [durs count];
        end
    end
    i = i+1;
end

mean_duration = mean(durs);

times = (0:mean_duration-1)/256;
peaktime = median(times); % seconds
width    = .09;
sinefreq = 7; % for sine wave

% create Gaussian taper
gaus = 0.1*exp( -(times-peaktime).^2 / (2*width^2) );
figure
plot(gaus)

sim_eyeblink = zeros(size(aa,1),1);
mean_duration = floor(mean_duration);
half_duration = floor(mean_duration/2);

for i=1:size(indexes,2)
    idx = indexes(i);
    sim_eyeblink(idx-half_duration:idx+half_duration) = sim_eyeblink(idx-half_duration:idx+half_duration) + gaus';
end


cov(sim_eyeblink,data_frontal);
data_central = EEG_data(frontal(2),:);
ts = cov(sim_eyeblink,data_central)
%%
maxcov = 0;
maxcov_idx = 0;
for i=1:size(icaEEG,1)
    comp_eyeblink_cov = cov(icaEEG(i,:),sim_eyeblink);
    if comp_eyeblink_cov(1,2) > maxcov
        maxcov = comp_eyeblink_cov(1,2);
        maxcov_idx = i;
    end
end
icaEEG2 = icaEEG;
icaEEG2(maxcov_idx,:) = 0;
maxcov = 0;
maxcov_idx = 0;
for i=1:size(icaEEG,1)
    comp_eyeblink_cov = cov(icaEEG2(i,:),sim_eyeblink);
    if comp_eyeblink_cov(1,2) > maxcov
        maxcov = comp_eyeblink_cov(1,2);
        maxcov_idx = i;
    end
end
ta = build_eyeblink_component(data_frontal);
%%
[weights,sphere] = runica(EEG_data);
xex = weights * EEG_data;
plotmts(xex(1:15,:)',1,1,K,1/Fs,[],1)
plotmts(xex(16:30,:)',1,1,K,1/Fs,[],2)
plotmts(xex(31:45,:)',1,1,K,1/Fs,[],3)
plotmts(xex(45:end,:)',1,1,K,1/Fs,[],4)
%%
strip.a = [1 33 34 2 3 37 36 35]; %front of the head strip
strip.b = [4:11 38:47];
strip.c = [12:19 32 48:56];
strip.d = [20:26 30 31 57:63];
strip.e = [27 28 29 64]; %back of the head strip
names = ['a' 'b' 'c' 'd' 'e'];
eyeblink_noise_strips = zeros(size(EEG_data));
data_frontal = EEG_data(frontal(1),:);
sim_eyeblink =  icaEEG(noiseindex,:);
for i=1:5
    current_strip = names(i);
    current_strip = strip.(current_strip);
    EEG_strip = EEG_data(current_strip,:);
    [icaEEG_strip, A_strip, W_strip] = fastica(EEG_strip,'stabilization','on','verbose','on'); 
    [eyeblink_idx_strip best_match_strip] = find_eyeblink_comp(sim_eyeblink,icaEEG_strip);

    noiseindex_strip = best_match_strip;
    eyeblink_noise_strips(current_strip,:) = A_strip(:,noiseindex_strip)*icaEEG_strip(noiseindex_strip,:);
    figure

    plotmts2(icaEEG_strip')
end
    
EEG_strips = EEG_data - eyeblink_noise_strips;

strips_sd_calm = std(eyeblink_noise_strips(:,13960:16700)');
figure
plot(eyeblink_noise_strips(1,13960:16700))

figure
to_plot = EEG_data - elec_noise;
plotmts2(to_plot(54:64,:)')
figure
plotmts2(EEG_strips(54:64,:)')
%%

ICA_std = std(ICA_noise)
strips_std = std(strips_noise)
figure
plotmts2(clean_data(strip.d,:)')
figure
plotmts2(EEG_strips(strip.d,:)')
figure
plotmts2(EEG_data(strip.d,:)')