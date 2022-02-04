clear all
%load emptyEEG %contains channel location info
load dsi_empty_EEG
headmodel = load('dsi_hm.mat');
projection = load('eyeblink_projection.mat'); %load the projection matrix for the eyeblink component
proj_arrangement = projection.channels;
projection_matrix = projection.mean_ch_dist;
true_proj_norm = projection_matrix'/max(abs(projection_matrix));

%The indexing of electrodes between the projection matrix and our EEGlab struct is the
%same, so there's no need to rearrange them.
leadfield(:,1,:) = headmodel.K;
%leadfield = -1*leadfield;
%%
kee = EEG;
chanlocs = kee.chanlocs;
eeglab
%% after loading .mat into eeglab, update some values
EEG = kee;


%%
Fs = 256;

sources_to_keep = 1:7; %keep 7 sources
K = size(sources_to_keep,2);
n = 5000;
frequency_bands = {'D'; 'Th'; 'A'; 'B'; 'low'};
%sim_eyeblink = load('artif_blink.mat').whole; %load simulated eyeblink component (WN + spikes)
sim_eyeblink = load('sim_eyeblink_2.mat').sim_eyeblink; %no WN
sim_eyeblink = downsample(sim_eyeblink,2); % downsample to reduce to 5000 samples
EEG_times = (0:n-1)/Fs;
chan_num = EEG.nbchan;
%dipole locations manually chosen
regions = headmodel.atlas.colorTable;
active_regions = [1 2 5 6 7 9 33 37 55 58 60 61 64];
active_regions = active_regions(1:size(sources_to_keep,2));
inactive_regions = 1:68;
inactive_regions(active_regions) = [];
dipoles_per_region = 10;

%Uncomment for VARP data
[xM,Adj_mat] = VAR1RingStructure(n,K); %generate simulated data
xM = xM(1:n,sources_to_keep);
Adj_mat = Adj_mat > 0.01;
Adj_mat = Adj_mat';
Adj_mat_ext = extend_mat(Adj_mat,active_regions,68);
%divide by 2 to account for downsampling
calm_start = 4000/2; %start of calm period (no eyeblinks)
calm_end = 6000/2; %end of calm period
act_start = 7000/2; %start of active period
act_end = 8700/2; %end of active period
alpha = 0.1;
alpha2 = 0.5;
fftfreqs = linspace(0,Fs/2,floor(n/2)+1); %frequency resolution of FFT
coh_freq = linspace(0,Fs/2,257); %frequency resolution of coherence

xM = xM(1:n,1:K);
xM_ext = extend_ts(xM,active_regions,68);
%compute CGCI for the whole period, the calm period and the active period
[CGCIM,pCGCIM] = CGCinall(xM,1,1);
[CGCIM_calm,pCGCIM_calm] = CGCinall(xM(calm_start:calm_end,:),1,1);
[CGCIM_act,pCGCIM_act] = CGCinall(xM(act_start:act_end,:),1,1);

%[CGCIM_ext,pCGCIM_ext] = CGCinall(xM_ext,1,1);
% [CGCIM_calm_ext,pCGCIM_calm_ext] = CGCinall(xM_ext(calm_start:calm_end,:),1,1);
% [CGCIM_act_ext,pCGCIM_act_ext] = CGCinall(xM_ext(act_start:act_end,:),1,1);
% pmax = 2;
% maketest = 1;
% [aaRCGCIM2,aaa]=mBTSCGCImatrix(xM_ext,pmax,maketest);
tic
figure('Name','Dipole/source timeseries')
plotmts2(xM)
xM_metrics = compute_source_metrics(xM',xM',Fs,calm_start,calm_end,act_start,act_end);
xM_ext_metrics = compute_source_metrics_ext(xM_ext',Fs,calm_start,calm_end,act_start,act_end);

[xM_metrics.CGCIM_missclassifications xM_metrics.CGCIM_confmat xM_metrics.CGCIM_confmat_metrics] = compute_network_variance(Adj_mat ,xM_metrics.CGCIM_network);
[xM_metrics.CGCIM_missclassifications_calm xM_metrics.CGCIM_confmat_calm xM_metrics.CGCIM_confmat_calm_metrics] = compute_network_variance(Adj_mat,xM_metrics.CGCIM_network_calm);
[xM_metrics.CGCIM_missclassifications_act xM_metrics.CGCIM_confmat_act xM_metrics.CGCIM_confmat_act_metrics] = compute_network_variance(Adj_mat,xM_metrics.CGCIM_network_act);

% [xM_ext_metrics.CGCIM_missclassifications xM_ext_metrics.CGCIM_confmat xM_ext_metrics.CGCIM_confmat_metrics] = compute_network_variance(Adj_mat_ext ,xM_ext_metrics.CGCIM_network);
% [xM_ext_metrics.CGCIM_missclassifications_calm xM_ext_metrics.CGCIM_confmat_calm xM_ext_metrics.CGCIM_confmat_calm_metrics] = compute_network_variance(Adj_mat_ext,xM_ext_metrics.CGCIM_network_calm);
% [xM_ext_metrics.CGCIM_missclassifications_act xM_ext_metrics.CGCIM_confmat_act xM_ext_metrics.CGCIM_confmat_act_metrics] = compute_network_variance(Adj_mat_ext,xM_ext_metrics.CGCIM_network_act);

[xM_ext_metrics.RCGCIM_missclassifications xM_ext_metrics.RCGCIM_confmat xM_ext_metrics.RCGCIM_confmat_metrics] = compute_network_variance(Adj_mat_ext ,xM_ext_metrics.RCGCIM_network);
[xM_ext_metrics.RCGCIM_missclassifications_calm xM_ext_metrics.RCGCIM_confmat_calm xM_ext_metrics.RCGCIM_confmat_calm_metrics] = compute_network_variance(Adj_mat_ext,xM_ext_metrics.RCGCIM_network_calm);
[xM_ext_metrics.RCGCIM_missclassifications_act xM_ext_metrics.RCGCIM_confmat_act xM_ext_metrics.RCGCIM_confmat_act_metrics] = compute_network_variance(Adj_mat_ext,xM_ext_metrics.RCGCIM_network_act);

for i=1:size(frequency_bands,1)
band = frequency_bands(i);
band = string(band);
[xM_metrics.DTF_missclassifications.(band) xM_metrics.DTF_confmat.(band) xM_metrics.DTF_confmat_metrics.(band)] = compute_network_variance(Adj_mat,xM_metrics.DTF_network.(band));
[xM_metrics.DTF_missclassifications_calm.(band) xM_metrics.DTF_confmat_calm.(band) xM_metrics.DTF_confmat_calm_metrics.(band)] = compute_network_variance(Adj_mat,xM_metrics.DTF_network_calm.(band));
[xM_metrics.DTF_missclassifications_act.(band) xM_metrics.DTF_confmat_act.(band) xM_metrics.DTF_confmat_act_metrics.(band)] = compute_network_variance(Adj_mat,xM_metrics.DTF_network_act.(band));

[xM_ext_metrics.DTF_missclassifications_act.(band) xM_ext_metrics.DTF_confmat_act.(band) xM_ext_metrics.DTF_confmat_act_metrics.(band)] = compute_network_variance(Adj_mat_ext,xM_ext_metrics.DTF_network_act.(band));

[xM_metrics.PDC_missclassifications.(band) xM_metrics.PDC_confmat.(band) xM_metrics.PDC_confmat_metrics.(band)] = compute_network_variance(Adj_mat,xM_metrics.PDC_network.(band));
[xM_metrics.PDC_missclassifications_calm.(band) xM_metrics.PDC_confmat_calm.(band) xM_metrics.PDC_confmat_calm_metrics.(band)] = compute_network_variance(Adj_mat,xM_metrics.PDC_network_calm.(band));
[xM_metrics.PDC_missclassifications_act.(band) xM_metrics.PDC_confmat_act.(band) xM_metrics.PDC_confmat_act_metrics.(band)] = compute_network_variance(Adj_mat,xM_metrics.PDC_network_act.(band));

[xM_ext_metrics.PDC_missclassifications_act.(band) xM_ext_metrics.PDC_confmat_act.(band) xM_ext_metrics.PDC_confmat_act_metrics.(band)] = compute_network_variance(Adj_mat_ext,xM_ext_metrics.PDC_network_act.(band));
end
toc
%% Dipole projection to regions


for i=1:K
    tmpidx = find(regions == active_regions(i));
    diploc(i,:) = tmpidx(floor(size(tmpidx,1)/2):floor(size(tmpidx,1)/2)+dipoles_per_region-1);
end

for i=1:68-K
    tmpidx = find(regions == inactive_regions(i));
    diploc_pinknoise(i,:) = tmpidx(randsample(size(tmpidx,1),5));
    %diploc(i,:) = tmpidx(floor(size(tmpidx,1)/2):floor(size(tmpidx,1)/2)+dipoles_per_region-1);
end

% initialize all dipole data
dipole_data = zeros(size(leadfield,3),n);

%the simulated neural components are added to manually selected dipoles
for i = 2:K+1 
    dipole_data(diploc(i-1,:),:) = repmat(xM(:,i-1)',10,1);
end

% add pink noise to all other ROIs (maybe comment this out?)
for i=2:68-K+1
    dipole_data(diploc_pinknoise(i-1,:),:) = repmat(pinknoise(n)',5,1);
end

EEG_clean(:,:) = squeeze(leadfield(:,1,:))*dipole_data; %project dipole data to electrodes (no eyeblink component)
sim_eyeblink = sim_eyeblink(1:n); %n must not be higher than 10000

%add eyeblink component projection (noise) to EEG data
amp = 400;
EEG_noisy = EEG_clean + projection_matrix'*sim_eyeblink/amp;
figure('Name','Noisy data, first 10 EEG channels')
plotmts2(EEG_noisy(1:10,:)')    %plotmts(EEG_noisy(1:10,:)',1,1,chan_num,1/Fs)
figure('Name','Clean data, first 10 EEG channels')
plotmts2(EEG_clean(1:10,:)')    %plotmts(EEG_clean(50:60,:)',1,1,chan_num,1/Fs,[],2)

my_dipoles_gain = leadfield(:,1,diploc);
my_dipoles_data = dipole_data([diploc],:);%store the timeseries of active dipoles
EEG = kee;
EEG.data = EEG_clean;

%% clean data source reconstruction
tic
saveFull = false;
account4artifacts = false; 
src2roiReductionType = 'sum';
solverType = 'bsbl'
updateFreq = 5;
asec = pop_rsbl(EEG, saveFull, account4artifacts, src2roiReductionType, solverType, updateFreq);

%all 68 source
ROI_acts = asec.etc.src.act;
%ROI_acts = ROI_acts(active_regions,:);

ROI_acts = ROI_acts/dipoles_per_region;


channel_metrics_clean = Compute_channel_metrics(EEG_clean,EEG_clean,calm_start,calm_end,act_start,act_end);

[channel_metrics_clean.RCGCIM_missclassifications channel_metrics_clean.RCGCIM_confmat channel_metrics_clean.RCGCIM_confmat_metrics] = compute_network_variance(channel_metrics_clean.RCGCIM_network ,channel_metrics_clean.RCGCIM_network);
[channel_metrics_clean.RCGCIM_missclassifications_calm channel_metrics_clean.RCGCIM_confmat_calm channel_metrics_clean.RCGCIM_confmat_calm_metrics] = compute_network_variance(channel_metrics_clean.RCGCIM_network_calm,channel_metrics_clean.RCGCIM_network_calm);
[channel_metrics_clean.RCGCIM_missclassifications_act channel_metrics_clean.RCGCIM_confmat_act channel_metrics_clean.RCGCIM_confmat_act_metrics] = compute_network_variance(channel_metrics_clean.RCGCIM_network_act,channel_metrics_clean.RCGCIM_network_act);

source_metrics_clean_ext = compute_source_metrics_ext(ROI_acts,Fs,calm_start,calm_end,act_start,act_end);

% [source_metrics_clean_ext.CGCIM_missclassifications source_metrics_clean_ext.CGCIM_confmat source_metrics_clean_ext.CGCIM_confmat_metrics] = compute_network_variance(Adj_mat_ext ,source_metrics_clean_ext.CGCIM_network);
% [source_metrics_clean_ext.CGCIM_missclassifications_calm source_metrics_clean_ext.CGCIM_confmat_calm source_metrics_clean_ext.CGCIM_confmat_calm_metrics] = compute_network_variance(Adj_mat_ext,source_metrics_clean_ext.CGCIM_network_calm);
% [source_metrics_clean_ext.CGCIM_missclassifications_act source_metrics_clean_ext.CGCIM_confmat_act source_metrics_clean_ext.CGCIM_confmat_act_metrics] = compute_network_variance(Adj_mat_ext,source_metrics_clean_ext.CGCIM_network_act);
[source_metrics_clean_ext.RCGCIM_missclassifications source_metrics_clean_ext.RCGCIM_confmat source_metrics_clean_ext.RCGCIM_confmat_metrics] = compute_network_variance(Adj_mat_ext ,source_metrics_clean_ext.RCGCIM_network);
[source_metrics_clean_ext.RCGCIM_missclassifications_calm source_metrics_clean_ext.RCGCIM_confmat_calm source_metrics_clean_ext.RCGCIM_confmat_calm_metrics] = compute_network_variance(Adj_mat_ext,source_metrics_clean_ext.RCGCIM_network_calm);
[source_metrics_clean_ext.RCGCIM_missclassifications_act source_metrics_clean_ext.RCGCIM_confmat_act source_metrics_clean_ext.RCGCIM_confmat_act_metrics] = compute_network_variance(Adj_mat_ext,source_metrics_clean_ext.RCGCIM_network_act);




for j=1:size(frequency_bands,1)
band = frequency_bands(j);
band = string(band);

[source_metrics_clean_ext.DTF_missclassifications_act.(band) source_metrics_clean_ext.DTF_confmat_act.(band) source_metrics_clean_ext.DTF_confmat_act_metrics.(band)] = compute_network_variance(Adj_mat_ext,source_metrics_clean_ext.DTF_network_act.(band));

[source_metrics_clean_ext.PDC_missclassifications_act.(band) source_metrics_clean_ext.PDC_confmat_act.(band) source_metrics_clean_ext.PDC_confmat_act_metrics.(band)] = compute_network_variance(Adj_mat_ext,source_metrics_clean_ext.PDC_network_act.(band));


[channel_metrics_clean.coh_missclassifications.(band) channel_metrics_clean.coh_confmat.(band) channel_metrics_clean.coh_confmat_metrics.(band) channel_metrics_clean.coh_error_network.(band)] = compute_network_variance(channel_metrics_clean.coh_network.(band),channel_metrics_clean.coh_network.(band));
[channel_metrics_clean.coh_missclassifications_calm.(band) channel_metrics_clean.coh_confmat_calm.(band) channel_metrics_clean.coh_confmat_calm_metrics.(band) channel_metrics_clean.coh_error_network_calm.(band)] = compute_network_variance(channel_metrics_clean.coh_network_calm.(band),channel_metrics_clean.coh_network_calm.(band));
[channel_metrics_clean.coh_missclassifications_act.(band) channel_metrics_clean.coh_confmat_act.(band) channel_metrics_clean.coh_confmat_act_metric.(band) channel_metrics_clean.coh_error_network_act.(band)] = compute_network_variance(channel_metrics_clean.coh_network_act.(band),channel_metrics_clean.coh_network_act.(band));
end



%8 sources only
ROI_acts = asec.etc.src.act;
ROI_acts = ROI_acts(active_regions,:);


ROI_acts = ROI_acts/dipoles_per_region;

source_metrics_clean = compute_source_metrics(ROI_acts,xM',Fs,calm_start,calm_end,act_start,act_end);

[source_metrics_clean.CGCIM_missclassifications source_metrics_clean.CGCIM_confmat source_metrics_clean.CGCIM_confmat_metrics] = compute_network_variance(Adj_mat ,source_metrics_clean.CGCIM_network);
[source_metrics_clean.CGCIM_missclassifications_calm source_metrics_clean.CGCIM_confmat_calm source_metrics_clean.CGCIM_confmat_calm_metrics] = compute_network_variance(Adj_mat,source_metrics_clean.CGCIM_network_calm);
[source_metrics_clean.CGCIM_missclassifications_act source_metrics_clean.CGCIM_confmat_act source_metrics_clean.CGCIM_confmat_act_metrics] = compute_network_variance(Adj_mat,source_metrics_clean.CGCIM_network_act);

for j=1:size(frequency_bands,1)
band = frequency_bands(j);
band = string(band);
[source_metrics_clean.DTF_missclassifications.(band) source_metrics_clean.DTF_confmat.(band) source_metrics_clean.DTF_confmat_metrics.(band)] = compute_network_variance(Adj_mat,source_metrics_clean.DTF_network.(band));
[source_metrics_clean.DTF_missclassifications_calm.(band) source_metrics_clean.DTF_confmat_calm.(band) source_metrics_clean.DTF_confmat_calm_metrics.(band)] = compute_network_variance(Adj_mat,source_metrics_clean.DTF_network_calm.(band));
[source_metrics_clean.DTF_missclassifications_act.(band) source_metrics_clean.DTF_confmat_act.(band) source_metrics_clean.DTF_confmat_act_metrics.(band)] = compute_network_variance(Adj_mat,source_metrics_clean.DTF_network_act.(band));

[source_metrics_clean.PDC_missclassifications.(band) source_metrics_clean.PDC_confmat.(band) source_metrics_clean.PDC_confmat_metrics.(band)] = compute_network_variance(Adj_mat,source_metrics_clean.PDC_network.(band));
[source_metrics_clean.PDC_missclassifications_calm.(band) source_metrics_clean.PDC_confmat_calm.(band) source_metrics_clean.PDC_confmat_calm_metrics.(band)] = compute_network_variance(Adj_mat,source_metrics_clean.PDC_network_calm.(band));
[source_metrics_clean.PDC_missclassifications_act.(band) source_metrics_clean.PDC_confmat_act.(band) source_metrics_clean.PDC_confmat_act_metrics.(band)] = compute_network_variance(Adj_mat,source_metrics_clean.PDC_network_act.(band));

%already did that
% [channel_metrics_clean.coh_missclassifications.(band) channel_metrics_clean.coh_confmat.(band) channel_metrics_clean.coh_confmat_metrics.(band) channel_metrics_clean.coh_error_network.(band)] = compute_network_variance(channel_metrics_clean.coh_network.(band),channel_metrics_clean.coh_network.(band));
% [channel_metrics_clean.coh_missclassifications_calm.(band) channel_metrics_clean.coh_confmat_calm.(band) channel_metrics_clean.coh_confmat_calm_metrics.(band) channel_metrics_clean.coh_error_network_calm.(band)] = compute_network_variance(channel_metrics_clean.coh_network_calm.(band),channel_metrics_clean.coh_network_calm.(band));
% [channel_metrics_clean.coh_missclassifications_act.(band) channel_metrics_clean.coh_confmat_act.(band) channel_metrics_clean.coh_confmat_act_metric.(band) channel_metrics_clean.coh_error_network_act.(band)] = compute_network_variance(channel_metrics_clean.coh_network_act.(band),channel_metrics_clean.coh_network_act.(band));
end



figure
plotmts2(xM)
figure
plotmts2(-ROI_acts')

toc

%% noisy data source reconstruction

EEG.data = EEG_noisy;

asec = pop_rsbl(EEG, saveFull, account4artifacts, src2roiReductionType, solverType, updateFreq);

%all 68 source
ROI_acts = asec.etc.src.act;
%ROI_acts = ROI_acts(active_regions,:);

ROI_acts = ROI_acts/dipoles_per_region;
a = 2
tic
channel_metrics_noisy = Compute_channel_metrics(EEG_noisy,EEG_clean,calm_start,calm_end,act_start,act_end);

[channel_metrics_noisy.RCGCIM_missclassifications channel_metrics_noisy.RCGCIM_confmat channel_metrics_noisy.RCGCIM_confmat_metrics channel_metrics_noisy.RCGCIM_error_network] = compute_network_variance(channel_metrics_clean.RCGCIM_network ,channel_metrics_noisy.RCGCIM_network);
[channel_metrics_noisy.RCGCIM_missclassifications_calm channel_metrics_noisy.RCGCIM_confmat_calm channel_metrics_noisy.RCGCIM_confmat_calm_metrics channel_metrics_noisy.RCGCIM_error_network_calm] = compute_network_variance(channel_metrics_clean.RCGCIM_network_calm,channel_metrics_noisy.RCGCIM_network_calm);
[channel_metrics_noisy.RCGCIM_missclassifications_act channel_metrics_noisy.RCGCIM_confmat_act channel_metrics_noisy.RCGCIM_confmat_act_metrics channel_metrics_noisy.RCGCIM_error_network_act] = compute_network_variance(channel_metrics_clean.RCGCIM_network_act,channel_metrics_noisy.RCGCIM_network_act);
toc
a = 3
source_metrics_noisy_ext = compute_source_metrics_ext(ROI_acts,Fs,calm_start,calm_end,act_start,act_end);
a = 4
tic
% [source_metrics_noisy_ext.CGCIM_missclassifications source_metrics_noisy_ext.CGCIM_confmat source_metrics_noisy_ext.CGCIM_confmat_metrics] = compute_network_variance(Adj_mat_ext ,source_metrics_noisy_ext.CGCIM_network);
% [source_metrics_noisy_ext.CGCIM_missclassifications_calm source_metrics_noisy_ext.CGCIM_confmat_calm source_metrics_noisy_ext.CGCIM_confmat_calm_metrics] = compute_network_variance(Adj_mat_ext,source_metrics_noisy_ext.CGCIM_network_calm);
% [source_metrics_noisy_ext.CGCIM_missclassifications_act source_metrics_noisy_ext.CGCIM_confmat_act source_metrics_noisy_ext.CGCIM_confmat_act_metrics] = compute_network_variance(Adj_mat_ext,source_metrics_noisy_ext.CGCIM_network_act);
[source_metrics_noisy_ext.RCGCIM_missclassifications source_metrics_noisy_ext.RCGCIM_confmat source_metrics_noisy_ext.RCGCIM_confmat_metrics source_metrics_noisy_ext.RCGCIM_error_network] = compute_network_variance(Adj_mat_ext ,source_metrics_noisy_ext.RCGCIM_network);
[source_metrics_noisy_ext.RCGCIM_missclassifications_calm source_metrics_noisy_ext.RCGCIM_confmat_calm source_metrics_noisy_ext.RCGCIM_confmat_calm_metrics source_metrics_noisy_ext.RCGCIM_error_network_calm] = compute_network_variance(Adj_mat_ext,source_metrics_noisy_ext.RCGCIM_network_calm);
[source_metrics_noisy_ext.RCGCIM_missclassifications_act source_metrics_noisy_ext.RCGCIM_confmat_act source_metrics_noisy_ext.RCGCIM_confmat_act_metrics source_metrics_noisy_ext.RCGCIM_error_network_act] = compute_network_variance(Adj_mat_ext,source_metrics_noisy_ext.RCGCIM_network_act);

%compare vs clean data metrics
[source_metrics_noisy_ext.RCGCIM_missclassifications_vs_clean source_metrics_noisy_ext.RCGCIM_confmat_vs_clean source_metrics_noisy_ext.RCGCIM_confmat_metrics_vs_clean] = compute_network_variance(source_metrics_clean_ext.RCGCIM_network ,source_metrics_noisy_ext.RCGCIM_network);
[source_metrics_noisy_ext.RCGCIM_missclassifications_calm_vs_clean source_metrics_noisy_ext.RCGCIM_confmat_calm_vs_clean source_metrics_noisy_ext.RCGCIM_confmat_calm_metrics_vs_clean] = compute_network_variance(source_metrics_clean_ext.RCGCIM_network_calm,source_metrics_noisy_ext.RCGCIM_network_calm);
[source_metrics_noisy_ext.RCGCIM_missclassifications_act_vs_clean source_metrics_noisy_ext.RCGCIM_confmat_act_vs_clean source_metrics_noisy_ext.RCGCIM_confmat_act_metrics_vs_clean] = compute_network_variance(source_metrics_clean_ext.RCGCIM_network_act,source_metrics_noisy_ext.RCGCIM_network_act);


for j=1:size(frequency_bands,1)
band = frequency_bands(j);
band = string(band);

[source_metrics_noisy_ext.DTF_missclassifications_act.(band) source_metrics_noisy_ext.DTF_confmat_act.(band) source_metrics_noisy_ext.DTF_confmat_act_metrics.(band) source_metrics_noisy_ext.DTF_error_network_act.(band)] = compute_network_variance(Adj_mat_ext,source_metrics_noisy_ext.DTF_network_act.(band));

[source_metrics_noisy_ext.PDC_missclassifications_act.(band) source_metrics_noisy_ext.PDC_confmat_act.(band) source_metrics_noisy_ext.PDC_confmat_act_metrics.(band) source_metrics_noisy_ext.PDC_error_network_act.(band)] = compute_network_variance(Adj_mat_ext,source_metrics_noisy_ext.PDC_network_act.(band));
[source_metrics_noisy_ext.PDC_missclassifications_act_vs_clean.(band) source_metrics_noisy_ext.PDC_confmat_act_vs_clean.(band) source_metrics_noisy_ext.PDC_confmat_act_metrics_vs_clean.(band) source_metrics_noisy_ext.PDC_error_network_act_vs_clean.(band)] = compute_network_variance(source_metrics_clean_ext.PDC_network_act.(band),source_metrics_noisy_ext.PDC_network_act.(band));


[channel_metrics_noisy.coh_missclassifications.(band) channel_metrics_noisy.coh_confmat.(band) channel_metrics_noisy.coh_confmat_metrics.(band) channel_metrics_noisy.coh_error_network.(band)] = compute_network_variance(channel_metrics_clean.coh_network.(band),channel_metrics_noisy.coh_network.(band));
[channel_metrics_noisy.coh_missclassifications_calm.(band) channel_metrics_noisy.coh_confmat_calm.(band) channel_metrics_noisy.coh_confmat_calm_metrics.(band) channel_metrics_noisy.coh_error_network_calm.(band)] = compute_network_variance(channel_metrics_clean.coh_network_calm.(band),channel_metrics_noisy.coh_network_calm.(band));
[channel_metrics_noisy.coh_missclassifications_act.(band) channel_metrics_noisy.coh_confmat_act.(band) channel_metrics_noisy.coh_confmat_act_metric.(band) channel_metrics_noisy.coh_error_network_act.(band)] = compute_network_variance(channel_metrics_clean.coh_network_act.(band),channel_metrics_noisy.coh_network_act.(band));
end




%8 source

ROI_acts = asec.etc.src.act;

ROI_acts = ROI_acts(active_regions,:);
ROI_acts = ROI_acts/dipoles_per_region;

source_metrics_noisy = compute_source_metrics(ROI_acts,xM',Fs,calm_start,calm_end,act_start,act_end);
[source_metrics_noisy.CGCIM_missclassifications source_metrics_noisy.CGCIM_confmat source_metrics_noisy.CGCIM_confmat_metrics source_metrics_noisy.CGCIM_error_network] = compute_network_variance(Adj_mat,source_metrics_noisy.CGCIM_network);
[source_metrics_noisy.CGCIM_missclassifications_calm source_metrics_noisy.CGCIM_confmat_calm source_metrics_noisy.CGCIM_confmat_calm_metrics source_metrics_noisy.CGCIM_error_network_calm] = compute_network_variance(Adj_mat,source_metrics_noisy.CGCIM_network_calm);
[source_metrics_noisy.CGCIM_missclassifications_act source_metrics_noisy.CGCIM_confmat_act source_metrics_noisy.CGCIM_confmat_act_metrics source_metrics_noisy.CGCIM_error_network_act] = compute_network_variance(Adj_mat,source_metrics_noisy.CGCIM_network_act);


for i=1:size(frequency_bands,1)
band = frequency_bands(i);
band = string(band);
[source_metrics_noisy.DTF_missclassifications.(band) source_metrics_noisy.DTF_confmat.(band) source_metrics_noisy.DTF_confmat_metrics.(band) source_metrics_noisy.DTF_error_network.(band)] = compute_network_variance(Adj_mat,source_metrics_noisy.DTF_network.(band));
[source_metrics_noisy.DTF_missclassifications_calm.(band) source_metrics_noisy.DTF_confmat_calm.(band) source_metrics_noisy.DTF_confmat_calm_metrics.(band) source_metrics_noisy.DTF_error_network_calm.(band)] = compute_network_variance(Adj_mat,source_metrics_noisy.DTF_network_calm.(band));
[source_metrics_noisy.DTF_missclassifications_act.(band) source_metrics_noisy.DTF_confmat_act.(band) source_metrics_noisy.DTF_confmat_act_metrics.(band) source_metrics_noisy.DTF_error_network_act.(band)] = compute_network_variance(Adj_mat,source_metrics_noisy.DTF_network_act.(band));

[source_metrics_noisy.PDC_missclassifications.(band) source_metrics_noisy.PDC_confmat.(band) source_metrics_noisy.PDC_confmat_metrics.(band) source_metrics_noisy.PDC_error_network.(band)] = compute_network_variance(Adj_mat,source_metrics_noisy.PDC_network.(band));
[source_metrics_noisy.PDC_missclassifications_calm.(band) source_metrics_noisy.PDC_confmat_calm.(band) source_metrics_noisy.PDC_confmat_calm_metrics.(band) source_metrics_noisy.PDC_error_network_calm.(band)] = compute_network_variance(Adj_mat,source_metrics_noisy.PDC_network_calm.(band));
[source_metrics_noisy.PDC_missclassifications_act.(band) source_metrics_noisy.PDC_confmat_act.(band) source_metrics_noisy.PDC_confmat_act_metrics.(band) source_metrics_noisy.PDC_error_network_act.(band)] = compute_network_variance(Adj_mat,source_metrics_noisy.PDC_network_act.(band));

%already did that
% [channel_metrics_noisy.coh_missclassifications.(band) channel_metrics_noisy.coh_confmat.(band) channel_metrics_noisy.coh_confmat_metrics.(band) channel_metrics_noisy.coh_error_network.(band)] = compute_network_variance(channel_metrics_clean.coh_network.(band),channel_metrics_noisy.coh_network.(band));
% [channel_metrics_noisy.coh_missclassifications_calm.(band) channel_metrics_noisy.coh_confmat_calm.(band) channel_metrics_noisy.coh_confmat_calm_metrics.(band) channel_metrics_noisy.coh_error_network_calm.(band)] = compute_network_variance(channel_metrics_clean.coh_network_calm.(band),channel_metrics_noisy.coh_network_calm.(band));
% [channel_metrics_noisy.coh_missclassifications_act.(band) channel_metrics_noisy.coh_confmat_act.(band) channel_metrics_noisy.coh_confmat_act_metrics.(band) channel_metrics_noisy.coh_error_network_act.(band)] = compute_network_variance(channel_metrics_clean.coh_network_act.(band),channel_metrics_noisy.coh_network_act.(band));
end

figure
plotmts2(ROI_acts')
toc

%% Find independent components (include ICA folder)
% Here for ICA I use fastICA algorithm (package included, for details see 
% the corresponding functions). Note: in the original paper runica from EEGLAB
% was used. You can also test other ICA algorithms at this step.
%
% Note, the use of long (in time) data sets REDUCES the quality of artifact
% suppression (for details see the abovementioned paper).
% Split long files into segments and clean them separately.

[icaEEG, A, W] = fastica(EEG_noisy,'stabilization','on','verbose','on'); 
[eyeblink_idx best_match] = find_eyeblink_comp(sim_eyeblink,icaEEG);
noiseindex = best_match %change this manually to the noise component index
icaEEG3 = icaEEG;
icaEEG3(noiseindex,:) = 0; %remove eyeblink component
EEG_ICA = A*icaEEG3; %calculate clean channel data
figure
plotmts2(icaEEG(1:11,:)')
% figure
% plotmts2(icaEEG(23:28,:)')
% figure
% plotmts2(icaEEG(15:22,:)')
%% ICA data source reconstruction
tic
EEG.data = EEG_ICA;

asec = pop_rsbl(EEG, saveFull, account4artifacts, src2roiReductionType, solverType, updateFreq);

%all 68 source
ROI_acts = asec.etc.src.act;
%ROI_acts = ROI_acts(active_regions,:);

ROI_acts = ROI_acts/dipoles_per_region;

channel_metrics_ICA = Compute_channel_metrics(EEG_ICA,EEG_clean,calm_start,calm_end,act_start,act_end);

[channel_metrics_ICA.RCGCIM_missclassifications channel_metrics_ICA.RCGCIM_confmat channel_metrics_ICA.RCGCIM_confmat_metrics channel_metrics_ICA.RCGCIM_error_network] = compute_network_variance(channel_metrics_clean.RCGCIM_network ,channel_metrics_ICA.RCGCIM_network);
[channel_metrics_ICA.RCGCIM_missclassifications_calm channel_metrics_ICA.RCGCIM_confmat_calm channel_metrics_ICA.RCGCIM_confmat_calm_metrics channel_metrics_ICA.RCGCIM_error_network_calm] = compute_network_variance(channel_metrics_clean.RCGCIM_network_calm,channel_metrics_ICA.RCGCIM_network_calm);
[channel_metrics_ICA.RCGCIM_missclassifications_act channel_metrics_ICA.RCGCIM_confmat_act channel_metrics_ICA.RCGCIM_confmat_act_metrics channel_metrics_ICA.RCGCIM_error_network_act] = compute_network_variance(channel_metrics_clean.RCGCIM_network_act,channel_metrics_ICA.RCGCIM_network_act);


source_metrics_ICA_ext = compute_source_metrics_ext(ROI_acts,Fs,calm_start,calm_end,act_start,act_end);

% [source_metrics_ICA_ext.CGCIM_missclassifications source_metrics_ICA_ext.CGCIM_confmat source_metrics_ICA_ext.CGCIM_confmat_metrics] = compute_network_variance(Adj_mat_ext ,source_metrics_ICA_ext.CGCIM_network);
% [source_metrics_ICA_ext.CGCIM_missclassifications_calm source_metrics_ICA_ext.CGCIM_confmat_calm source_metrics_ICA_ext.CGCIM_confmat_calm_metrics] = compute_network_variance(Adj_mat_ext,source_metrics_ICA_ext.CGCIM_network_calm);
% [source_metrics_ICA_ext.CGCIM_missclassifications_act source_metrics_ICA_ext.CGCIM_confmat_act source_metrics_ICA_ext.CGCIM_confmat_act_metrics] = compute_network_variance(Adj_mat_ext,source_metrics_ICA_ext.CGCIM_network_act);
[source_metrics_ICA_ext.RCGCIM_missclassifications source_metrics_ICA_ext.RCGCIM_confmat source_metrics_ICA_ext.RCGCIM_confmat_metrics source_metrics_ICA_ext.RCGCIM_error_network] = compute_network_variance(Adj_mat_ext ,source_metrics_ICA_ext.RCGCIM_network);
[source_metrics_ICA_ext.RCGCIM_missclassifications_calm source_metrics_ICA_ext.RCGCIM_confmat_calm source_metrics_ICA_ext.RCGCIM_confmat_calm_metrics source_metrics_ICA_ext.RCGCIM_error_network_calm] = compute_network_variance(Adj_mat_ext,source_metrics_ICA_ext.RCGCIM_network_calm);
[source_metrics_ICA_ext.RCGCIM_missclassifications_act source_metrics_ICA_ext.RCGCIM_confmat_act source_metrics_ICA_ext.RCGCIM_confmat_act_metrics source_metrics_ICA_ext.RCGCIM_error_network_act] = compute_network_variance(Adj_mat_ext,source_metrics_ICA_ext.RCGCIM_network_act);



for j=1:size(frequency_bands,1)
band = frequency_bands(j);
band = string(band);

[source_metrics_ICA_ext.DTF_missclassifications_act.(band) source_metrics_ICA_ext.DTF_confmat_act.(band) source_metrics_ICA_ext.DTF_confmat_act_metrics.(band) source_metrics_ICA_ext.DTF_error_network_act.(band)] = compute_network_variance(Adj_mat_ext,source_metrics_ICA_ext.DTF_network_act.(band));

[source_metrics_ICA_ext.PDC_missclassifications_act.(band) source_metrics_ICA_ext.PDC_confmat_act.(band) source_metrics_ICA_ext.PDC_confmat_act_metrics.(band) source_metrics_ICA_ext.PDC_error_network_act.(band)] = compute_network_variance(Adj_mat_ext,source_metrics_ICA_ext.PDC_network_act.(band));


[channel_metrics_ICA.coh_missclassifications.(band) channel_metrics_ICA.coh_confmat.(band) channel_metrics_ICA.coh_confmat_metrics.(band) channel_metrics_ICA.coh_error_network.(band)] = compute_network_variance(channel_metrics_clean.coh_network.(band),channel_metrics_ICA.coh_network.(band));
[channel_metrics_ICA.coh_missclassifications_calm.(band) channel_metrics_ICA.coh_confmat_calm.(band) channel_metrics_ICA.coh_confmat_calm_metrics.(band) channel_metrics_ICA.coh_error_network_calm.(band)] = compute_network_variance(channel_metrics_clean.coh_network_calm.(band),channel_metrics_ICA.coh_network_calm.(band));
[channel_metrics_ICA.coh_missclassifications_act.(band) channel_metrics_ICA.coh_confmat_act.(band) channel_metrics_ICA.coh_confmat_act_metric.(band) channel_metrics_ICA.coh_error_network_act.(band)] = compute_network_variance(channel_metrics_clean.coh_network_act.(band),channel_metrics_ICA.coh_network_act.(band));
%[metrics_clean.RCGCIM_missclassifications_act metrics_clean.RCGCIM_confmat_act] = compute_network_variance(metrics_clean.RCGCIM_network_act,metrics_clean.RCGCIM_network_act);
end




% 7 sources
ROI_acts = asec.etc.src.act;

ROI_acts = ROI_acts(active_regions,:);
ROI_acts = ROI_acts/dipoles_per_region;

source_metrics_ICA = compute_source_metrics(ROI_acts,xM',Fs,calm_start,calm_end,act_start,act_end);
toc
tic
[source_metrics_ICA.CGCIM_missclassifications source_metrics_ICA.CGCIM_confmat source_metrics_ICA.CGCIM_confmat_metrics source_metrics_ICA.CGCIM_error_network] = compute_network_variance(Adj_mat,source_metrics_ICA.CGCIM_network);
[source_metrics_ICA.CGCIM_missclassifications_calm source_metrics_ICA.CGCIM_confmat_calm source_metrics_ICA.CGCIM_confmat_calm_metrics source_metrics_ICA.CGCIM_error_network_calm] = compute_network_variance(Adj_mat,source_metrics_ICA.CGCIM_network_calm);
[source_metrics_ICA.CGCIM_missclassifications_act source_metrics_ICA.CGCIM_confmat_act source_metrics_ICA.CGCIM_confmat_act_metrics source_metrics_ICA.CGCIM_error_network_act] = compute_network_variance(Adj_mat,source_metrics_ICA.CGCIM_network_act);


for i=1:size(frequency_bands,1)
band = frequency_bands(i);
band = string(band);
[source_metrics_ICA.DTF_missclassifications.(band) source_metrics_ICA.DTF_confmat.(band) source_metrics_ICA.DTF_confmat_metrics.(band) source_metrics_ICA.DTF_error_network.(band)] = compute_network_variance(Adj_mat,source_metrics_ICA.DTF_network.(band));
[source_metrics_ICA.DTF_missclassifications_calm.(band) source_metrics_ICA.DTF_confmat_calm.(band) source_metrics_ICA.DTF_confmat_calm_metrics.(band) source_metrics_ICA.DTF_error_network_calm.(band)] = compute_network_variance(Adj_mat,source_metrics_ICA.DTF_network_calm.(band));
[source_metrics_ICA.DTF_missclassifications_act.(band) source_metrics_ICA.DTF_confmat_act.(band) source_metrics_ICA.DTF_confmat_act_metrics.(band) source_metrics_ICA.DTF_error_network_act.(band)] = compute_network_variance(Adj_mat,source_metrics_ICA.DTF_network_act.(band));

[source_metrics_ICA.PDC_missclassifications.(band) source_metrics_ICA.PDC_confmat.(band) source_metrics_ICA.PDC_confmat_metrics.(band) source_metrics_ICA.PDC_error_network.(band)] = compute_network_variance(Adj_mat,source_metrics_ICA.PDC_network.(band));
[source_metrics_ICA.PDC_missclassifications_calm.(band) source_metrics_ICA.PDC_confmat_calm.(band) source_metrics_ICA.PDC_confmat_calm_metrics.(band) source_metrics_ICA.PDC_error_network_calm.(band)] = compute_network_variance(Adj_mat,source_metrics_ICA.PDC_network_calm.(band));
[source_metrics_ICA.PDC_missclassifications_act.(band) source_metrics_ICA.PDC_confmat_act.(band) source_metrics_ICA.PDC_confmat_act_metrics.(band) source_metrics_ICA.PDC_error_network_act.(band)] = compute_network_variance(Adj_mat,source_metrics_ICA.PDC_network_act.(band));

%already did that
% [channel_metrics_ICA.coh_missclassifications.(band) channel_metrics_ICA.coh_confmat.(band) channel_metrics_ICA.coh_confmat_metrics.(band) channel_metrics_ICA.coh_error_network.(band)] = compute_network_variance(channel_metrics_clean.coh_network.(band),channel_metrics_ICA.coh_network.(band));
% [channel_metrics_ICA.coh_missclassifications_calm.(band) channel_metrics_ICA.coh_confmat_calm.(band) channel_metrics_ICA.coh_confmat_calm_metrics.(band) channel_metrics_ICA.coh_error_network_calm.(band)] = compute_network_variance(channel_metrics_clean.coh_network_calm.(band),channel_metrics_ICA.coh_network_calm.(band));
% [channel_metrics_ICA.coh_missclassifications_act.(band) channel_metrics_ICA.coh_confmat_act.(band) channel_metrics_ICA.coh_confmat_act_metrics.(band) channel_metrics_ICA.coh_error_network_act.(band)] = compute_network_variance(channel_metrics_clean.coh_network_act.(band),channel_metrics_ICA.coh_network_act.(band));
end

keep = noiseindex; %change this manually to the kept component index
icaEEG4 = zeros(size(icaEEG)); %zero out all ICA components
icaEEG4(keep,:) = icaEEG(keep,:); %restore eyeblink component
%icaEEG4(keep,:) = strong_comp; 

single_comp_data = A*icaEEG4; %calculate single-component channel data


% ICA_eyeblink_proj = A(:,noiseindex);
% ICA_eyeblink_proj = ICA_eyeblink_proj'/max(abs(ICA_eyeblink_proj));
ICA_eyeblink_proj = single_comp_data(:,act_start+40) ;
ICA_eyeblink_proj = ICA_eyeblink_proj/max(abs(ICA_eyeblink_proj));
if ICA_eyeblink_proj(33) * projection_matrix(33) < 0
    ICA_eyeblink_proj = -ICA_eyeblink_proj;
end


true_proj_norm = projection_matrix'/max(abs(projection_matrix));


local_error = true_proj_norm - ICA_eyeblink_proj;
channel_metrics_ICA.localization_error = local_error;
toc
%% wICA artifact rejection (include wICA folder)
% NOTE: For better artifact suppression, provide manually the numbers 
% of components to be processed. You can also tune the other arguments
% for your data set.
%

nICs = eyeblink_idx; % Components to be processed, e.g. [1, 4:7]
Kthr = 1.15;             % Tolerance for cleaning artifacts, try: 1, 1.15,...
ArtefThreshold = 3;      % Threshold for detection of ICs with artefacts
                         % Set lower values if you manually select ICs with 
                         % artifacts. Otherwise increase
verbose = 'on';          % print some intermediate results                         
icaEEG2 = RemoveStrongArtifacts(icaEEG, nICs, Kthr, ArtefThreshold, Fs, verbose);

strong_comp = icaEEG(nICs,:) - icaEEG2(nICs,:);
figure('Name','wICA noise component')
plot(strong_comp')
figure
plotmts2(icaEEG2(1:10,:)')
EEG_wICA = A*icaEEG2; %cleaned channel data
%%

EEG.data = EEG_wICA;
asec = pop_rsbl(EEG, saveFull, account4artifacts, src2roiReductionType, solverType, updateFreq);

%all 68 source
ROI_acts = asec.etc.src.act;
%ROI_acts = ROI_acts(active_regions,:);

ROI_acts = ROI_acts/dipoles_per_region;

channel_metrics_wICA = Compute_channel_metrics(EEG_wICA,EEG_clean,calm_start,calm_end,act_start,act_end);

[channel_metrics_wICA.RCGCIM_missclassifications channel_metrics_wICA.RCGCIM_confmat channel_metrics_wICA.RCGCIM_confmat_metrics channel_metrics_wICA.RCGCIM_error_network] = compute_network_variance(channel_metrics_clean.RCGCIM_network ,channel_metrics_wICA.RCGCIM_network);
[channel_metrics_wICA.RCGCIM_missclassifications_calm channel_metrics_wICA.RCGCIM_confmat_calm channel_metrics_wICA.RCGCIM_confmat_calm_metrics channel_metrics_wICA.RCGCIM_error_network_calm] = compute_network_variance(channel_metrics_clean.RCGCIM_network_calm,channel_metrics_wICA.RCGCIM_network_calm);
[channel_metrics_wICA.RCGCIM_missclassifications_act channel_metrics_wICA.RCGCIM_confmat_act channel_metrics_wICA.RCGCIM_confmat_act_metrics channel_metrics_wICA.RCGCIM_error_network_act] = compute_network_variance(channel_metrics_clean.RCGCIM_network_act,channel_metrics_wICA.RCGCIM_network_act);


source_metrics_wICA_ext = compute_source_metrics_ext(ROI_acts,Fs,calm_start,calm_end,act_start,act_end);

% [source_metrics_wICA_ext.CGCIM_missclassifications source_metrics_wICA_ext.CGCIM_confmat source_metrics_wICA_ext.CGCIM_confmat_metrics] = compute_network_variance(Adj_mat_ext ,source_metrics_wICA_ext.CGCIM_network);
% [source_metrics_wICA_ext.CGCIM_missclassifications_calm source_metrics_wICA_ext.CGCIM_confmat_calm source_metrics_wICA_ext.CGCIM_confmat_calm_metrics] = compute_network_variance(Adj_mat_ext,source_metrics_wICA_ext.CGCIM_network_calm);
% [source_metrics_wICA_ext.CGCIM_missclassifications_act source_metrics_wICA_ext.CGCIM_confmat_act source_metrics_wICA_ext.CGCIM_confmat_act_metrics] = compute_network_variance(Adj_mat_ext,source_metrics_wICA_ext.CGCIM_network_act);
[source_metrics_wICA_ext.RCGCIM_missclassifications source_metrics_wICA_ext.RCGCIM_confmat source_metrics_wICA_ext.RCGCIM_confmat_metrics source_metrics_wICA_ext.RCGCIM_error_network] = compute_network_variance(Adj_mat_ext ,source_metrics_wICA_ext.RCGCIM_network);
[source_metrics_wICA_ext.RCGCIM_missclassifications_calm source_metrics_wICA_ext.RCGCIM_confmat_calm source_metrics_wICA_ext.RCGCIM_confmat_calm_metrics source_metrics_wICA_ext.RCGCIM_error_network_calm] = compute_network_variance(Adj_mat_ext,source_metrics_wICA_ext.RCGCIM_network_calm);
[source_metrics_wICA_ext.RCGCIM_missclassifications_act source_metrics_wICA_ext.RCGCIM_confmat_act source_metrics_wICA_ext.RCGCIM_confmat_act_metrics source_metrics_wICA_ext.RCGCIM_error_network_act] = compute_network_variance(Adj_mat_ext,source_metrics_wICA_ext.RCGCIM_network_act);



for j=1:size(frequency_bands,1)
band = frequency_bands(j);
band = string(band);

[source_metrics_wICA_ext.DTF_missclassifications_act.(band) source_metrics_wICA_ext.DTF_confmat_act.(band) source_metrics_wICA_ext.DTF_confmat_act_metrics.(band) source_metrics_wICA_ext.DTF_error_network_act.(band)] = compute_network_variance(Adj_mat_ext,source_metrics_wICA_ext.DTF_network_act.(band));

[source_metrics_wICA_ext.PDC_missclassifications_act.(band) source_metrics_wICA_ext.PDC_confmat_act.(band) source_metrics_wICA_ext.PDC_confmat_act_metrics.(band) source_metrics_wICA_ext.PDC_error_network_act.(band)] = compute_network_variance(Adj_mat_ext,source_metrics_wICA_ext.PDC_network_act.(band));

[channel_metrics_wICA.coh_missclassifications.(band) channel_metrics_wICA.coh_confmat.(band) channel_metrics_wICA.coh_confmat_metrics.(band) channel_metrics_wICA.coh_error_network.(band)] = compute_network_variance(channel_metrics_clean.coh_network.(band),channel_metrics_wICA.coh_network.(band));
[channel_metrics_wICA.coh_missclassifications_calm.(band) channel_metrics_wICA.coh_confmat_calm.(band) channel_metrics_wICA.coh_confmat_calm_metrics.(band) channel_metrics_wICA.coh_error_network_calm.(band)] = compute_network_variance(channel_metrics_clean.coh_network_calm.(band),channel_metrics_wICA.coh_network_calm.(band));
[channel_metrics_wICA.coh_missclassifications_act.(band) channel_metrics_wICA.coh_confmat_act.(band) channel_metrics_wICA.coh_confmat_act_metric.(band) channel_metrics_wICA.coh_error_network_act.(band)] = compute_network_variance(channel_metrics_clean.coh_network_act.(band),channel_metrics_wICA.coh_network_act.(band));
end


% 7 sources
ROI_acts = asec.etc.src.act;

ROI_acts = ROI_acts(active_regions,:);
ROI_acts = ROI_acts/dipoles_per_region;

source_metrics_wICA = compute_source_metrics(ROI_acts,xM',Fs,calm_start,calm_end,act_start,act_end);

[source_metrics_wICA.CGCIM_missclassifications source_metrics_wICA.CGCIM_confmat source_metrics_wICA.CGCIM_confmat_metrics source_metrics_wICA.CGCIM_error_network] = compute_network_variance(Adj_mat,source_metrics_wICA.CGCIM_network);
[source_metrics_wICA.CGCIM_missclassifications_calm source_metrics_wICA.CGCIM_confmat_calm source_metrics_wICA.CGCIM_confmat_calm_metrics source_metrics_wICA.CGCIM_error_network_calm] = compute_network_variance(Adj_mat,source_metrics_wICA.CGCIM_network_calm);
[source_metrics_wICA.CGCIM_missclassifications_act source_metrics_wICA.CGCIM_confmat_act source_metrics_wICA.CGCIM_confmat_act_metrics source_metrics_wICA.CGCIM_error_network_act] = compute_network_variance(Adj_mat,source_metrics_wICA.CGCIM_network_act);


for i=1:size(frequency_bands,1)
band = frequency_bands(i);
band = string(band);
[source_metrics_wICA.DTF_missclassifications.(band) source_metrics_wICA.DTF_confmat.(band) source_metrics_wICA.DTF_confmat_metrics.(band) source_metrics_wICA.DTF_error_network.(band)] = compute_network_variance(Adj_mat,source_metrics_wICA.DTF_network.(band));
[source_metrics_wICA.DTF_missclassifications_calm.(band) source_metrics_wICA.DTF_confmat_calm.(band) source_metrics_wICA.DTF_confmat_calm_metrics.(band) source_metrics_wICA.DTF_error_network_calm.(band)] = compute_network_variance(Adj_mat,source_metrics_wICA.DTF_network_calm.(band));
[source_metrics_wICA.DTF_missclassifications_act.(band) source_metrics_wICA.DTF_confmat_act.(band) source_metrics_wICA.DTF_confmat_act_metrics.(band) source_metrics_wICA.DTF_error_network_act.(band)] = compute_network_variance(Adj_mat,source_metrics_wICA.DTF_network_act.(band));

[source_metrics_wICA.PDC_missclassifications.(band) source_metrics_wICA.PDC_confmat.(band) source_metrics_wICA.PDC_confmat_metrics.(band) source_metrics_wICA.PDC_error_network.(band)] = compute_network_variance(Adj_mat,source_metrics_wICA.PDC_network.(band));
[source_metrics_wICA.PDC_missclassifications_calm.(band) source_metrics_wICA.PDC_confmat_calm.(band) source_metrics_wICA.PDC_confmat_calm_metrics.(band) source_metrics_wICA.PDC_error_network_calm.(band)] = compute_network_variance(Adj_mat,source_metrics_wICA.PDC_network_calm.(band));
[source_metrics_wICA.PDC_missclassifications_act.(band) source_metrics_wICA.PDC_confmat_act.(band) source_metrics_wICA.PDC_confmat_act_metrics.(band) source_metrics_wICA.PDC_error_network_act.(band)] = compute_network_variance(Adj_mat,source_metrics_wICA.PDC_network_act.(band));

[channel_metrics_wICA.coh_missclassifications.(band) channel_metrics_wICA.coh_confmat.(band) channel_metrics_wICA.coh_confmat_metrics.(band) channel_metrics_wICA.coh_error_network.(band)] = compute_network_variance(channel_metrics_clean.coh_network.(band),channel_metrics_wICA.coh_network.(band));
[channel_metrics_wICA.coh_missclassifications_calm.(band) channel_metrics_wICA.coh_confmat_calm.(band) channel_metrics_wICA.coh_confmat_calm_metrics.(band) channel_metrics_wICA.coh_error_network_calm.(band)] = compute_network_variance(channel_metrics_clean.coh_network_calm.(band),channel_metrics_wICA.coh_network_calm.(band));
[channel_metrics_wICA.coh_missclassifications_act.(band) channel_metrics_wICA.coh_confmat_act.(band) channel_metrics_wICA.coh_confmat_act_metrics.(band) channel_metrics_wICA.coh_error_network_act.(band)] = compute_network_variance(channel_metrics_clean.coh_network_act.(band),channel_metrics_wICA.coh_network_act.(band));
end

%% Strip by strip ICA

strip.a = [1 33 34 2 3 37 36 35]; %front of the head strip
strip.b = [4:11 38:47];
strip.c = [12:19 32 48:56];
strip.d = [20:26 30 31 57:63];
strip.e = [27 28 29 64]; %back of the head strip
names = ['a' 'b' 'c' 'd' 'e'];
eyeblink_noise_strips = zeros(size(EEG_clean));
for i=1:5
    current_strip = names(i);
    current_strip = strip.(current_strip);
    EEG_strip = EEG_noisy(current_strip,:);
    [icaEEG_strip, A_strip, W_strip] = fastica(EEG_strip,'stabilization','on','verbose','on'); 
    [eyeblink_idx_strip best_match_strip] = find_eyeblink_comp(sim_eyeblink,icaEEG_strip);
    noiseindex_strip = best_match_strip;
    eyeblink_noise_strips(current_strip,:) = A_strip(:,noiseindex_strip)*icaEEG_strip(noiseindex_strip,:);
    figure
    plotmts2(icaEEG_strip')
end
    
EEG_strips = EEG_noisy - eyeblink_noise_strips;
%calculate locality from a single sample during eyeblink
strips_local = eyeblink_noise_strips(:,act+40)/max(abs(eyeblink_noise_strips(:,act+40)));
local_error = true_proj_norm - strips_local';
figure
plotmts2(eyeblink_noise_strips(1:15,:)')
figure
plotmts2(EEG_noisy(1:15,:)')
%% ICA strips source reconstruction

EEG.data = EEG_strips;
asec = pop_rsbl(EEG, saveFull, account4artifacts, src2roiReductionType, solverType, updateFreq);

%all 68 source
ROI_acts = asec.etc.src.act;
%ROI_acts = ROI_acts(active_regions,:);

ROI_acts = ROI_acts/dipoles_per_region;

channel_metrics_ICA_strips = Compute_channel_metrics(EEG_strips,EEG_clean,calm_start,calm_end,act_start,act_end);

[channel_metrics_ICA_strips.RCGCIM_missclassifications channel_metrics_ICA_strips.RCGCIM_confmat channel_metrics_ICA_strips.RCGCIM_confmat_metrics channel_metrics_ICA_strips.RCGCIM_error_network] = compute_network_variance(channel_metrics_clean.RCGCIM_network ,channel_metrics_ICA_strips.RCGCIM_network);
[channel_metrics_ICA_strips.RCGCIM_missclassifications_calm channel_metrics_ICA_strips.RCGCIM_confmat_calm channel_metrics_ICA_strips.RCGCIM_confmat_calm_metrics channel_metrics_ICA_strips.RCGCIM_error_network_calm] = compute_network_variance(channel_metrics_clean.RCGCIM_network_calm,channel_metrics_ICA_strips.RCGCIM_network_calm);
[channel_metrics_ICA_strips.RCGCIM_missclassifications_act channel_metrics_ICA_strips.RCGCIM_confmat_act channel_metrics_ICA_strips.RCGCIM_confmat_act_metrics channel_metrics_ICA_strips.RCGCIM_error_network_act] = compute_network_variance(channel_metrics_clean.RCGCIM_network_act,channel_metrics_ICA_strips.RCGCIM_network_act);


source_metrics_ICA_strips_ext = compute_source_metrics_ext(ROI_acts,Fs,calm_start,calm_end,act_start,act_end);

% [source_metrics_ICA_strips_ext.CGCIM_missclassifications source_metrics_ICA_strips_ext.CGCIM_confmat source_metrics_ICA_strips_ext.CGCIM_confmat_metrics] = compute_network_variance(Adj_mat_ext ,source_metrics_ICA_strips_ext.CGCIM_network);
% [source_metrics_ICA_strips_ext.CGCIM_missclassifications_calm source_metrics_ICA_strips_ext.CGCIM_confmat_calm source_metrics_ICA_strips_ext.CGCIM_confmat_calm_metrics] = compute_network_variance(Adj_mat_ext,source_metrics_ICA_strips_ext.CGCIM_network_calm);
% [source_metrics_ICA_strips_ext.CGCIM_missclassifications_act source_metrics_ICA_strips_ext.CGCIM_confmat_act source_metrics_ICA_strips_ext.CGCIM_confmat_act_metrics] = compute_network_variance(Adj_mat_ext,source_metrics_ICA_strips_ext.CGCIM_network_act);
[source_metrics_ICA_strips_ext.RCGCIM_missclassifications source_metrics_ICA_strips_ext.RCGCIM_confmat source_metrics_ICA_strips_ext.RCGCIM_confmat_metrics source_metrics_ICA_strips_ext.RCGCIM_error_network] = compute_network_variance(Adj_mat_ext ,source_metrics_ICA_strips_ext.RCGCIM_network);
[source_metrics_ICA_strips_ext.RCGCIM_missclassifications_calm source_metrics_ICA_strips_ext.RCGCIM_confmat_calm source_metrics_ICA_strips_ext.RCGCIM_confmat_calm_metricssource_metrics_ICA_strips_ext.RCGCIM_error_network_calm] = compute_network_variance(Adj_mat_ext,source_metrics_ICA_strips_ext.RCGCIM_network_calm);
[source_metrics_ICA_strips_ext.RCGCIM_missclassifications_act source_metrics_ICA_strips_ext.RCGCIM_confmat_act source_metrics_ICA_strips_ext.RCGCIM_confmat_act_metrics source_metrics_ICA_strips_ext.RCGCIM_error_network_act] = compute_network_variance(Adj_mat_ext,source_metrics_ICA_strips_ext.RCGCIM_network_act);



for j=1:size(frequency_bands,1)
band = frequency_bands(j);
band = string(band);

[source_metrics_ICA_strips_ext.DTF_missclassifications_act.(band) source_metrics_ICA_strips_ext.DTF_confmat_act.(band) source_metrics_ICA_strips_ext.DTF_confmat_act_metrics.(band) source_metrics_ICA_strips_ext.DTF_error_network_act.(band)] = compute_network_variance(Adj_mat_ext,source_metrics_ICA_strips_ext.DTF_network_act.(band));

[source_metrics_ICA_strips_ext.PDC_missclassifications_act.(band) source_metrics_ICA_strips_ext.PDC_confmat_act.(band) source_metrics_ICA_strips_ext.PDC_confmat_act_metrics.(band) source_metrics_ICA_strips_ext.PDC_error_netowork_act.(band)] = compute_network_variance(Adj_mat_ext,source_metrics_ICA_strips_ext.PDC_network_act.(band));

[channel_metrics_ICA_strips.coh_missclassifications.(band) channel_metrics_ICA_strips.coh_confmat.(band) channel_metrics_ICA_strips.coh_confmat_metrics.(band) channel_metrics_ICA_strips.coh_error_network.(band)] = compute_network_variance(channel_metrics_clean.coh_network.(band),channel_metrics_ICA_strips.coh_network.(band));
[channel_metrics_ICA_strips.coh_missclassifications_calm.(band) channel_metrics_ICA_strips.coh_confmat_calm.(band) channel_metrics_ICA_strips.coh_confmat_calm_metrics.(band) channel_metrics_ICA_strips.coh_error_network_calm.(band)] = compute_network_variance(channel_metrics_clean.coh_network_calm.(band),channel_metrics_ICA_strips.coh_network_calm.(band));
[channel_metrics_ICA_strips.coh_missclassifications_act.(band) channel_metrics_ICA_strips.coh_confmat_act.(band) channel_metrics_ICA_strips.coh_confmat_act_metric.(band) channel_metrics_ICA_strips.coh_error_network_act.(band)] = compute_network_variance(channel_metrics_clean.coh_network_act.(band),channel_metrics_ICA_strips.coh_network_act.(band));
end


% 7 sources
ROI_acts = asec.etc.src.act;

ROI_acts = ROI_acts(active_regions,:);
ROI_acts = ROI_acts/dipoles_per_region;

source_metrics_ICA_strips = compute_source_metrics(ROI_acts,xM',Fs,calm_start,calm_end,act_start,act_end);

[source_metrics_ICA_strips.CGCIM_missclassifications source_metrics_ICA_strips.CGCIM_confmat source_metrics_ICA_strips.CGCIM_confmat_metrics source_metrics_ICA_strips.CGCIM_error_network] = compute_network_variance(Adj_mat,source_metrics_ICA_strips.CGCIM_network);
[source_metrics_ICA_strips.CGCIM_missclassifications_calm source_metrics_ICA_strips.CGCIM_confmat_calm source_metrics_ICA_strips.CGCIM_confmat_calm_metrics source_metrics_ICA_strips.CGCIM_error_network_calm] = compute_network_variance(Adj_mat,source_metrics_ICA_strips.CGCIM_network_calm);
[source_metrics_ICA_strips.CGCIM_missclassifications_act source_metrics_ICA_strips.CGCIM_confmat_act source_metrics_ICA_strips.CGCIM_confmat_act_metrics source_metrics_ICA_strips.CGCIM_error_network_act] = compute_network_variance(Adj_mat,source_metrics_ICA_strips.CGCIM_network_act);


for i=1:size(frequency_bands,1)
band = frequency_bands(i);
band = string(band);
[source_metrics_ICA_strips.DTF_missclassifications.(band) source_metrics_ICA_strips.DTF_confmat.(band) source_metrics_ICA_strips.DTF_confmat_metrics.(band) source_metrics_ICA_strips.DTF_error_network.(band)] = compute_network_variance(Adj_mat,source_metrics_ICA_strips.DTF_network.(band));
[source_metrics_ICA_strips.DTF_missclassifications_calm.(band) source_metrics_ICA_strips.DTF_confmat_calm.(band) source_metrics_ICA_strips.DTF_confmat_calm_metrics.(band) source_metrics_ICA_strips.DTF_error_network_calm.(band)] = compute_network_variance(Adj_mat,source_metrics_ICA_strips.DTF_network_calm.(band));
[source_metrics_ICA_strips.DTF_missclassifications_act.(band) source_metrics_ICA_strips.DTF_confmat_act.(band) source_metrics_ICA_strips.DTF_confmat_act_metrics.(band) source_metrics_ICA_strips.DTF_error_network_act.(band)] = compute_network_variance(Adj_mat,source_metrics_ICA_strips.DTF_network_act.(band));

[source_metrics_ICA_strips.PDC_missclassifications.(band) source_metrics_ICA_strips.PDC_confmat.(band) source_metrics_ICA_strips.PDC_confmat_metrics.(band) source_metrics_ICA_strips.PDC_error_network.(band)] = compute_network_variance(Adj_mat,source_metrics_ICA_strips.PDC_network.(band));
[source_metrics_ICA_strips.PDC_missclassifications_calm.(band) source_metrics_ICA_strips.PDC_confmat_calm.(band) source_metrics_ICA_strips.PDC_confmat_calm_metrics.(band) source_metrics_ICA_strips.PDC_error_network_calm.(band)] = compute_network_variance(Adj_mat,source_metrics_ICA_strips.PDC_network_calm.(band));
[source_metrics_ICA_strips.PDC_missclassifications_act.(band) source_metrics_ICA_strips.PDC_confmat_act.(band) source_metrics_ICA_strips.PDC_confmat_act_metrics.(band) source_metrics_ICA_strips.PDC_error_network_act.(band)] = compute_network_variance(Adj_mat,source_metrics_ICA_strips.PDC_network_act.(band));

[channel_metrics_ICA_strips.coh_missclassifications.(band) channel_metrics_wICA.coh_confmat.(band) channel_metrics_ICA_strips.coh_confmat_metrics.(band) channel_metrics_ICA_strips.coh_error_network.(band)] = compute_network_variance(channel_metrics_clean.coh_network.(band),channel_metrics_ICA_strips.coh_network.(band));
[channel_metrics_ICA_strips.coh_missclassifications_calm.(band) channel_metrics_ICA_strips.coh_confmat_calm.(band) channel_metrics_ICA_strips.coh_confmat_calm_metrics.(band) channel_metrics_ICA_strips.coh_error_network_calm.(band)] = compute_network_variance(channel_metrics_clean.coh_network_calm.(band),channel_metrics_ICA_strips.coh_network_calm.(band));
[channel_metrics_ICA_strips.coh_missclassifications_act.(band) channel_metrics_ICA_strips.coh_confmat_act.(band) channel_metrics_ICA_strips.coh_confmat_act_metrics.(band) channel_metrics_ICA_strips.coh_error_network_act.(band)] = compute_network_variance(channel_metrics_clean.coh_network_act.(band),channel_metrics_ICA_strips.coh_network_act.(band));
end

channel_metrics_ICA_strips.localization_error = local_error;
%% Rebuild eyeblink B


% current_strip = strip.a; %frontal strip
% eyeblink_noise_2 = zeros(size(EEG_clean));
% reference = 33;
% reference_index = find(current_strip == 33); % index (in strip) of electrode 33
% EEG_strip = EEG_noisy(current_strip,:);
% [icaEEG_strip, A_strip, W_strip] = fastica(EEG_strip,'stabilization','on','verbose','on'); 
% [eyeblink_idx_strip best_match_strip] = find_eyeblink_comp(sim_eyeblink,icaEEG_strip);
% noiseindex_strip = best_match_strip;
% eyeblink_noise_strips(current_strip,:) = A_strip(:,noiseindex_strip)*icaEEG_strip(noiseindex_strip,:);


reference_idx = 33; % a frontal electrode
frontal_electrodes = [1 33 34 2 3 37 36 35];
eyeblink_noise_reference = eyeblink_noise_strips(reference_idx,:); %get reference noise from ICA_strips, not ICA

eyeblink_dev = 0.3;
data_frontal = EEG_noisy(reference_idx,:);
[distances indexes] = find_peaks_ts(data_frontal,eyeblink_dev); %find eyeblinks

%find average eyeblink "peak to peak" distance for all electrodes.
%%%%%%%%%%%%%%%%%%%%%%%
window_mean = 3;  %%%%% parameters
window = 40;      %%%%%
%%%%%%%%%%%%%%%%%%%%%%%

dist = zeros(size(EEG_noisy,1),size(indexes,2));
for ch=1:size(EEG_noisy,1)
    for i = 1:size(indexes,2)
        dist(ch,i) = mean(EEG_noisy(ch,indexes(i)-window_mean:indexes(i)+window_mean)) - mean(EEG_noisy(ch,indexes(i)-window_mean-window:indexes(i)-window+window_mean));
    end
end
mean_ch_dist = mean(dist');
mean_ch_dist_norm = mean_ch_dist/max(abs(mean_ch_dist));

%definitely same polarity, no need to check for mean_ch_dist_norm sign
local_error = true_proj_norm - mean_ch_dist_norm';

eyeblink_noise_2 = mean_ch_dist_norm'*eyeblink_noise_reference ;
EEG_A_rebuilt = EEG_noisy - eyeblink_noise_2;
% activity = zeros(n,1);
% for i = 1:length(indexes)
%     activity(indexes(i)-50:indexes(i)+150) = 1;
% end
% figure
% plot(activity)
% activity = logical(activity);
% EEG_affected = EEG_noisy(:,activity);
% EEG_affected = [EEG_affected EEG_affected EEG_affected EEG_affected EEG_affected];
% figure
% plotmts2(EEG_affected(1:10,:)')
% [icaEEG_affected, A_affected, W_affected] = fastica(EEG_affected,'stabilization','on','verbose','on'); 
% [eyeblink_idx_affected best_match_affected] = find_eyeblink_comp([sim_eyeblink(activity) sim_eyeblink(activity)],icaEEG_affected);
% noiseindex_affected = best_match_affected %change this manually to the noise component index
% figure
% plotmts2(icaEEG_affected(1:12,:)')
% A_correction = A_affected(:,noiseindex_affected);
% 
% 
% A_correction = A_correction/max(abs(A_correction));
% A_reference = A(33,noiseindex); %%any frontal electrode will do - assuming ICA removes a decent % of noise at that electrode
% if A_correction(33)*A_reference < 0 % ??
%     A_correction = - A_correction;
% end
% A_rebuilt = A;
% A_rebuilt(:,noiseindex) = A_reference*A_correction;
% A_rebuilt(frontal_electrodes,noiseindex) = A(frontal_electrodes,noiseindex); % don't make any changes to frontal electrodes
% %A_rebuilt(:,noiseindex) = (A_rebuilt(:,noiseindex) + A(:,noiseindex))/2;

%% A rebuilt source reconstruction

EEG.data = EEG_A_rebuilt;
asec = pop_rsbl(EEG, saveFull, account4artifacts, src2roiReductionType, solverType, updateFreq);

%all 68 source
ROI_acts = asec.etc.src.act;
%ROI_acts = ROI_acts(active_regions,:);

ROI_acts = ROI_acts/dipoles_per_region;

channel_metrics_ICA_A_rebuilt = Compute_channel_metrics(EEG_strips,EEG_clean,calm_start,calm_end,act_start,act_end);

[channel_metrics_ICA_A_rebuilt.RCGCIM_missclassifications channel_metrics_ICA_A_rebuilt.RCGCIM_confmat channel_metrics_ICA_A_rebuilt.RCGCIM_confmat_metrics channel_metrics_ICA_A_rebuilt.RCGCIM_error_network] = compute_network_variance(channel_metrics_clean.RCGCIM_network ,channel_metrics_ICA_A_rebuilt.RCGCIM_network);
[channel_metrics_ICA_A_rebuilt.RCGCIM_missclassifications_calm channel_metrics_ICA_A_rebuilt.RCGCIM_confmat_calm channel_metrics_ICA_A_rebuilt.RCGCIM_confmat_calm_metrics channel_metrics_ICA_A_rebuilt.RCGCIM_error_network_calm] = compute_network_variance(channel_metrics_clean.RCGCIM_network_calm,channel_metrics_ICA_A_rebuilt.RCGCIM_network_calm);
[channel_metrics_ICA_A_rebuilt.RCGCIM_missclassifications_act channel_metrics_ICA_A_rebuilt.RCGCIM_confmat_act channel_metrics_ICA_A_rebuilt.RCGCIM_confmat_act_metrics channel_metrics_ICA_A_rebuilt.RCGCIM_error_network_act] = compute_network_variance(channel_metrics_clean.RCGCIM_network_act,channel_metrics_ICA_A_rebuilt.RCGCIM_network_act);


source_metrics_ICA_A_rebuilt_ext = compute_source_metrics_ext(ROI_acts,Fs,calm_start,calm_end,act_start,act_end);

% [source_metrics_ICA_strips_ext.CGCIM_missclassifications source_metrics_ICA_strips_ext.CGCIM_confmat source_metrics_ICA_strips_ext.CGCIM_confmat_metrics] = compute_network_variance(Adj_mat_ext ,source_metrics_ICA_strips_ext.CGCIM_network);
% [source_metrics_ICA_strips_ext.CGCIM_missclassifications_calm source_metrics_ICA_strips_ext.CGCIM_confmat_calm source_metrics_ICA_strips_ext.CGCIM_confmat_calm_metrics] = compute_network_variance(Adj_mat_ext,source_metrics_ICA_strips_ext.CGCIM_network_calm);
% [source_metrics_ICA_strips_ext.CGCIM_missclassifications_act source_metrics_ICA_strips_ext.CGCIM_confmat_act source_metrics_ICA_strips_ext.CGCIM_confmat_act_metrics] = compute_network_variance(Adj_mat_ext,source_metrics_ICA_strips_ext.CGCIM_network_act);
[source_metrics_ICA_A_rebuilt_ext.RCGCIM_missclassifications source_metrics_ICA_A_rebuilt_ext.RCGCIM_confmat source_metrics_ICA_A_rebuilt_ext.RCGCIM_confmat_metrics source_metrics_ICA_A_rebuilt_ext.RCGCIM_error_network] = compute_network_variance(Adj_mat_ext ,source_metrics_ICA_A_rebuilt_ext.RCGCIM_network);
[source_metrics_ICA_A_rebuilt_ext.RCGCIM_missclassifications_calm source_metrics_ICA_A_rebuilt_ext.RCGCIM_confmat_calm source_metrics_ICA_A_rebuilt_ext.RCGCIM_confmat_calm_metricssource_metrics_ICA_strips_ext.RCGCIM_error_network_calm] = compute_network_variance(Adj_mat_ext,source_metrics_ICA_A_rebuilt_ext.RCGCIM_network_calm);
[source_metrics_ICA_A_rebuilt_ext.RCGCIM_missclassifications_act source_metrics_ICA_A_rebuilt_ext.RCGCIM_confmat_act source_metrics_ICA_A_rebuilt_ext.RCGCIM_confmat_act_metrics source_metrics_ICA_A_rebuilt_ext.RCGCIM_error_network_act] = compute_network_variance(Adj_mat_ext,source_metrics_ICA_A_rebuilt_ext.RCGCIM_network_act);



for j=1:size(frequency_bands,1)
band = frequency_bands(j);
band = string(band);

[source_metrics_ICA_A_rebuilt_ext.DTF_missclassifications_act.(band) source_metrics_ICA_A_rebuilt_ext.DTF_confmat_act.(band) source_metrics_ICA_A_rebuilt_ext.DTF_confmat_act_metrics.(band) source_metrics_ICA_A_rebuilt_ext.DTF_error_network_act.(band)] = compute_network_variance(Adj_mat_ext,source_metrics_ICA_A_rebuilt_ext.DTF_network_act.(band));

[source_metrics_ICA_A_rebuilt_ext.PDC_missclassifications_act.(band) source_metrics_ICA_A_rebuilt_ext.PDC_confmat_act.(band) source_metrics_ICA_A_rebuilt_ext.PDC_confmat_act_metrics.(band) source_metrics_ICA_A_rebuilt_ext.PDC_error_netowork_act.(band)] = compute_network_variance(Adj_mat_ext,source_metrics_ICA_A_rebuilt_ext.PDC_network_act.(band));

[channel_metrics_ICA_A_rebuilt.coh_missclassifications.(band) channel_metrics_ICA_A_rebuilt.coh_confmat.(band) channel_metrics_ICA_A_rebuilt.coh_confmat_metrics.(band) channel_metrics_ICA_A_rebuilt.coh_error_network.(band)] = compute_network_variance(channel_metrics_clean.coh_network.(band),channel_metrics_ICA_A_rebuilt.coh_network.(band));
[channel_metrics_ICA_A_rebuilt.coh_missclassifications_calm.(band) channel_metrics_ICA_A_rebuilt.coh_confmat_calm.(band) channel_metrics_ICA_A_rebuilt.coh_confmat_calm_metrics.(band) channel_metrics_ICA_A_rebuilt.coh_error_network_calm.(band)] = compute_network_variance(channel_metrics_clean.coh_network_calm.(band),channel_metrics_ICA_A_rebuilt.coh_network_calm.(band));
[channel_metrics_ICA_A_rebuilt.coh_missclassifications_act.(band) channel_metrics_ICA_A_rebuilt.coh_confmat_act.(band) channel_metrics_ICA_A_rebuilt.coh_confmat_act_metric.(band) channel_metrics_ICA_A_rebuilt.coh_error_network_act.(band)] = compute_network_variance(channel_metrics_clean.coh_network_act.(band),channel_metrics_ICA_A_rebuilt.coh_network_act.(band));
end


% 7 sources
ROI_acts = asec.etc.src.act;

ROI_acts = ROI_acts(active_regions,:);
ROI_acts = ROI_acts/dipoles_per_region;

source_metrics_ICA_A_rebuilt = compute_source_metrics(ROI_acts,xM',Fs,calm_start,calm_end,act_start,act_end);

[source_metrics_ICA_A_rebuilt.CGCIM_missclassifications source_metrics_ICA_A_rebuilt.CGCIM_confmat source_metrics_ICA_A_rebuilt.CGCIM_confmat_metrics source_metrics_ICA_A_rebuilt.CGCIM_error_network] = compute_network_variance(Adj_mat,source_metrics_ICA_A_rebuilt.CGCIM_network);
[source_metrics_ICA_A_rebuilt.CGCIM_missclassifications_calm source_metrics_ICA_A_rebuilt.CGCIM_confmat_calm source_metrics_ICA_A_rebuilt.CGCIM_confmat_calm_metrics source_metrics_ICA_A_rebuilt.CGCIM_error_network_calm] = compute_network_variance(Adj_mat,source_metrics_ICA_A_rebuilt.CGCIM_network_calm);
[source_metrics_ICA_A_rebuilt.CGCIM_missclassifications_act source_metrics_ICA_A_rebuilt.CGCIM_confmat_act source_metrics_ICA_A_rebuilt.CGCIM_confmat_act_metrics source_metrics_ICA_A_rebuilt.CGCIM_error_network_act] = compute_network_variance(Adj_mat,source_metrics_ICA_A_rebuilt.CGCIM_network_act);


for i=1:size(frequency_bands,1)
band = frequency_bands(i);
band = string(band);
[source_metrics_ICA_A_rebuilt.DTF_missclassifications.(band) source_metrics_ICA_A_rebuilt.DTF_confmat.(band) source_metrics_ICA_A_rebuilt.DTF_confmat_metrics.(band) source_metrics_ICA_A_rebuilt.DTF_error_network.(band)] = compute_network_variance(Adj_mat,source_metrics_ICA_A_rebuilt.DTF_network.(band));
[source_metrics_ICA_A_rebuilt.DTF_missclassifications_calm.(band) source_metrics_ICA_A_rebuilt.DTF_confmat_calm.(band) source_metrics_ICA_A_rebuilt.DTF_confmat_calm_metrics.(band) source_metrics_ICA_A_rebuilt.DTF_error_network_calm.(band)] = compute_network_variance(Adj_mat,source_metrics_ICA_A_rebuilt.DTF_network_calm.(band));
[source_metrics_ICA_A_rebuilt.DTF_missclassifications_act.(band) source_metrics_ICA_A_rebuilt.DTF_confmat_act.(band) source_metrics_ICA_A_rebuilt.DTF_confmat_act_metrics.(band) source_metrics_ICA_A_rebuilt.DTF_error_network_act.(band)] = compute_network_variance(Adj_mat,source_metrics_ICA_A_rebuilt.DTF_network_act.(band));

[source_metrics_ICA_A_rebuilt.PDC_missclassifications.(band) source_metrics_ICA_A_rebuilt.PDC_confmat.(band) source_metrics_ICA_A_rebuilt.PDC_confmat_metrics.(band) source_metrics_ICA_A_rebuilt.PDC_error_network.(band)] = compute_network_variance(Adj_mat,source_metrics_ICA_A_rebuilt.PDC_network.(band));
[source_metrics_ICA_A_rebuilt.PDC_missclassifications_calm.(band) source_metrics_ICA_A_rebuilt.PDC_confmat_calm.(band) source_metrics_ICA_A_rebuilt.PDC_confmat_calm_metrics.(band) source_metrics_ICA_A_rebuilt.PDC_error_network_calm.(band)] = compute_network_variance(Adj_mat,source_metrics_ICA_A_rebuilt.PDC_network_calm.(band));
[source_metrics_ICA_A_rebuilt.PDC_missclassifications_act.(band) source_metrics_ICA_A_rebuilt.PDC_confmat_act.(band) source_metrics_ICA_A_rebuilt.PDC_confmat_act_metrics.(band) source_metrics_ICA_A_rebuilt.PDC_error_network_act.(band)] = compute_network_variance(Adj_mat,source_metrics_ICA_A_rebuilt.PDC_network_act.(band));

% already did that
% [channel_metrics_ICA_A_rebuilt.coh_missclassifications.(band) channel_metrics_wICA.coh_confmat.(band) channel_metrics_ICA_A_rebuilt.coh_confmat_metrics.(band) channel_metrics_ICA_A_rebuilt.coh_error_network.(band)] = compute_network_variance(channel_metrics_clean.coh_network.(band),channel_metrics_ICA_A_rebuilt.coh_network.(band));
% [channel_metrics_ICA_A_rebuilt.coh_missclassifications_calm.(band) channel_metrics_ICA_A_rebuilt.coh_confmat_calm.(band) channel_metrics_ICA_A_rebuilt.coh_confmat_calm_metrics.(band) channel_metrics_ICA_A_rebuilt.coh_error_network_calm.(band)] = compute_network_variance(channel_metrics_clean.coh_network_calm.(band),channel_metrics_ICA_A_rebuilt.coh_network_calm.(band));
% [channel_metrics_ICA_A_rebuilt.coh_missclassifications_act.(band) channel_metrics_ICA_A_rebuilt.coh_confmat_act.(band) channel_metrics_ICA_A_rebuilt.coh_confmat_act_metrics.(band) channel_metrics_ICA_A_rebuilt.coh_error_network_act.(band)] = compute_network_variance(channel_metrics_clean.coh_network_act.(band),channel_metrics_ICA_A_rebuilt.coh_network_act.(band));
end

channel_metrics_ICA_A_rebuilt.localization_error = local_error;


%% plot source topography
Fs = 256;
%leadfield = -leadfield;
sources_to_keep = 1:7;
K = size(sources_to_keep,2);

%dipole locations manually chosen
regions = headmodel.atlas.colorTable;
active_regions = [1 2 5 6 7 9 33 37 55 58 60 61 64];
active_regions = active_regions(1:size(sources_to_keep,2));
figure()

for i=1:K
    figno = 100+K*10+i;

    tmpidx = find(regions == active_regions(i));
    dipoles_in_region(i) = length(tmpidx);
    ROI_projection = leadfield(:,1,tmpidx(1:5));
    ROI_projection = squeeze(ROI_projection);
    for j=1:5
        ROI_projection(:,j) = ROI_projection(:,j)./max(abs(ROI_projection(:,j)));
    end
    ROI_pro_mean = mean(ROI_projection,2);
    subplot(figno)
    topoplotIndie(ROI_pro_mean, EEG.chanlocs,'numcontour',0,'electrodes','numbers','shading','interp');
    set(gca,'clim',[min(ROI_pro_mean) max(ROI_pro_mean)])
    colorbar('southoutside')
    title(['ROI projection ' num2str(i,'%02d')])
end
%%
metric_to_plot = 'CGCIM_network';
metric_to_plot2 = 'CGCIM_error_network';
source_metrics = source_metrics_ICA;
network = source_metrics.(metric_to_plot);
error_network = source_metrics.(metric_to_plot2);
figure
h = plotnetworktitle(network);
figure
h2 = plotnetworktitle(error_network);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%END OF CODE%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Rebuild eyeblink A - DONT RUN THIS
%might be promising
% didn't work very well because a)find activity doesn't work very well and
% b) ICA eyeblink component isn't very accurate or consistent
activity = find_activity(EEG_noisy(reference_idx,:),0.3);
figure
plot(activity)
hold on
plot(sim_eyeblink)
figure
plot(EEG_noisy(reference_idx,:))
frontal_electrodes = [1 33 34 2 3 37 36 35];
%frontal_electrodes = [1 33 34 2 3 37 36 35 3 4 5 6 7 38 39 40 41 42 8 9 10 11 43 44 45 46 47 12 13 14 15 48 49 50 51 52];;


EEG_affected = EEG_noisy(:,activity);
figure
plotmts2(EEG_affected(1:10,:)')
[icaEEG_affected, A_affected, W_affected] = fastica(EEG_affected,'stabilization','on','verbose','on'); 
[eyeblink_idx_affected best_match_affected] = find_eyeblink_comp(sim_eyeblink(activity),icaEEG_affected);
noiseindex_affected = best_match_affected %change this manually to the noise component index
figure
plotmts2(icaEEG_affected(1:11,:)')
A_correction = A_affected(:,noiseindex_affected);


A_correction = A_correction/max(abs(A_correction));
A_reference = A(33,noiseindex); %%any frontal electrode will do - assuming ICA removes a decent % of noise at that electrode
if A_correction(33)*A_reference < 0 % ??
    A_correction = - A_correction;
end
A_rebuilt = A;
A_rebuilt(:,noiseindex) = A_reference*A_correction;
A_rebuilt(frontal_electrodes,noiseindex) = A(frontal_electrodes,noiseindex); % don't make any changes to frontal electrodes
%A_rebuilt(:,noiseindex) = (A_rebuilt(:,noiseindex) + A(:,noiseindex))/2;
%%
eyeblink_noise = A_rebuilt(:,noiseindex)*icaEEG(noiseindex,:); %calculate clean channel data
EEG_ICA_affected = EEG_noisy - eyeblink_noise;
real_noise = projection_matrix'*sim_eyeblink/amp;

tic
EEG.data = EEG_ICA_affected;
asec = pop_rsbl(EEG, saveFull, account4artifacts, src2roiReductionType, solverType, updateFreq);


ROI_acts = asec.etc.src.act;

ROI_acts = ROI_acts(active_regions,:);
ROI_acts = ROI_acts/dipoles_per_region;

source_metrics_ICA_affected = compute_source_metrics(ROI_acts,xM',Fs,calm_start,calm_end,act_start,act_end);
toc

tic
[source_metrics_ICA_affected.CGCIM_missclassifications source_metrics_ICA_affected.CGCIM_confmat source_metrics_ICA_affected.CGCIM_confmat_metrics source_metrics_ICA_affected.CGCIM_error_network] = compute_network_variance(Adj_mat,source_metrics_ICA_affected.CGCIM_network);
[source_metrics_ICA_affected.CGCIM_missclassifications_calm source_metrics_ICA_affected.CGCIM_confmat_calm source_metrics_ICA_affected.CGCIM_confmat_calm_metrics source_metrics_ICA_affected.CGCIM_error_network_calm] = compute_network_variance(Adj_mat,source_metrics_ICA_affected.CGCIM_network_calm);
[source_metrics_ICA_affected.CGCIM_missclassifications_act source_metrics_ICA_affected.CGCIM_confmat_act source_metrics_ICA_affected.CGCIM_confmat_act_metrics source_metrics_ICA_affected.CGCIM_error_network_act] = compute_network_variance(Adj_mat,source_metrics_ICA_affected.CGCIM_network_act);

channel_metrics_ICA_affected = Compute_channel_metrics(EEG_ICA_affected,EEG_clean,calm_start,calm_end,act_start,act_end);

for i=1:size(frequency_bands,1)
band = frequency_bands(i);
band = string(band);
[source_metrics_ICA_affected.DTF_missclassifications.(band) source_metrics_ICA_affected.DTF_confmat.(band) source_metrics_ICA_affected.DTF_confmat_metrics.(band) source_metrics_ICA_affected.DTF_error_network.(band)] = compute_network_variance(Adj_mat,source_metrics_ICA_affected.DTF_network.(band));
[source_metrics_ICA_affected.DTF_missclassifications_calm.(band) source_metrics_ICA_affected.DTF_confmat_calm.(band) source_metrics_ICA_affected.DTF_confmat_calm_metrics.(band) source_metrics_ICA_affected.DTF_error_network_calm.(band)] = compute_network_variance(Adj_mat,source_metrics_ICA_affected.DTF_network_calm.(band));
[source_metrics_ICA_affected.DTF_missclassifications_act.(band) source_metrics_ICA_affected.DTF_confmat_act.(band) source_metrics_ICA_affected.DTF_confmat_act_metrics.(band) source_metrics_ICA_affected.DTF_error_network_act.(band)] = compute_network_variance(Adj_mat,source_metrics_ICA_affected.DTF_network_act.(band));

[source_metrics_ICA_affected.PDC_missclassifications.(band) source_metrics_ICA_affected.PDC_confmat.(band) source_metrics_ICA_affected.PDC_confmat_metrics.(band) source_metrics_ICA_affected.PDC_error_network.(band)] = compute_network_variance(Adj_mat,source_metrics_ICA_affected.PDC_network.(band));
[source_metrics_ICA_affected.PDC_missclassifications_calm.(band) source_metrics_ICA_affected.PDC_confmat_calm.(band) source_metrics_ICA_affected.PDC_confmat_calm_metrics.(band) source_metrics_ICA_affected.PDC_error_network_calm.(band)] = compute_network_variance(Adj_mat,source_metrics_ICA_affected.PDC_network_calm.(band));
[source_metrics_ICA_affected.PDC_missclassifications_act.(band) source_metrics_ICA_affected.PDC_confmat_act.(band) source_metrics_ICA_affected.PDC_confmat_act_metrics.(band) source_metrics_ICA_affected.PDC_error_network_act.(band)] = compute_network_variance(Adj_mat,source_metrics_ICA_affected.PDC_network_act.(band));

[channel_metrics_ICA_affected.coh_missclassifications.(band) channel_metrics_ICA_affected.coh_confmat.(band) channel_metrics_ICA_affected.coh_confmat_metrics.(band) channel_metrics_ICA_affected.coh_error_network.(band)] = compute_network_variance(channel_metrics_clean.coh_network.(band),channel_metrics_ICA_affected.coh_network.(band));
[channel_metrics_ICA_affected.coh_missclassifications_calm.(band) channel_metrics_ICA_affected.coh_confmat_calm.(band) channel_metrics_ICA_affected.coh_confmat_calm_metrics.(band) channel_metrics_ICA_affected.coh_error_network_calm.(band)] = compute_network_variance(channel_metrics_clean.coh_network_calm.(band),channel_metrics_ICA_affected.coh_network_calm.(band));
[channel_metrics_ICA_affected.coh_missclassifications_act.(band) channel_metrics_ICA_affected.coh_confmat_act.(band) channel_metrics_ICA_affected.coh_confmat_act_metrics.(band) channel_metrics_ICA_affected.coh_error_network_act.(band)] = compute_network_variance(channel_metrics_clean.coh_network_act.(band),channel_metrics_ICA_affected.coh_network_act.(band));
end
toc


ICA_eyeblink_proj = A_rebuilt(:,noiseindex);
ICA_eyeblink_proj = ICA_eyeblink_proj'/max(abs(ICA_eyeblink_proj));
projection_error = true_proj_norm - ICA_eyeblink_proj';
channel_metrics_ICA_affected.projection_error_norm = projection_error;
channel_metrics_ICA_affected.projection_error = projection_error.*projection_matrix';


%% Remove all but one component

%The purpose of this section is to demonstrate that ICA localizes the
%eyeblink component incorrectly. For this, we discard all ICA components
%except the eyeblink, then return back to channel level data. We then
%compare the colormaps between our eyeblink channel level data and the
%projection matrix we originally used to generate the eyeblink noise.

keep = eyeblink_idx; %change this manually to the kept component index
icaEEG4 = zeros(size(icaEEG)); %zero out all ICA components
icaEEG4(keep,:) = icaEEG(keep,:); %restore eyeblink component
%icaEEG4(keep,:) = strong_comp; 

single_comp_data = A*icaEEG4; %calculate single-component channel data


% ICA_eyeblink_proj = A(:,noiseindex);
% ICA_eyeblink_proj = ICA_eyeblink_proj'/max(abs(ICA_eyeblink_proj));
ICA_eyeblink_proj = single_comp_data(:,act_start+40) ;
ICA_eyeblink_proj = ICA_eyeblink_proj/max(abs(ICA_eyeblink_proj));
if ICA_eyeblink_proj(33) * projection_matrix(33) < 0
    ICA_eyeblink_proj = -ICA_eyeblink_proj;
end


true_proj_norm = projection_matrix'/max(abs(projection_matrix));


local_error = true_proj_norm - ICA_eyeblink_proj;
channel_metrics_ICA.localization_error = local_error;
figure()
subplot(131)
topoplotIndie(ICA_eyeblink_proj, EEG.chanlocs,'numcontour',0,'electrodes','numbers','shading','interp');
set(gca,'clim',[0 1])
colorbar('southoutside')
title('ICA eyeblink projection')

subplot(132)
topoplotIndie(true_proj_norm, EEG.chanlocs,'numcontour',0,'electrodes','numbers','shading','interp');
set(gca,'clim',[0 1])

title('real eyeblink component projection')
colorbar('southoutside')

subplot(133)
topoplotIndie(diff, EEG.chanlocs,'numcontour',0,'electrodes','numbers','shading','interp');
set(gca,'clim',[min(diff) max(diff)])
title('projection variance')
colorbar('southoutside')



%% second ICA test

[icaEEG11, A11, W11] = fastica(EEG_ICA,'stabilization','on','verbose','on'); 
[eyeblink_idx2 best_match2] = find_eyeblink_comp(sim_eyeblink,icaEEG11);

figure
plotmts2(EEG_ICA(1:12,:)')
%% Strip by strip ICA
strip.a = [1 33 34 2 3 37 36 35]; %front of the head strip
strip.b = [4:11 38:47];
strip.c = [12:19 32 48:56];
strip.d = [20:26 30 31 57:63];
strip.e = [27 28 29 64]; %back of the head strip
names = ['a' 'b' 'c' 'd' 'e'];
eyeblink_noise = zeros(size(EEG_clean));
for i=1:5
    current_strip = names(i);
    current_strip = strip.(current_strip);
    EEG_strip = EEG_noisy(current_strip,:);
    [icaEEG_strip, A_strip, W_strip] = fastica(EEG_strip,'stabilization','on','verbose','on'); 
    eyeblink_idx_strip = find_eyeblink_comp(sim_eyeblink,icaEEG_strip);
    noiseindex_strip = eyeblink_idx_strip
    eyeblink_noise(current_strip,:) = A_strip(:,noiseindex_strip)*icaEEG_strip(noiseindex_strip,:);
    figure
    plotmts2(icaEEG_strip')
end
    
EEG_rebuilt = EEG_noisy - eyeblink_noise;
figure
plotmts2(EEG_rebuilt(1:15,:)')
figure
plotmts2(EEG_noisy(1:15,:)')

%%


comparison = compare_channel_metrics(channel_metrics_ICA,channel_metrics_clean,act_start,act_end);


ts = comparison.relative_error;
% ts = ts > 0;
figure()
subplot(141)
topoplotIndie(ts, EEG.chanlocs,'numcontour',0,'electrodes','numbers','shading','interp');
set(gca,'clim',[min(ts) max(ts)])
colorbar('southoutside')
title('relative error, whole timeseries')

ts = comparison.relative_error_calm;
%ts = ts > 0;

subplot(142)
topoplotIndie(ts, EEG.chanlocs,'numcontour',0,'electrodes','numbers','shading','interp');
set(gca,'clim',[min(ts) max(ts)])
title('relative error, no eyeblink')
colorbar('southoutside')

ts = comparison.relative_error_act;
%ts = ts > 0;

subplot(143)
topoplotIndie(ts, EEG.chanlocs,'numcontour',0,'electrodes','numbers','shading','interp');
set(gca,'clim',[min(ts) max(ts)])
title('relative error with eyeblink')
colorbar('southoutside')


ts = comparison.error_pct;
%ts = ts.*projection_matrix';
%ts = ts > 0;

subplot(144)
topoplotIndie(ts, EEG.chanlocs,'numcontour',0,'electrodes','numbers','shading','interp');
set(gca,'clim',[0 101])
title('pct error remaining')
colorbar('southoutside')


figure()
h = plotnetworktitle(source_metrics_noisy_ext.RCGCIM_network)
%%
je = EEG_noisy(1,1:4500);

[imf,residual] = emd(je)
for i=1:8
    aa = imf(:,i);
    figure
    plot(aa)
end
%% plot eyeblink topography


figure()
topoplotIndie(true_proj_norm, EEG.chanlocs,'numcontour',0,'electrodes','numbers','shading','interp');
set(gca,'clim',[min(true_proj_norm) 1])
title('Simulated eyeblink projection topography')
colorbar('southoutside')