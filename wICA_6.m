clear all
%load emptyEEG %contains channel location info
load dsi_empty_EEG
headmodel = load('dsi_hm.mat');
projection = load('eyeblink_projection.mat'); %load the projection matrix for the eyeblink component
proj_arrangement = projection.channels;
projection_matrix = projection.mean_ch_dist;
true_proj_norm = projection_matrix'/max(abs(projection_matrix));
sim_eyeblink = load('sim_eyeblink_2.mat').sim_eyeblink; %no WN
sim_eyeblink = downsample(sim_eyeblink,2); % downsample to reduce to 5000 samples
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
EEG.etc = kee.etc;
EEG.data = rand(64,5000);
%%
for seed = 6:13
Fs = 256;
%seed = 3;
rng(seed+50);

%randomize projection matrix 
projection_matrix = projection.mean_ch_dist;
rand_coef = 1 + (rand(1,64) - 0.5)/20;
projection_matrix = projection_matrix.*rand_coef;
true_proj_norm = projection_matrix'/max(abs(projection_matrix));

sources_to_keep = 1:7; %keep 7 sources
K = size(sources_to_keep,2);
n = 5000;
frequency_bands = {'D'; 'Th'; 'A'; 'B'; 'low'};
%sim_eyeblink = load('artif_blink.mat').whole; %load simulated eyeblink component (WN + spikes)

EEG_times = (0:n-1)/Fs;
chan_num = 64;
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



%%Dipole projection to regions


for i=1:K
    tmpidx = find(regions == active_regions(i));
    diploc(i,:) = tmpidx(floor(size(tmpidx,1)/2):floor(size(tmpidx,1)/2)+dipoles_per_region-1);
end

for i=1:68-K
    tmpidx = find(regions == inactive_regions(i));
    diploc_pinknoise(i,:) = tmpidx(randsample(size(tmpidx,1),5));
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


my_dipoles_gain = leadfield(:,1,diploc);
my_dipoles_data = dipole_data([diploc],:);%store the timeseries of active dipoles

EEG.data = EEG_clean;

%%clean data source reconstruction

saveFull = false;
account4artifacts = false; 
src2roiReductionType = 'sum';
solverType = 'bsbl'
updateFreq = 5;
asec = pop_rsbl(EEG, saveFull, account4artifacts, src2roiReductionType, solverType, updateFreq);

%all 68 source
ROI_acts = asec.etc.src.act;
ROI_acts = ROI_acts/dipoles_per_region; %get mean activity per region


channel_metrics_clean = Compute_channel_metrics(EEG_clean,EEG_clean,calm_start,calm_end,act_start,act_end);

[channel_metrics_clean.RCGCIM_missclassifications channel_metrics_clean.RCGCIM_confmat channel_metrics_clean.RCGCIM_confmat_metrics] = compute_network_variance(channel_metrics_clean.RCGCIM_network ,channel_metrics_clean.RCGCIM_network);
[channel_metrics_clean.RCGCIM_missclassifications_calm channel_metrics_clean.RCGCIM_confmat_calm channel_metrics_clean.RCGCIM_confmat_calm_metrics] = compute_network_variance(channel_metrics_clean.RCGCIM_network_calm,channel_metrics_clean.RCGCIM_network_calm);
[channel_metrics_clean.RCGCIM_missclassifications_act channel_metrics_clean.RCGCIM_confmat_act channel_metrics_clean.RCGCIM_confmat_act_metrics] = compute_network_variance(channel_metrics_clean.RCGCIM_network_act,channel_metrics_clean.RCGCIM_network_act);

source_metrics_clean_ext = compute_source_metrics_ext(ROI_acts,Fs,calm_start,calm_end,act_start,act_end);

[source_metrics_clean_ext.RCGCIM_missclassifications source_metrics_clean_ext.RCGCIM_confmat source_metrics_clean_ext.RCGCIM_confmat_metrics] = compute_network_variance(Adj_mat_ext ,source_metrics_clean_ext.RCGCIM_network);
[source_metrics_clean_ext.RCGCIM_missclassifications_calm source_metrics_clean_ext.RCGCIM_confmat_calm source_metrics_clean_ext.RCGCIM_confmat_calm_metrics] = compute_network_variance(Adj_mat_ext,source_metrics_clean_ext.RCGCIM_network_calm);
[source_metrics_clean_ext.RCGCIM_missclassifications_act source_metrics_clean_ext.RCGCIM_confmat_act source_metrics_clean_ext.RCGCIM_confmat_act_metrics] = compute_network_variance(Adj_mat_ext,source_metrics_clean_ext.RCGCIM_network_act);


%7 sources only
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

end




%%noisy data source reconstruction

EEG.data = EEG_noisy;

asec = pop_rsbl(EEG, saveFull, account4artifacts, src2roiReductionType, solverType, updateFreq);

%all 68 source
ROI_acts = asec.etc.src.act;
ROI_acts = ROI_acts/dipoles_per_region;

channel_metrics_noisy = Compute_channel_metrics(EEG_noisy,EEG_clean,calm_start,calm_end,act_start,act_end);

[channel_metrics_noisy.RCGCIM_missclassifications channel_metrics_noisy.RCGCIM_confmat channel_metrics_noisy.RCGCIM_confmat_metrics channel_metrics_noisy.RCGCIM_error_network] = compute_network_variance(channel_metrics_clean.RCGCIM_network ,channel_metrics_noisy.RCGCIM_network);
[channel_metrics_noisy.RCGCIM_missclassifications_calm channel_metrics_noisy.RCGCIM_confmat_calm channel_metrics_noisy.RCGCIM_confmat_calm_metrics channel_metrics_noisy.RCGCIM_error_network_calm] = compute_network_variance(channel_metrics_clean.RCGCIM_network_calm,channel_metrics_noisy.RCGCIM_network_calm);
[channel_metrics_noisy.RCGCIM_missclassifications_act channel_metrics_noisy.RCGCIM_confmat_act channel_metrics_noisy.RCGCIM_confmat_act_metrics channel_metrics_noisy.RCGCIM_error_network_act] = compute_network_variance(channel_metrics_clean.RCGCIM_network_act,channel_metrics_noisy.RCGCIM_network_act);

source_metrics_noisy_ext = compute_source_metrics_ext(ROI_acts,Fs,calm_start,calm_end,act_start,act_end);


[source_metrics_noisy_ext.RCGCIM_missclassifications source_metrics_noisy_ext.RCGCIM_confmat source_metrics_noisy_ext.RCGCIM_confmat_metrics source_metrics_noisy_ext.RCGCIM_error_network] = compute_network_variance(Adj_mat_ext ,source_metrics_noisy_ext.RCGCIM_network);
[source_metrics_noisy_ext.RCGCIM_missclassifications_calm source_metrics_noisy_ext.RCGCIM_confmat_calm source_metrics_noisy_ext.RCGCIM_confmat_calm_metrics source_metrics_noisy_ext.RCGCIM_error_network_calm] = compute_network_variance(Adj_mat_ext,source_metrics_noisy_ext.RCGCIM_network_calm);
[source_metrics_noisy_ext.RCGCIM_missclassifications_act source_metrics_noisy_ext.RCGCIM_confmat_act source_metrics_noisy_ext.RCGCIM_confmat_act_metrics source_metrics_noisy_ext.RCGCIM_error_network_act] = compute_network_variance(Adj_mat_ext,source_metrics_noisy_ext.RCGCIM_network_act);

%compare vs clean data metrics
[source_metrics_noisy_ext.RCGCIM_missclassifications_vs_clean source_metrics_noisy_ext.RCGCIM_confmat_vs_clean source_metrics_noisy_ext.RCGCIM_confmat_metrics_vs_clean] = compute_network_variance(source_metrics_clean_ext.RCGCIM_network ,source_metrics_noisy_ext.RCGCIM_network);
[source_metrics_noisy_ext.RCGCIM_missclassifications_calm_vs_clean source_metrics_noisy_ext.RCGCIM_confmat_calm_vs_clean source_metrics_noisy_ext.RCGCIM_confmat_calm_metrics_vs_clean] = compute_network_variance(source_metrics_clean_ext.RCGCIM_network_calm,source_metrics_noisy_ext.RCGCIM_network_calm);
[source_metrics_noisy_ext.RCGCIM_missclassifications_act_vs_clean source_metrics_noisy_ext.RCGCIM_confmat_act_vs_clean source_metrics_noisy_ext.RCGCIM_confmat_act_metrics_vs_clean] = compute_network_variance(source_metrics_clean_ext.RCGCIM_network_act,source_metrics_noisy_ext.RCGCIM_network_act);



%7 sources

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

end




%%Find independent components (include ICA folder)
% Here for ICA I use fastICA algorithm (package included, for details see 
% the corresponding functions). Note: in the original paper runica from EEGLAB
% was used. You can also test other ICA algorithms at this step.
%
% Note, the use of long (in time) data sets REDUCES the quality of artifact
% suppression (for details see the abovementioned paper).
% Split long files into segments and clean them separately.

[icaEEG, A, W] = fastica(EEG_noisy,'stabilization','on','verbose','on'); 
[eyeblink_idx best_match] = find_eyeblink_comp(sim_eyeblink,icaEEG);
noiseindex = eyeblink_idx; %change this manually to the noise component index
icaEEG3 = icaEEG;
icaEEG3(noiseindex,:) = 0; %remove eyeblink component
EEG_ICA = A*icaEEG3; %calculate clean channel data
% figure
% plotmts2(icaEEG(1:10,:)');

          %ICA data source reconstruction
EEG.data = EEG_ICA;

asec = pop_rsbl(EEG, saveFull, account4artifacts, src2roiReductionType, solverType, updateFreq);

%all 68 source
ROI_acts = asec.etc.src.act;

ROI_acts = ROI_acts/dipoles_per_region;

channel_metrics_ICA = Compute_channel_metrics(EEG_ICA,EEG_clean,calm_start,calm_end,act_start,act_end);

[channel_metrics_ICA.RCGCIM_missclassifications channel_metrics_ICA.RCGCIM_confmat channel_metrics_ICA.RCGCIM_confmat_metrics channel_metrics_ICA.RCGCIM_error_network] = compute_network_variance(channel_metrics_clean.RCGCIM_network ,channel_metrics_ICA.RCGCIM_network);
[channel_metrics_ICA.RCGCIM_missclassifications_calm channel_metrics_ICA.RCGCIM_confmat_calm channel_metrics_ICA.RCGCIM_confmat_calm_metrics channel_metrics_ICA.RCGCIM_error_network_calm] = compute_network_variance(channel_metrics_clean.RCGCIM_network_calm,channel_metrics_ICA.RCGCIM_network_calm);
[channel_metrics_ICA.RCGCIM_missclassifications_act channel_metrics_ICA.RCGCIM_confmat_act channel_metrics_ICA.RCGCIM_confmat_act_metrics channel_metrics_ICA.RCGCIM_error_network_act] = compute_network_variance(channel_metrics_clean.RCGCIM_network_act,channel_metrics_ICA.RCGCIM_network_act);


source_metrics_ICA_ext = compute_source_metrics_ext(ROI_acts,Fs,calm_start,calm_end,act_start,act_end);

[source_metrics_ICA_ext.RCGCIM_missclassifications source_metrics_ICA_ext.RCGCIM_confmat source_metrics_ICA_ext.RCGCIM_confmat_metrics source_metrics_ICA_ext.RCGCIM_error_network] = compute_network_variance(Adj_mat_ext ,source_metrics_ICA_ext.RCGCIM_network);
[source_metrics_ICA_ext.RCGCIM_missclassifications_calm source_metrics_ICA_ext.RCGCIM_confmat_calm source_metrics_ICA_ext.RCGCIM_confmat_calm_metrics source_metrics_ICA_ext.RCGCIM_error_network_calm] = compute_network_variance(Adj_mat_ext,source_metrics_ICA_ext.RCGCIM_network_calm);
[source_metrics_ICA_ext.RCGCIM_missclassifications_act source_metrics_ICA_ext.RCGCIM_confmat_act source_metrics_ICA_ext.RCGCIM_confmat_act_metrics source_metrics_ICA_ext.RCGCIM_error_network_act] = compute_network_variance(Adj_mat_ext,source_metrics_ICA_ext.RCGCIM_network_act);


% 7 sources
ROI_acts = asec.etc.src.act;
ROI_acts = ROI_acts(active_regions,:);
ROI_acts = ROI_acts/dipoles_per_region;

source_metrics_ICA = compute_source_metrics(ROI_acts,xM',Fs,calm_start,calm_end,act_start,act_end);


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
end


keep = noiseindex; 
icaEEG4 = zeros(size(icaEEG)); %zero out all ICA components
icaEEG4(keep,:) = icaEEG(keep,:); %restore eyeblink component

single_comp_data = A*icaEEG4; %calculate single-component channel data


% ICA_eyeblink_proj = A(:,noiseindex);
% ICA_eyeblink_proj = ICA_eyeblink_proj'/max(abs(ICA_eyeblink_proj));
ICA_eyeblink_proj = std(single_comp_data(:,act_start:act_end)') ;
ICA_eyeblink_proj = ICA_eyeblink_proj/max(abs(ICA_eyeblink_proj));
if ICA_eyeblink_proj(33) * projection_matrix(33) < 0
    ICA_eyeblink_proj = -ICA_eyeblink_proj;
end


true_proj_norm = projection_matrix'/max(abs(projection_matrix));


local_error = true_proj_norm - ICA_eyeblink_proj;
channel_metrics_ICA.localization_error = local_error;

%%wICA artifact rejection (include wICA folder)
% NOTE: For better artifact suppression, provide manually the numbers 
% of components to be processed. You can also tune the other arguments
% for your data set.
%

nICs = noiseindex; % Components to be processed, e.g. [1, 4:7]
Kthr = 1.15;             % Tolerance for cleaning artifacts, try: 1, 1.15,...
ArtefThreshold = 3;      % Threshold for detection of ICs with artefacts
                         % Set lower values if you manually select ICs with 
                         % artifacts. Otherwise increase
verbose = 'on';          % print some intermediate results                         
icaEEG2 = RemoveStrongArtifacts(icaEEG, nICs, Kthr, ArtefThreshold, Fs, verbose);

strong_comp = icaEEG(nICs,:) - icaEEG2(nICs,:);
figure('Name','wICA noise component')
plot(strong_comp')

EEG_wICA = A*icaEEG2; %cleaned channel data
   %%wICA data reconstruction

EEG.data = EEG_wICA;
asec = pop_rsbl(EEG, saveFull, account4artifacts, src2roiReductionType, solverType, updateFreq);

%all 68 source
ROI_acts = asec.etc.src.act;

ROI_acts = ROI_acts/dipoles_per_region;

channel_metrics_wICA = Compute_channel_metrics(EEG_wICA,EEG_clean,calm_start,calm_end,act_start,act_end);

[channel_metrics_wICA.RCGCIM_missclassifications channel_metrics_wICA.RCGCIM_confmat channel_metrics_wICA.RCGCIM_confmat_metrics channel_metrics_wICA.RCGCIM_error_network] = compute_network_variance(channel_metrics_clean.RCGCIM_network ,channel_metrics_wICA.RCGCIM_network);
[channel_metrics_wICA.RCGCIM_missclassifications_calm channel_metrics_wICA.RCGCIM_confmat_calm channel_metrics_wICA.RCGCIM_confmat_calm_metrics channel_metrics_wICA.RCGCIM_error_network_calm] = compute_network_variance(channel_metrics_clean.RCGCIM_network_calm,channel_metrics_wICA.RCGCIM_network_calm);
[channel_metrics_wICA.RCGCIM_missclassifications_act channel_metrics_wICA.RCGCIM_confmat_act channel_metrics_wICA.RCGCIM_confmat_act_metrics channel_metrics_wICA.RCGCIM_error_network_act] = compute_network_variance(channel_metrics_clean.RCGCIM_network_act,channel_metrics_wICA.RCGCIM_network_act);


source_metrics_wICA_ext = compute_source_metrics_ext(ROI_acts,Fs,calm_start,calm_end,act_start,act_end);

[source_metrics_wICA_ext.RCGCIM_missclassifications source_metrics_wICA_ext.RCGCIM_confmat source_metrics_wICA_ext.RCGCIM_confmat_metrics source_metrics_wICA_ext.RCGCIM_error_network] = compute_network_variance(Adj_mat_ext ,source_metrics_wICA_ext.RCGCIM_network);
[source_metrics_wICA_ext.RCGCIM_missclassifications_calm source_metrics_wICA_ext.RCGCIM_confmat_calm source_metrics_wICA_ext.RCGCIM_confmat_calm_metrics source_metrics_wICA_ext.RCGCIM_error_network_calm] = compute_network_variance(Adj_mat_ext,source_metrics_wICA_ext.RCGCIM_network_calm);
[source_metrics_wICA_ext.RCGCIM_missclassifications_act source_metrics_wICA_ext.RCGCIM_confmat_act source_metrics_wICA_ext.RCGCIM_confmat_act_metrics source_metrics_wICA_ext.RCGCIM_error_network_act] = compute_network_variance(Adj_mat_ext,source_metrics_wICA_ext.RCGCIM_network_act);


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

end

%%Strip by strip ICA

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

end

EEG_strips = EEG_noisy - eyeblink_noise_strips;
%calculate locality from a single sample during eyeblink
strips_local = eyeblink_noise_strips(:,act_start+40)/max(abs(eyeblink_noise_strips(:,act_start+40)));
local_error = true_proj_norm - strips_local';

%%ICA strips source reconstruction

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

[source_metrics_ICA_strips_ext.RCGCIM_missclassifications source_metrics_ICA_strips_ext.RCGCIM_confmat source_metrics_ICA_strips_ext.RCGCIM_confmat_metrics source_metrics_ICA_strips_ext.RCGCIM_error_network] = compute_network_variance(Adj_mat_ext ,source_metrics_ICA_strips_ext.RCGCIM_network);
[source_metrics_ICA_strips_ext.RCGCIM_missclassifications_calm source_metrics_ICA_strips_ext.RCGCIM_confmat_calm source_metrics_ICA_strips_ext.RCGCIM_confmat_calm_metrics source_metrics_ICA_strips_ext.RCGCIM_error_network_calm] = compute_network_variance(Adj_mat_ext,source_metrics_ICA_strips_ext.RCGCIM_network_calm);
[source_metrics_ICA_strips_ext.RCGCIM_missclassifications_act source_metrics_ICA_strips_ext.RCGCIM_confmat_act source_metrics_ICA_strips_ext.RCGCIM_confmat_act_metrics source_metrics_ICA_strips_ext.RCGCIM_error_network_act] = compute_network_variance(Adj_mat_ext,source_metrics_ICA_strips_ext.RCGCIM_network_act);


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

end

channel_metrics_ICA_strips.localization_error = local_error;

%%Approximate un - mixing matrix

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

%%A rebuilt source reconstruction

EEG.data = EEG_A_rebuilt;
asec = pop_rsbl(EEG, saveFull, account4artifacts, src2roiReductionType, solverType, updateFreq);

%all 68 source
ROI_acts = asec.etc.src.act;

ROI_acts = ROI_acts/dipoles_per_region;

channel_metrics_ICA_A_rebuilt = Compute_channel_metrics(EEG_A_rebuilt,EEG_clean,calm_start,calm_end,act_start,act_end);

[channel_metrics_ICA_A_rebuilt.RCGCIM_missclassifications channel_metrics_ICA_A_rebuilt.RCGCIM_confmat channel_metrics_ICA_A_rebuilt.RCGCIM_confmat_metrics channel_metrics_ICA_A_rebuilt.RCGCIM_error_network] = compute_network_variance(channel_metrics_clean.RCGCIM_network ,channel_metrics_ICA_A_rebuilt.RCGCIM_network);
[channel_metrics_ICA_A_rebuilt.RCGCIM_missclassifications_calm channel_metrics_ICA_A_rebuilt.RCGCIM_confmat_calm channel_metrics_ICA_A_rebuilt.RCGCIM_confmat_calm_metrics channel_metrics_ICA_A_rebuilt.RCGCIM_error_network_calm] = compute_network_variance(channel_metrics_clean.RCGCIM_network_calm,channel_metrics_ICA_A_rebuilt.RCGCIM_network_calm);
[channel_metrics_ICA_A_rebuilt.RCGCIM_missclassifications_act channel_metrics_ICA_A_rebuilt.RCGCIM_confmat_act channel_metrics_ICA_A_rebuilt.RCGCIM_confmat_act_metrics channel_metrics_ICA_A_rebuilt.RCGCIM_error_network_act] = compute_network_variance(channel_metrics_clean.RCGCIM_network_act,channel_metrics_ICA_A_rebuilt.RCGCIM_network_act);


source_metrics_ICA_A_rebuilt_ext = compute_source_metrics_ext(ROI_acts,Fs,calm_start,calm_end,act_start,act_end);

[source_metrics_ICA_A_rebuilt_ext.RCGCIM_missclassifications source_metrics_ICA_A_rebuilt_ext.RCGCIM_confmat source_metrics_ICA_A_rebuilt_ext.RCGCIM_confmat_metrics source_metrics_ICA_A_rebuilt_ext.RCGCIM_error_network] = compute_network_variance(Adj_mat_ext ,source_metrics_ICA_A_rebuilt_ext.RCGCIM_network);
[source_metrics_ICA_A_rebuilt_ext.RCGCIM_missclassifications_calm source_metrics_ICA_A_rebuilt_ext.RCGCIM_confmat_calm source_metrics_ICA_A_rebuilt_ext.RCGCIM_confmat_calm_metrics source_metrics_ICA_A_rebuilt_ext.RCGCIM_error_network_calm] = compute_network_variance(Adj_mat_ext,source_metrics_ICA_A_rebuilt_ext.RCGCIM_network_calm);
[source_metrics_ICA_A_rebuilt_ext.RCGCIM_missclassifications_act source_metrics_ICA_A_rebuilt_ext.RCGCIM_confmat_act source_metrics_ICA_A_rebuilt_ext.RCGCIM_confmat_act_metrics source_metrics_ICA_A_rebuilt_ext.RCGCIM_error_network_act] = compute_network_variance(Adj_mat_ext,source_metrics_ICA_A_rebuilt_ext.RCGCIM_network_act);


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

end

channel_metrics_ICA_A_rebuilt.localization_error = local_error;

%%savefile
filename = ['Pink_data_analysis_' num2str(seed,'%02d') '.mat']
save(['Pink_data_analysis_' num2str(seed,'%02d') '.mat'],'channel_metrics_clean','channel_metrics_noisy','channel_metrics_ICA','channel_metrics_wICA','channel_metrics_ICA_strips','channel_metrics_ICA_A_rebuilt','source_metrics_clean','source_metrics_clean_ext','source_metrics_noisy','source_metrics_noisy_ext','source_metrics_ICA','source_metrics_ICA_ext','source_metrics_wICA','source_metrics_wICA_ext','source_metrics_ICA_strips','source_metrics_ICA_strips_ext','source_metrics_ICA_A_rebuilt','source_metrics_ICA_A_rebuilt_ext','projection_matrix');
end
%% DON'T RUN THESE!
A1 = load('Pink_data_analysis_16');
A2 = load('Pink_data_analysis_15');

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
source_metrics = A1.source_metrics_ICA_strips;
network = source_metrics.(metric_to_plot);
error_network = source_metrics.(metric_to_plot2);
figure
h = plotnetworktitle(network);
figure
h2 = plotnetworktitle(error_network);
%% Compare one value across all experiments
clear metric_to_count
level = 'channel';
cleaning = 'ICA_strips';
metric = 'channel_mean_sq_diff';
metrics_from = [level '_metrics_' cleaning];
for i = 17:31
    filename = ['Pink_data_analysis_' num2str(i,'%02d')]
    A3 = load(filename);
    metric_to_count(i,:) = A3.(metrics_from).(metric);
end
figure
boxplot(metric_to_count)
%% Compare one value across all cleaning methods and all experiments

clear metric_to_count
level = 'source'
cleaning = {'clean'; 'noisy'; 'ICA'; 'wICA'; 'ICA_strips'; 'ICA_A_rebuilt'};
metric = 'CGCIM_missclassifications_act';
metric2 = 'sensitivity'
for i = 17:31
    filename = ['Pink_data_analysis_' num2str(i,'%02d')]
    A3 = load(filename);
    for j = 1:length(cleaning)
        metrics_from = [level '_metrics_' cleaning(j)];
        metrics_from = string([metrics_from{1} metrics_from{2} metrics_from{3}]);
        metric_to_count(i,j,:,:) = A3.(metrics_from).(metric);
    end
end
metric_to_count2 = squeeze(mean(metric_to_count,1)); 
%% compare second level metrics
clear metric_to_count
level = 'source'
cleaning = { 'noisy'; 'ICA'; 'wICA'; 'ICA_strips'; 'ICA_A_rebuilt'};
metric = 'CGCIM_confmat_metrics';
metric2 = 'sensitivity'
for i = 17:31
    filename = ['Pink_data_analysis_' num2str(i,'%02d')]
    A3 = load(filename);
    for j = 1:length(cleaning)
        metrics_from = [level '_metrics_' cleaning(j)];
        metrics_from = string([metrics_from{1} metrics_from{2} metrics_from{3}]);
        metric_to_count(i,j,:,:) = A3.(metrics_from).(metric).(metric2);
    end
end
%% fix localization error
A1 = load('Pink_data_analysis_16');
act_start = 3520;
act_end = 3720;

true_proj = A1.projection_matrix;
true_proj = true_proj/ max(abs(true_proj));

noisy_data = A1.channel_metrics_noisy.data;
ICA_data = A1.channel_metrics_ICA.data;
diff = noisy_data - ICA_data;

ICA_proj = mean_dev(diff(:,act_start:act_end)');
ICA_proj = ICA_proj/max(abs(ICA_proj));

proj_diff = ICA_proj - proj;
%% compute causal metric distance
clear all
[~,Adj_mat] = VAR1RingStructure(5000,7); %generate simulated data
Adj_mat = Adj_mat > 0.01;
Adj_mat = Adj_mat';
A1 = load('Pink_data_analysis_06');

CGC = A1.source_metrics_noisy.CGCIM_act;
diff = CGC - Adj_mat;
non_con_diff = diff.*not(Adj_mat);
con_diff = diff.*Adj_mat;

all_non_con = nonzeros(non_con_diff);
index=isnan(all_non_con)
all_non_con(index) = [];

all_con = nonzeros(con_diff);
index=isnan(all_con)
all_con(index) = [];

max_non_con = max(abs(all_non_con))
min_con = 1 - max(abs(all_con))
%% compute causal metric distance for 100 experiments
clear all
[~,Adj_mat] = VAR1RingStructure(5000,7); %generate simulated data
Adj_mat = Adj_mat > 0.01;
Adj_mat = Adj_mat';

for i = 1:100
    filename = ['Pink_data_analysis_' num2str(i,'%02d')];
    A1 = load(filename);

    CGC = A1.source_metrics_wICA.CGCIM_act;
    diff = CGC - Adj_mat;
    non_con_diff = diff.*not(Adj_mat);
    con_diff = diff.*Adj_mat;

    all_non_con = nonzeros(non_con_diff);
    index=isnan(all_non_con);
    all_non_con(index) = [];

    all_con = nonzeros(con_diff);
    index=isnan(all_con);
    all_con(index) = [];

    max_non_con(i) = max(abs(all_non_con));
    min_con(i) = 1 - max(abs(all_con));
end

max(max_non_con)
min(min_con)
%%
clear all
active_regions = [1 2 5 6 7 9 33];
[~,Adj_mat] = VAR1RingStructure(5000,7); %generate simulated data
Adj_mat = Adj_mat > 0.01;
Adj_mat = Adj_mat';
Adj_mat_ext = extend_mat(Adj_mat,active_regions,68);
all_non_con_everywhere = [];
all_con_everywhere = [];
for i = 1:100
    filename = ['Pink_data_analysis_' num2str(i,'%02d')];
    A1 = load(filename);

    CGC = A1.source_metrics_clean_ext.RCGCIM;
    diff = Adj_mat_ext - CGC;
    non_con_diff = CGC.*not(Adj_mat_ext);
    con_diff = CGC.*Adj_mat_ext;
    
    
    all_non_con = nonzeros(non_con_diff);
    index=isnan(all_non_con);
    all_non_con(index) = [];
    all_non_con_everywhere = [all_non_con_everywhere; all_non_con];
    all_con = nonzeros(con_diff);
    index=isnan(all_con);
    all_con(index) = [];
    all_con_everywhere = [all_con_everywhere; all_con];
    max_non_con(i) = max(abs(all_non_con));
    min_con(i) = 1 - max(abs(all_con));
end

sorted_all_cons = sort(all_con_everywhere);
sorted_all_non_cons = sort(all_non_con_everywhere);

p = 0:1:100;
percentiles_non_con = prctile(sorted_all_non_cons,p);
percentiles_con = prctile(sorted_all_cons,p);