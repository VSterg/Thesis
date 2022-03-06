clear all
%load emptyEEG %contains channel location info
load dsi_empty_EEG
headmodel = load('dsi_hm.mat');
projection = load('eyeblink_projection.mat'); %load the projection matrix for the eyeblink component
proj_arrangement = projection.channels;
projection_matrix = projection.mean_ch_dist;
true_proj_norm = projection_matrix'/max(abs(projection_matrix));
sim_eyeblink = load('sim_eyeblink_2.mat').sim_eyeblink; %no WN
sim_eyeblink = downsample(sim_eyeblink,2);
%The indexing of electrodes between the projection matrix and our EEGlab struct is the
%same, so there's no need to rearrange them.
leadfield(:,1,:) = headmodel.K;
%leadfield = -1*leadfield;
%load (eyes open) EEG dat
EEG_clean = load('S05_restingPre_EC.mat').dataRest;



%%
kee = EEG;
chanlocs = kee.chanlocs;
eeglab
%% after loading .mat into eeglab, update some values
EEG = kee;
EEG.etc = kee.etc;
EEG.data = rand(64,38401);

%%

Fs = 256;
calm_start = 4000/2; %start of calm period (no eyeblinks)
calm_end = 6000/2; %end of calm period
act_start = 7000/2; %start of active period
act_end = 8700/2; %end of active period
%channel names are copied from database documentation
samples_no = size(EEG_clean,2);
while length(sim_eyeblink) < samples_no
    sim_eyeblink = [sim_eyeblink  sim_eyeblink];
end
samples_to_cut = length(sim_eyeblink) - samples_no - 1;
sim_eyeblink(end-samples_to_cut:end) = [];
eyeblink_samples_no = length(sim_eyeblink);

channels = ["Fp1" "AF7" 'AF3' 'F1' 'F3' 'F5' 'F7' 'FT7' 'FC5' 'FC3' 'FC1' 'C1' 'C3' 'C5' 'T7' 'TP7' 'CP5' 'CP3' 'CP1' 'P1' 'P3' 'P5' 'P7' 'P9' 'PO7' 'PO3' 'O1' 'Iz' 'Oz' 'POz' 'Pz' 'CPz' 'Fpz' 'Fp2' 'AF8' 'AF4' 'Afz' 'Fz' 'F2' 'F4' 'F6' 'F8' 'FT8' 'FC6' 'FC4' 'FC2' 'FCz' 'Cz' 'C2' 'C4' 'C6' 'T8' 'TP8' 'CP6' 'CP4' 'CP2' 'P2' 'P4' 'P6' 'P8' 'P10' 'PO8' 'PO4' 'O2'];
EEG_clean(65:end,:) = [];
channels(65:end) = [];
frontal = [find(channels == 'Fp1') find(channels == 'Fp2')];
central = [find(channels == 'C1') find(channels == 'C2')];
occipital = [find(channels == 'O1') find(channels == 'O2')];

figure
plotmts2(EEG_clean(1:15,:)')

amp = 20;
EEG_noisy = EEG_clean + projection_matrix'*sim_eyeblink/amp;
figure
plotmts2(EEG_noisy(1:15,:)')

channel_metrics_clean = Compute_channel_metrics(EEG_clean,EEG_clean,calm_start,calm_end,act_start,act_end);

[channel_metrics_clean.RCGCIM_missclassifications channel_metrics_clean.RCGCIM_confmat channel_metrics_clean.RCGCIM_confmat_metrics] = compute_network_variance(channel_metrics_clean.RCGCIM_network ,channel_metrics_clean.RCGCIM_network);
[channel_metrics_clean.RCGCIM_missclassifications_calm channel_metrics_clean.RCGCIM_confmat_calm channel_metrics_clean.RCGCIM_confmat_calm_metrics] = compute_network_variance(channel_metrics_clean.RCGCIM_network_calm,channel_metrics_clean.RCGCIM_network_calm);
[channel_metrics_clean.RCGCIM_missclassifications_act channel_metrics_clean.RCGCIM_confmat_act channel_metrics_clean.RCGCIM_confmat_act_metrics] = compute_network_variance(channel_metrics_clean.RCGCIM_network_act,channel_metrics_clean.RCGCIM_network_act);


EEG.data = EEG_clean;

%% clean data source reconstruction

saveFull = false;
account4artifacts = false; 
src2roiReductionType = 'mean';
solverType = 'bsbl'
updateFreq = 5;
asec = pop_rsbl(EEG, saveFull, account4artifacts, src2roiReductionType, solverType, updateFreq);
%all 68 source
ROI_acts = asec.etc.src.act;

source_metrics_clean_ext = compute_source_metrics_ext(ROI_acts,Fs,calm_start,calm_end,act_start,act_end);

[source_metrics_clean_ext.RCGCIM_missclassifications source_metrics_clean_ext.RCGCIM_confmat source_metrics_clean_ext.RCGCIM_confmat_metrics] = compute_network_variance(source_metrics_clean_ext.RCGCIM_network ,source_metrics_clean_ext.RCGCIM_network);
[source_metrics_clean_ext.RCGCIM_missclassifications_calm source_metrics_clean_ext.RCGCIM_confmat_calm source_metrics_clean_ext.RCGCIM_confmat_calm_metrics] = compute_network_variance(source_metrics_clean_ext.RCGCIM_network_calm,source_metrics_clean_ext.RCGCIM_network_calm);
[source_metrics_clean_ext.RCGCIM_missclassifications_act source_metrics_clean_ext.RCGCIM_confmat_act source_metrics_clean_ext.RCGCIM_confmat_act_metrics] = compute_network_variance(source_metrics_clean_ext.RCGCIM_network_act,source_metrics_clean_ext.RCGCIM_network_act);




channel_metrics_noisy = Compute_channel_metrics(EEG_noisy,EEG_clean,calm_start,calm_end,act_start,act_end);

[channel_metrics_noisy.RCGCIM_missclassifications channel_metrics_noisy.RCGCIM_confmat channel_metrics_noisy.RCGCIM_confmat_metrics channel_metrics_noisy.RCGCIM_error_network] = compute_network_variance(channel_metrics_clean.RCGCIM_network ,channel_metrics_noisy.RCGCIM_network);
[channel_metrics_noisy.RCGCIM_missclassifications_calm channel_metrics_noisy.RCGCIM_confmat_calm channel_metrics_noisy.RCGCIM_confmat_calm_metrics channel_metrics_noisy.RCGCIM_error_network_calm] = compute_network_variance(channel_metrics_clean.RCGCIM_network_calm,channel_metrics_noisy.RCGCIM_network_calm);
[channel_metrics_noisy.RCGCIM_missclassifications_act channel_metrics_noisy.RCGCIM_confmat_act channel_metrics_noisy.RCGCIM_confmat_act_metrics channel_metrics_noisy.RCGCIM_error_network_act] = compute_network_variance(channel_metrics_clean.RCGCIM_network_act,channel_metrics_noisy.RCGCIM_network_act);


%noisy data source reconstruction
EEG.data = EEG_noisy;

asec = pop_rsbl(EEG, saveFull, account4artifacts, src2roiReductionType, solverType, updateFreq);

%all 68 source
ROI_acts = asec.etc.src.act;

source_metrics_noisy_ext = compute_source_metrics_ext(ROI_acts,Fs,calm_start,calm_end,act_start,act_end);


[source_metrics_noisy_ext.RCGCIM_missclassifications source_metrics_noisy_ext.RCGCIM_confmat source_metrics_noisy_ext.RCGCIM_confmat_metrics source_metrics_noisy_ext.RCGCIM_error_network] = compute_network_variance(source_metrics_clean_ext.RCGCIM_network ,source_metrics_noisy_ext.RCGCIM_network);
[source_metrics_noisy_ext.RCGCIM_missclassifications_calm source_metrics_noisy_ext.RCGCIM_confmat_calm source_metrics_noisy_ext.RCGCIM_confmat_calm_metrics source_metrics_noisy_ext.RCGCIM_error_network_calm] = compute_network_variance(source_metrics_clean_ext.RCGCIM_network_calm,source_metrics_noisy_ext.RCGCIM_network_calm);
[source_metrics_noisy_ext.RCGCIM_missclassifications_act source_metrics_noisy_ext.RCGCIM_confmat_act source_metrics_noisy_ext.RCGCIM_confmat_act_metrics source_metrics_noisy_ext.RCGCIM_error_network_act] = compute_network_variance(source_metrics_clean_ext.RCGCIM_network_act,source_metrics_noisy_ext.RCGCIM_network_act);


%%
K = size(EEG_clean,1);
[icaEEG, A, W] = fastica(EEG_noisy,'stabilization','on','verbose','on'); 
[eyeblink_idx best_match] = find_eyeblink_comp(sim_eyeblink,icaEEG);
noiseindex = eyeblink_idx %change this manually to the noise component index
icaEEG3 = icaEEG;
icaEEG3(noiseindex,:) = mean(icaEEG3(noiseindex,:)); %remove eyeblink component
EEG_ICA = A*icaEEG3; %calculate clean channel data

%plot ICA components to manually find eyeblink (noise) component
figure
plotmts2(icaEEG(7:15,:)')
figure
plotmts2(EEG_clean(1:15,:)')
diff = EEG_clean - EEG_ICA;
figure
plotmts2(diff(1:15,:)')
figure
plotmts2(icaEEG(8:15,:)')

channel_metrics_ICA = Compute_channel_metrics(EEG_ICA,EEG_clean,calm_start,calm_end,act_start,act_end);

[channel_metrics_ICA.RCGCIM_missclassifications channel_metrics_ICA.RCGCIM_confmat channel_metrics_ICA.RCGCIM_confmat_metrics channel_metrics_ICA.RCGCIM_error_network] = compute_network_variance(channel_metrics_clean.RCGCIM_network ,channel_metrics_ICA.RCGCIM_network);
[channel_metrics_ICA.RCGCIM_missclassifications_calm channel_metrics_ICA.RCGCIM_confmat_calm channel_metrics_ICA.RCGCIM_confmat_calm_metrics channel_metrics_ICA.RCGCIM_error_network_calm] = compute_network_variance(channel_metrics_clean.RCGCIM_network_calm,channel_metrics_ICA.RCGCIM_network_calm);
[channel_metrics_ICA.RCGCIM_missclassifications_act channel_metrics_ICA.RCGCIM_confmat_act channel_metrics_ICA.RCGCIM_confmat_act_metrics channel_metrics_ICA.RCGCIM_error_network_act] = compute_network_variance(channel_metrics_clean.RCGCIM_network_act,channel_metrics_ICA.RCGCIM_network_act);

EEG.data = EEG_ICA;

asec = pop_rsbl(EEG, saveFull, account4artifacts, src2roiReductionType, solverType, updateFreq);

%all 68 source
ROI_acts = asec.etc.src.act;

source_metrics_ICA_ext = compute_source_metrics_ext(ROI_acts,Fs,calm_start,calm_end,act_start,act_end);

[source_metrics_ICA_ext.RCGCIM_missclassifications source_metrics_ICA_ext.RCGCIM_confmat source_metrics_ICA_ext.RCGCIM_confmat_metrics source_metrics_ICA_ext.RCGCIM_error_network] = compute_network_variance(source_metrics_clean_ext.RCGCIM_network ,source_metrics_ICA_ext.RCGCIM_network);
[source_metrics_ICA_ext.RCGCIM_missclassifications_calm source_metrics_ICA_ext.RCGCIM_confmat_calm source_metrics_ICA_ext.RCGCIM_confmat_calm_metrics source_metrics_ICA_ext.RCGCIM_error_network_calm] = compute_network_variance(source_metrics_clean_ext.RCGCIM_network_calm,source_metrics_ICA_ext.RCGCIM_network_calm);
[source_metrics_ICA_ext.RCGCIM_missclassifications_act source_metrics_ICA_ext.RCGCIM_confmat_act source_metrics_ICA_ext.RCGCIM_confmat_act_metrics source_metrics_ICA_ext.RCGCIM_error_network_act] = compute_network_variance(source_metrics_clean_ext.RCGCIM_network_act,source_metrics_ICA_ext.RCGCIM_network_act);

%ROI_acts = ROI_acts/dipoles_per_region;

%%
Fs = 256;
nICs = noiseindex; % Components to be processed, e.g. [1, 4:7]
Kthr = 1.15;             % Tolerance for cleaning artifacts, try: 1, 1.15,...
ArtefThreshold = 4;      % Threshold for detection of ICs with artefacts
                         % Set lower values if you manually select ICs with 
                         % artifacts. Otherwise increase
verbose = 'on';          % print some intermediate results                         
icaEEG2 = RemoveStrongArtifacts(icaEEG, nICs, Kthr, ArtefThreshold, Fs, verbose);

strong_comp = icaEEG(nICs,:) - icaEEG2(nICs,:);
figure('Name','wICA noise component')
plot(strong_comp')

EEG_wICA = A*icaEEG2; %cleaned channel data

figure
plot(icaEEG(noiseindex,:))
%%



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
    eyeblink_noise_strips(current_strip,:) = A_strip(:,noiseindex_strip)*(icaEEG_strip(noiseindex_strip,:) - mean(icaEEG_strip(noiseindex_strip,:)));
    figure
    plotmts2(icaEEG_strip(1:end,:)');
end

EEG_strips = EEG_noisy - eyeblink_noise_strips;
figure
plotmts2(EEG_strips(1:15,:)')
channel_metrics_ICA_strips = Compute_channel_metrics(EEG_strips - mean(EEG_strips,2),EEG_clean-mean(EEG_clean,2),calm_start,calm_end,act_start,act_end);

[channel_metrics_ICA_strips.RCGCIM_missclassifications channel_metrics_ICA_strips.RCGCIM_confmat channel_metrics_ICA_strips.RCGCIM_confmat_metrics channel_metrics_ICA_strips.RCGCIM_error_network] = compute_network_variance(channel_metrics_clean.RCGCIM_network ,channel_metrics_ICA_strips.RCGCIM_network);
[channel_metrics_ICA_strips.RCGCIM_missclassifications_calm channel_metrics_ICA_strips.RCGCIM_confmat_calm channel_metrics_ICA_strips.RCGCIM_confmat_calm_metrics channel_metrics_ICA_strips.RCGCIM_error_network_calm] = compute_network_variance(channel_metrics_clean.RCGCIM_network_calm,channel_metrics_ICA_strips.RCGCIM_network_calm);
[channel_metrics_ICA_strips.RCGCIM_missclassifications_act channel_metrics_ICA_strips.RCGCIM_confmat_act channel_metrics_ICA_strips.RCGCIM_confmat_act_metrics channel_metrics_ICA_strips.RCGCIM_error_network_act] = compute_network_variance(channel_metrics_clean.RCGCIM_network_act,channel_metrics_ICA_strips.RCGCIM_network_act);


EEG.data = EEG_strips;
asec = pop_rsbl(EEG, saveFull, account4artifacts, src2roiReductionType, solverType, updateFreq);

%all 68 source
ROI_acts = asec.etc.src.act;

source_metrics_ICA_strips_ext = compute_source_metrics_ext(ROI_acts,Fs,calm_start,calm_end,act_start,act_end);

[source_metrics_ICA_strips_ext.RCGCIM_missclassifications source_metrics_ICA_strips_ext.RCGCIM_confmat source_metrics_ICA_strips_ext.RCGCIM_confmat_metrics source_metrics_ICA_strips_ext.RCGCIM_error_network] = compute_network_variance(source_metrics_clean_ext.RCGCIM_network ,source_metrics_ICA_strips_ext.RCGCIM_network);
[source_metrics_ICA_strips_ext.RCGCIM_missclassifications_calm source_metrics_ICA_strips_ext.RCGCIM_confmat_calm source_metrics_ICA_strips_ext.RCGCIM_confmat_calm_metrics source_metrics_ICA_strips_ext.RCGCIM_error_network_calm] = compute_network_variance(source_metrics_clean_ext.RCGCIM_network_calm,source_metrics_ICA_strips_ext.RCGCIM_network_calm);
[source_metrics_ICA_strips_ext.RCGCIM_missclassifications_act source_metrics_ICA_strips_ext.RCGCIM_confmat_act source_metrics_ICA_strips_ext.RCGCIM_confmat_act_metrics source_metrics_ICA_strips_ext.RCGCIM_error_network_act] = compute_network_variance(source_metrics_clean_ext.RCGCIM_network_act,source_metrics_ICA_strips_ext.RCGCIM_network_act);


%%


reference_idx = 33; % a frontal electrode
frontal_electrodes = [1 33 34 2 3 37 36 35];
eyeblink_noise_reference = EEG_noisy(reference_idx,:) - EEG_ICA(reference_idx,:);%eyeblink_noise_strips(reference_idx,:); %get reference noise from ICA_strips, not ICA

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

figure
plotmts2(EEG_A_rebuilt(1:15,:)')
figure
plotmts2(EEG_clean(1:15,:)')

channel_metrics_ICA_A_rebuilt = Compute_channel_metrics(EEG_A_rebuilt - mean(EEG_A_rebuilt,2),EEG_clean-mean(EEG_clean,2),calm_start,calm_end,act_start,act_end);

[channel_metrics_ICA_A_rebuilt.RCGCIM_missclassifications channel_metrics_ICA_A_rebuilt.RCGCIM_confmat channel_metrics_ICA_A_rebuilt.RCGCIM_confmat_metrics channel_metrics_ICA_A_rebuilt.RCGCIM_error_network] = compute_network_variance(channel_metrics_clean.RCGCIM_network ,channel_metrics_ICA_A_rebuilt.RCGCIM_network);
[channel_metrics_ICA_A_rebuilt.RCGCIM_missclassifications_calm channel_metrics_ICA_A_rebuilt.RCGCIM_confmat_calm channel_metrics_ICA_A_rebuilt.RCGCIM_confmat_calm_metrics channel_metrics_ICA_A_rebuilt.RCGCIM_error_network_calm] = compute_network_variance(channel_metrics_clean.RCGCIM_network_calm,channel_metrics_ICA_A_rebuilt.RCGCIM_network_calm);
[channel_metrics_ICA_A_rebuilt.RCGCIM_missclassifications_act channel_metrics_ICA_A_rebuilt.RCGCIM_confmat_act channel_metrics_ICA_A_rebuilt.RCGCIM_confmat_act_metrics channel_metrics_ICA_A_rebuilt.RCGCIM_error_network_act] = compute_network_variance(channel_metrics_clean.RCGCIM_network_act,channel_metrics_ICA_A_rebuilt.RCGCIM_network_act);

EEG.data = EEG_A_rebuilt;
asec = pop_rsbl(EEG, saveFull, account4artifacts, src2roiReductionType, solverType, updateFreq);

%all 68 source
ROI_acts = asec.etc.src.act;

source_metrics_ICA_A_rebuilt_ext = compute_source_metrics_ext(ROI_acts,Fs,calm_start,calm_end,act_start,act_end);

[source_metrics_ICA_A_rebuilt_ext.RCGCIM_missclassifications source_metrics_ICA_A_rebuilt_ext.RCGCIM_confmat source_metrics_ICA_A_rebuilt_ext.RCGCIM_confmat_metrics source_metrics_ICA_A_rebuilt_ext.RCGCIM_error_network] = compute_network_variance(source_metrics_clean_ext.RCGCIM_network ,source_metrics_ICA_A_rebuilt_ext.RCGCIM_network);
[source_metrics_ICA_A_rebuilt_ext.RCGCIM_missclassifications_calm source_metrics_ICA_A_rebuilt_ext.RCGCIM_confmat_calm source_metrics_ICA_A_rebuilt_ext.RCGCIM_confmat_calm_metrics source_metrics_ICA_A_rebuilt_ext.RCGCIM_error_network_calm] = compute_network_variance(source_metrics_clean_ext.RCGCIM_network_calm,source_metrics_ICA_A_rebuilt_ext.RCGCIM_network_calm);
[source_metrics_ICA_A_rebuilt_ext.RCGCIM_missclassifications_act source_metrics_ICA_A_rebuilt_ext.RCGCIM_confmat_act source_metrics_ICA_A_rebuilt_ext.RCGCIM_confmat_act_metrics source_metrics_ICA_A_rebuilt_ext.RCGCIM_error_network_act] = compute_network_variance(source_metrics_clean_ext.RCGCIM_network_act,source_metrics_ICA_A_rebuilt_ext.RCGCIM_network_act);

%%

plot_these = 46:64
figure
plotmts2(EEG_ICA(plot_these,:)')
figure
plotmts2(EEG_strips(plot_these,:)')
figure
plotmts2(EEG_A_rebuilt(plot_these,:)')
figure
plotmts2(channel_metrics_ICA_A_rebuilt.channel_diff(1:15,:)')
figure
plotmts2(channel_metrics_ICA.channel_diff(1:15,:)')
%%
plot_this = EEG_ICA(1:5,:) - EEG_clean(1:5,:);
figure
plotmts2(plot_this')
plot_this_3 = channel_metrics_ICA.channel_diff(1:5,:);
figure
plotmts2(plot_this_3')
plot_this_2 = EEG_A_rebuilt(41:50,:) - EEG_clean(41:50,:);
figure
plotmts2(plot_this_2')
figure
plotmts2(channel_metrics_ICA_A_rebuilt.channel_diff(1:15,:)')
aa = diff(33,:) + channel_metrics_ICA.channel_diff(33,:);
%%

active_regions = [1 2 5 6 7 9 33];
ROI_mapping = [1:2:9 13:2:17 21:2:31 35 33 37:2:63 11 65 67 19 2:2:10 14:2:18 22:2:32 36 34 38:2:64 12 66 68 20];
Adj_mat_ext = a3;
%Adj_mat_ext = extend_mat(Adj_mat_ext,active_regions,68);
for i=1:68
    for j = 1:68
        Adj_mat_ext2(i,j) = Adj_mat_ext(ROI_mapping(i),ROI_mapping(j));
    end
end
writematrix(Adj_mat_ext2,'tabledata2.txt','Delimiter','\t')
%%
a1 = RCGCIM_network_noisy - RCGCIM_network_ICA;
a2 = RCGCIM_network_noisy - RCGCIM_network_strips;
a3 = RCGCIM_network_noisy - RCGCIM_network_A_rebuilt;
aa1 = abs(a1 - a2);
aa2 = abs(a1 - a3);
aaa1 = sum(sum(aa1));
aaa2 = sum(sum(aa2));
%%