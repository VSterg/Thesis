clear all
%load emptyEEG %contains channel location info
load dsi_empty_EEG
EEG_noisy = load('S05_restingPre_EO.mat').dataRest;
figure
plotmts2(EEG_noisy(1:10,:)')
data_frontal = EEG_noisy(33,:);
EEG_noisy(65:end,:) = [];
calm_start = 6450;
calm_end = 6900;
act_start = 5500;
act_end = 5750;
%%
kee = EEG;
chanlocs = kee.chanlocs;
eeglab
%% after loading .mat into eeglab, update some values
EEG = kee;
EEG.etc = kee.etc;
EEG.data = rand(64,38401);
%% noisy source reconstruction
EEG.data = EEG_noisy
saveFull = false;
account4artifacts = false; 
src2roiReductionType = 'mean';
solverType = 'bsbl'
updateFreq = 5;
asec = pop_rsbl(EEG, saveFull, account4artifacts, src2roiReductionType, solverType, updateFreq);
%all 68 source
ROI_acts = asec.etc.src.act;

beta = 0.25; %threshold for RCGCIM
alpha = 0.05; %threshold for rCGCIM

pmax = 1;
maketest = 1;

[RCGCIM_elec_noisy,pRCGCIM_elec_noisy]=mBTSCGCImatrix(EEG_noisy(:,5000:10000)',pmax,maketest);
[RCGCIM_elec_noisy_calm,pRCGCIM_elec_noisy_calm]=mBTSCGCImatrix(EEG_noisy(:,calm_start:calm_end)',pmax,maketest);
[RCGCIM_elec_noisy_act,pRCGCIM_elec_noisy_act]=mBTSCGCImatrix(EEG_noisy(:,act_start:act_end)',pmax,maketest);
RCGCIM_network_elec_noisy = RCGCIM_elec_noisy > alpha;
RCGCIM_network_elec_noisy_calm = RCGCIM_elec_noisy_calm > alpha;
RCGCIM_network_elec_noisy_act = RCGCIM_elec_noisy_act > alpha;


[RCGCIM_noisy,pRCGCIM_noisy]=mBTSCGCImatrix(ROI_acts(:,5000:10000)',pmax,maketest);
[RCGCIM_noisy_calm,pRCGCIM_noisy_calm]=mBTSCGCImatrix(ROI_acts(:,calm_start:calm_end)',pmax,maketest);
[RCGCIM_noisy_act,pRCGCIM_noisy_act]=mBTSCGCImatrix(ROI_acts(:,act_start:act_end)',pmax,maketest);
RCGCIM_network_noisy = RCGCIM_noisy > beta;
RCGCIM_network_noisy_calm = RCGCIM_noisy_calm > beta;
RCGCIM_network_noisy_act = RCGCIM_noisy_act > beta;
%% ICA
K = size(EEG_noisy,1);
[icaEEG, A, W] = fastica(EEG_noisy,'stabilization','on','verbose','on'); 
[eyeblink_idx best_match] = find_eyeblink_comp(data_frontal,icaEEG);
figure
plotmts2(icaEEG(1:15,:)')
figure
plotmts2(icaEEG(16:30,:)')
figure
plotmts2(icaEEG(31:45,:)')
figure
plotmts2(icaEEG(45:end,:)')
%%
noiseindex = [24];
figure
plotmts2(icaEEG(noiseindex,:)')
icaEEG3 = icaEEG;
for i = 1:length(noiseindex)
    icaEEG3(noiseindex(i),:) = mean(icaEEG(noiseindex(i),:),2);
end%remove eyeblink component
EEG_ICA = A*icaEEG3; %calculate clean channel data
figure
plotmts2(EEG_ICA(1:10,:)')
figure
plotmts2(EEG_noisy(1:10,:)')

%% noisy source reconstruction
EEG.data = EEG_ICA

asec = pop_rsbl(EEG, saveFull, account4artifacts, src2roiReductionType, solverType, updateFreq);
%all 68 source
ROI_acts = asec.etc.src.act;

beta = 0.25; %threshold for RCGCIM
alpha = 0.05; %threshold for rCGCIM

pmax = 1;
maketest = 1;

[RCGCIM_elec_ICA,pRCGCIM_elec_ICA]=mBTSCGCImatrix(EEG_ICA(:,5000:10000)',pmax,maketest);
[RCGCIM_elec_ICA_calm,pRCGCIM_elec_ICA_calm]=mBTSCGCImatrix(EEG_ICA(:,calm_start:calm_end)',pmax,maketest);
[RCGCIM_elec_ICA_act,pRCGCIM_elec_ICA_act]=mBTSCGCImatrix(EEG_ICA(:,act_start:act_end)',pmax,maketest);
RCGCIM_network_elec_ICA = RCGCIM_elec_ICA > alpha;
RCGCIM_network_elec_ICA_calm = RCGCIM_elec_ICA_calm > alpha;
RCGCIM_network_elec_ICA_act = RCGCIM_elec_ICA_act > alpha;


[RCGCIM_ICA,pRCGCIM_ICA]=mBTSCGCImatrix(ROI_acts(:,5000:10000)',pmax,maketest);
[RCGCIM_ICA_calm,pRCGCIM_ICA_calm]=mBTSCGCImatrix(ROI_acts(:,calm_start:calm_end)',pmax,maketest);
[RCGCIM_ICA_act,pRCGCIM_ICA_act]=mBTSCGCImatrix(ROI_acts(:,act_start:act_end)',pmax,maketest);
RCGCIM_network_ICA = RCGCIM_ICA> beta;
RCGCIM_network_ICA_calm = RCGCIM_ICA_calm > beta;
RCGCIM_network_ICA_act = RCGCIM_ICA_act > beta;
%%

strip.a = [1 33 34 2 3 37 36 35 4:11 38:47]; %front of the head strip
strip.b = [4:11 38:47];
strip.c = [12:19 32 48:56];
strip.d = [20:26 30 31 57:63];
strip.e = [27 28 29 64]; %back of the head strip
names = ['a' 'b' 'c' 'd' 'e'];
eyeblink_noise_strips = zeros(size(EEG_noisy));
sim_eyeblink = icaEEG(noiseindex(1),:);
for i=1:1
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
figure
plotmts2(eyeblink_noise_strips(1:10,:)')

EEG_strips = EEG_noisy - eyeblink_noise_strips;
figure
plotmts2(EEG_strips(20:26,:)')
%%
EEG.data = EEG_strips

asec = pop_rsbl(EEG, saveFull, account4artifacts, src2roiReductionType, solverType, updateFreq);
%all 68 source
ROI_acts = asec.etc.src.act;

beta = 0.25; %threshold for RCGCIM
alpha = 0.05; %threshold for rCGCIM

pmax = 1;
maketest = 1;

[RCGCIM_elec_strips,pRCGCIM_elec_strips]=mBTSCGCImatrix(EEG_strips(:,5000:10000)',pmax,maketest);
[RCGCIM_elec_strips_calm,pRCGCIM_elec_strips_calm]=mBTSCGCImatrix(EEG_strips(:,calm_start:calm_end)',pmax,maketest);
[RCGCIM_elec_strips_act,pRCGCIM_elec_strips_act]=mBTSCGCImatrix(EEG_strips(:,act_start:act_end)',pmax,maketest);
RCGCIM_network_elec_strips = RCGCIM_elec_strips > alpha;
RCGCIM_network_elec_strips_calm = RCGCIM_elec_strips_calm > alpha;
RCGCIM_network_elec_strips_act = RCGCIM_elec_strips_act > alpha;


[RCGCIM_strips,pRCGCIM_strips]=mBTSCGCImatrix(ROI_acts(:,5000:10000)',pmax,maketest);
[RCGCIM_strips_calm,pRCGCIM_strips_calm]=mBTSCGCImatrix(ROI_acts(:,calm_start:calm_end)',pmax,maketest);
[RCGCIM_strips_act,pRCGCIM_strips_act]=mBTSCGCImatrix(ROI_acts(:,act_start:act_end)',pmax,maketest);
RCGCIM_network_strips = RCGCIM_strips> beta;
RCGCIM_network_strips_calm = RCGCIM_strips_calm > beta;
RCGCIM_network_strips_act = RCGCIM_strips_act > beta;
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
%local_error = true_proj_norm - mean_ch_dist_norm';

eyeblink_noise_2 = mean_ch_dist_norm'*eyeblink_noise_reference ;

EEG_A_rebuilt = EEG_noisy - eyeblink_noise_2;
range = 46:64
figure
plotmts2(EEG_A_rebuilt(range,:)')
figure
plotmts2(EEG_ICA(range,:)')
%%
EEG.data = EEG_A_rebuilt

asec = pop_rsbl(EEG, saveFull, account4artifacts, src2roiReductionType, solverType, updateFreq);
%all 68 source
ROI_acts = asec.etc.src.act;

beta = 0.25; %threshold for RCGCIM
alpha = 0.05; %threshold for rCGCIM

pmax = 1;
maketest = 1;

[RCGCIM_elec_A_rebuilt,pRCGCIM_elec_A_rebuilt]=mBTSCGCImatrix(EEG_A_rebuilt(:,5000:10000)',pmax,maketest);
[RCGCIM_elec_A_rebuilt_calm,pRCGCIM_elec_A_rebuilt_calm]=mBTSCGCImatrix(EEG_A_rebuilt(:,calm_start:calm_end)',pmax,maketest);
[RCGCIM_elec_A_rebuilt_act,pRCGCIM_elec_A_rebuilt_act]=mBTSCGCImatrix(EEG_A_rebuilt(:,act_start:act_end)',pmax,maketest);
RCGCIM_network_elec_A_rebuilt = RCGCIM_elec_A_rebuilt > alpha;
RCGCIM_network_elec_A_rebuilt_calm = RCGCIM_elec_A_rebuilt_calm > alpha;
RCGCIM_network_elec_A_rebuilt_act = RCGCIM_elec_A_rebuilt_act > alpha;


[RCGCIM_A_rebuilt,pRCGCIM_A_rebuilt]=mBTSCGCImatrix(ROI_acts(:,5000:10000)',pmax,maketest);
[RCGCIM_A_rebuilt_calm,pRCGCIM_A_rebuilt_calm]=mBTSCGCImatrix(ROI_acts(:,calm_start:calm_end)',pmax,maketest);
[RCGCIM_A_rebuilt_act,pRCGCIM_A_rebuilt_act]=mBTSCGCImatrix(ROI_acts(:,act_start:act_end)',pmax,maketest);
RCGCIM_network_A_rebuilt = RCGCIM_A_rebuilt> beta;
RCGCIM_network_A_rebuilt_calm = RCGCIM_A_rebuilt_calm > beta;
RCGCIM_network_A_rebuilt_act = RCGCIM_A_rebuilt_act > beta;
%%

% figure
% plotmts2(EEG_A_rebuilt(1:15,:)')
% figure
% plotmts2(EEG_A_rebuilt(16:30,:)')
% figure
% plotmts2(EEG_A_rebuilt(31:45,:)')
% figure
% plotmts2(EEG_A_rebuilt(46:end,:)')
diff = EEG_A_rebuilt - EEG_ICA;
figure
plotmts2(diff(1:9,:)')
idx =62;
a1 = EEG_ICA(idx,:) - mean(EEG_ICA(idx,:));
a2 = EEG_A_rebuilt(idx,:) - mean(EEG_A_rebuilt(idx,:));
figure
plot(a1)
figure
plot(a2)
%%
load dsi_empty_EEG
topoplotIndie(A(:,noiseindex), EEG.chanlocs,'numcontour',0,'electrodes','numbers','shading','interp');
set(gca,'clim',[min(A(:,noiseindex)) max(A(:,noiseindex))])
colorbar('southoutside')
title(['ROI projection ' num2str(i,'%02d')])
%%
info = RCGCIM_network_elec_A_rebuilt;
figure
G = digraph(info);
p = plot(G,'layout','circle')
p.XData = [EEG.chanlocs(:).X];
p.YData = [EEG.chanlocs(:).Y];
p.ZData = [EEG.chanlocs(:).Z];