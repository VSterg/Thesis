clear all
load emptyEEG %contains channel location info
sim_eyeblink = load('artif_blink.mat').whole; %load simulated eyeblink component (WN + spikes)
xM = load('xM2.mat').xM; %load simulated sources

projection = load('eyeblink_projection.mat'); %load the projection matrix for the eyeblink component
proj_arrangement = projection.channels;
projection_matrix = projection.mean_ch_dist;
%The indexing of electrodes between the projection matrix and our EEGlab struct is the
%same, so there's no need to rearrange them.


%% dipole generation and projecting to channels

K = 12; %number of simulated sources
n = 10000; %number of time samples
Fs = 256; %sampling rate
%dipole locations manually chosen
diploc = [62 1918 744 199 223 255 304 888 925 1033 1034 1684 155 746 544 654 123 765 971 323 542 756 323 45 778 442 675 235 765 32 500 549 284 253 679 184 834 135 742 418 986 836 283 259 813 983 814 935 873 913 986 93 73 875 961 853 273 983 341 312 653 131 764 764 435 874]; 
diploc = diploc(1:K); %keep as many dipoles as we have (neural) sources
%reduce number of samples in signal
EEG.pnts  = n; 
EEG.times = (0:EEG.pnts-1)/EEG.srate;
chan_num = EEG.nbchan;
xM = xM(1:10000,:); 

% initialize all dipole data
dipole_data = zeros(size(lf.Gain,3),EEG.pnts);
for i = 2:K+1 %the simulated neural components are added to manually selected dipoles
    dipole_data(diploc(i-1),:) = xM(:,i-1);
end

EEG_clean(:,:) = squeeze(lf.Gain(:,1,:))*dipole_data; %project dipole data to electrodes (no eyeblink component)
sim_eyeblink = sim_eyeblink(1:n); %n must not be higher than 10000

%add eyeblink component projection to EEG data
EEG.data(:,:) = EEG_clean + projection_matrix'*sim_eyeblink/5000; 

figure('Name','Noisy data, first 10 EEG channels')
plotmts(EEG.data(1:10,:)',1,1,chan_num,1/Fs)
figure('Name','Clean data, first 10 EEG channels')
plotmts(EEG_clean(1:10,:)',1,1,chan_num,1/Fs,[],2)

my_dipoles_gain = lf.Gain(:,:,[diploc]);%store the gains of active dipoles
my_dipoles_data = dipole_data([diploc],:);%store the timeseries of active dipoles

%% Find independent components
% Here for ICA I use fastICA algorithm (package included, for details see 
% the corresponding functions). Note: in the original paper runica from EEGLAB
% was used. You can also test other ICA algorithms at this step.
%
% Note, the use of long (in time) data sets REDUCES the quality of artifact
% suppression (for details see the abovementioned paper).
% Split long files into segments and clean them separately.

[icaEEG, A, W] = fastica(EEG.data,'stabilization','on','verbose','on'); 
data = A*icaEEG;

figure('Name','ICA components timeseries')
plotmts(icaEEG(:,:)',1,0,chan_num,1/Fs,[],1) %find eyeblink component manually


%% Remove all but one component

%The purpose of this section is to demonstrate that ICA localizes the
%eyeblink component incorrectly. For this, we discard all ICA components
%except the eyeblink, then return back to channel level data. We then
%compare the ratio between our eyeblink channel lete data and the
%projection matrix we originally used to generate the eyeblink noise.

keep = 1; %change this manually to the kept component index
icaEEG2 = zeros(size(icaEEG)); %zero out all ICA components
icaEEG2(keep,:) = icaEEG(keep,:); %restore eyeblink component

clean_data = A*icaEEG2; %calculate single-component channel data

projection_ratio = std(clean_data');
projection_ratio = projection_ratio/max(projection_ratio);

true_projection_ratio = projection_matrix/max(projection_matrix);

projection_ratio_diff = projection_ratio - true_projection_ratio;
%%
figure
plot(clean_data