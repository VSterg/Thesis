function metrics = compute_channel_metrics(data,clean_data,calm_start,calm_end,act_start,act_end)
%This function computes a few channel level metrics for EEG data. Inputs
%are: 
%data : the EEG data we measure metrics for
%clean_data : the clean EEG data, the time
%indexes when a calm period starts and ends (calm_start - calm_end) and
%the time indexes when an active (noisy) period starts and ends (act_start - act_end)
%Outputs is metrics : a struct containing the EEG data, channel difference,
%mean square channel difference and standard deviation, power spectrum and power spectrum
%variance in regards to the clean signal and connectivity measures (Coherence, RCGCIM) as well
%as a corresponding network mixing matrix for each connectivity measure.
metrics.data = data;

metrics.channel_diff = data - clean_data;
metrics.channel_sq_diff = metrics.channel_diff.^2;
metrics.channel_sq_diff_calm = metrics.channel_sq_diff(:,calm_start:calm_end);
metrics.channel_sq_diff_act = metrics.channel_sq_diff(:,act_start:act_end);
metrics.channel_mean_sq_diff = mean(metrics.channel_sq_diff,2); %channel mean square difference
metrics.channel_mean_sq_diff_calm = mean(metrics.channel_sq_diff_calm,2);
metrics.channel_mean_sq_diff_act = mean(metrics.channel_sq_diff_act,2);
metrics.channel_diff_calm_dev = std(metrics.channel_diff(:,calm_start:calm_end)');
metrics.channel_dif_act_dev = std(metrics.channel_diff(:,act_start:act_end)');
metrics.channel_diff_overall_dev = std(metrics.channel_diff');
metrics.residual_var_norm = var(metrics.channel_diff')./var(clean_data');
metrics.residual_var_norm_calm = var(metrics.channel_diff(:,calm_start:calm_end)')./var(clean_data(:,calm_start:calm_end)');
metrics.residual_var_norm_act = var(metrics.channel_diff(:,act_start:act_end)')./var(clean_data(:,act_start:act_end)');


% metrics.channel_power = (2*abs(fft(data,[],2))/size(data,2)).^2;
% clean_channel_power = (2*abs(fft(clean_data,[],2))/size(clean_data,2)).^2;
% metrics.channel_power_diff = metrics.channel_power - clean_channel_power;
% 
% %connectivity measures
% 
% x1 = data'; 
% x1_calm = data(:,calm_start:calm_end)';
% x1_act = data(:,act_start:act_end)';
% for i=1:size(data,1)
%     x2 = data(i,:)';
%     x2_calm = data(i,calm_start:calm_end)';
%     x2_act = data(i,act_start:act_end)';
%     metrics.coherency(:,:,i) = mscohere(x1,x2,300,100);
%     metrics.coherency_calm(:,:,i) = mscohere(x1_calm,x2_calm);
%     metrics.coherency_act(:,:,i) = mscohere(x1_act,x2_act);
% end
% 
% metrics.coh.D = squeeze(mean(metrics.coherency(2:7,:,:)));
% metrics.coh.Th = squeeze(mean(metrics.coherency(7:17,:,:)));
% metrics.coh.A = squeeze(mean(metrics.coherency(17:25,:,:)));
% metrics.coh.B = squeeze(mean(metrics.coherency(25:77,:,:)));
% metrics.coh.low = squeeze(mean(metrics.coherency(2:77,:,:)));
% metrics.coh_network.D = metrics.coh.D > 0.75;
% metrics.coh_network.Th = metrics.coh.Th > 0.75;
% metrics.coh_network.A = metrics.coh.A > 0.75;
% metrics.coh_network.B = metrics.coh.B > 0.75;
% metrics.coh_network.low = metrics.coh.low > 0.75;
% 
% metrics.coh_calm.D = squeeze(mean(metrics.coherency_calm(2:7,:,:)));
% metrics.coh_calm.Th = squeeze(mean(metrics.coherency_calm(7:17,:,:)));
% metrics.coh_calm.A = squeeze(mean(metrics.coherency_calm(17:25,:,:)));
% metrics.coh_calm.B = squeeze(mean(metrics.coherency_calm(25:77,:,:)));
% metrics.coh_calm.low = squeeze(mean(metrics.coherency_calm(2:77,:,:)));
% metrics.coh_network_calm.D = metrics.coh_calm.D > 0.75;
% metrics.coh_network_calm.Th = metrics.coh_calm.Th > 0.75;
% metrics.coh_network_calm.A = metrics.coh_calm.A > 0.75;
% metrics.coh_network_calm.B = metrics.coh_calm.B > 0.75;
% metrics.coh_network_calm.low = metrics.coh_calm.low > 0.75;
% 
% 
% metrics.coh_act.D = squeeze(mean(metrics.coherency_act(2:7,:,:)));
% metrics.coh_act.Th = squeeze(mean(metrics.coherency_act(7:17,:,:)));
% metrics.coh_act.A = squeeze(mean(metrics.coherency_act(17:25,:,:)));
% metrics.coh_act.B = squeeze(mean(metrics.coherency_act(25:77,:,:)));
% metrics.coh_act.low = squeeze(mean(metrics.coherency_act(2:77,:,:)));
% metrics.coh_network_act.D = metrics.coh_act.D > 0.75;
% metrics.coh_network_act.Th = metrics.coh_act.Th > 0.75;
% metrics.coh_network_act.A = metrics.coh_act.A > 0.75;
% metrics.coh_network_act.B = metrics.coh_act.B > 0.75;
% metrics.coh_network_act.low = metrics.coh_act.low > 0.75;

% metrics.coh_means_low =  squeeze(mean(metrics.coherency(2:40,:,:)));
% coh_means_low_calm  = squeeze(mean(metrics.coherency_calm(2:40,:,:)));
% coh_means_low_act  = squeeze(mean(metrics.coherency_act(2:40,:,:)));
% metrics.coh_means_high = squeeze(mean(metrics.coherency(end-40:end,:,:)));
% coh_low_norm = mean(metrics.coh_means_low);   
% coh_high_norm = mean(metrics.coh_means_high);

% %Coherence networks
% coh_bool_low = metrics.coh_means_low > 0.85;
% coh_bool_high = metrics.coh_means_low - metrics.coh_means_high > 0.1; %??
% metrics.coh_network = coh_bool_low & coh_bool_high;
% metrics.coh_network_calm = coh_means_low_calm > 0.85;
% metrics.coh_network_act = coh_means_low_act > 0.85;
% metrics.coh_missclassifications = -1;
% metrics.coh_confmat = zeros(2);
% metrics.coh_missclassifications_calm = -1;
% metrics.coh_confmat_calm = zeros(2);
% metrics.coh_missclassifications_act = -1;
% metrics.coh_confmat_act = zeros(2);
% metrics.coh_error_network = -1*ones(size(data,1));
% metrics.coh_error_network = -1*ones(size(data,1));
% metrics.coh_error_network_calm = -1*ones(size(data,1));
% metrics.coh_error_network_act = -1*ones(size(data,1));



%Compute channel RCGCIM 
alpha = 0.05; %threshold for rCGCIM
pmax = 1;
maketest = 1;
[metrics.RCGCIM,metrics.pRCGCIM]=mBTSCGCImatrix(data(:,5000:10000)',pmax,maketest);
[metrics.RCGCIM_calm,metrics.pRCGCIM_calm]=mBTSCGCImatrix(data(:,calm_start:calm_end)',pmax,maketest);
[metrics.RCGCIM_act,metrics.pRCGCIM_act]=mBTSCGCImatrix(data(:,act_start:act_end)',pmax,maketest);
metrics.RCGCIM_network = metrics.RCGCIM > alpha;
metrics.RCGCIM_missclassifications = -1;
metrics.RCGCIM_confmat = zeros(2);
metrics.RCGCIM_network_calm = metrics.RCGCIM_calm > alpha;
metrics.RCGCIM_missclassifications_calm = -1;
metrics.RCGCIM_confmat_calm = zeros(2);
metrics.RCGCIM_network_act = metrics.RCGCIM_act > alpha;
metrics.RCGCIM_missclassifications_act = -1;
metrics.RCGCIM_confmat_act = zeros(2);

