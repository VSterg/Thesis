%% RCGCIM threshold compute
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

%% compare second level source metrics
clear metric_to_count
level = 'source';
%remove/add 'clean' to switch between 7-ROI and 68-ROI level metrics
cleaning = {'clean'; 'noisy'; 'ICA'; 'wICA'; 'ICA_strips'; 'ICA_A_rebuilt'};
%switch between RCGCIM - CGCIM and _act - _calm - 
metric = 'PDC_confmat_act_metrics';
%switch between sensitivity - specificity - accuracy
metric2 = 'sensitivity';
for i = 1:100
    filename = ['Pink_data_analysis_' num2str(i,'%02d')]
    A3 = load(filename);
    for j = 1:length(cleaning)
        metrics_from = [level '_metrics_' cleaning(j)];
        %remove/add _ext from line below to switch between 7-ROI and 68-ROI level metrics.
        metrics_from = string([metrics_from{1} metrics_from{2} metrics_from{3}]);
        metric_to_count(i,j,:,:) = A3.(metrics_from).(metric).low.(metric2);
    end
end
figure
boxplot(metric_to_count,'Labels',{'Clean','Noisy','ICA','wICA','sICA','MNA'})
%% Compare first level source metrics
clear metric_to_count
level = 'source';
cleaning = {'clean'; 'noisy'; 'ICA'; 'wICA'; 'ICA_strips'; 'ICA_A_rebuilt'};
metric = 'RCGCIM_confmat_act';
for i = 1:100
    filename = ['Pink_data_analysis_' num2str(i,'%02d')]
    A3 = load(filename);
    for j = 1:length(cleaning)
        metrics_from = [level '_metrics_' cleaning(j)];
        metrics_from = string([metrics_from{1} metrics_from{2} metrics_from{3} '_ext']);
        tmp = A3.(metrics_from).(metric)
        %metric_to_count(i,j,:,:) = A3.(metrics_from).(metric).low;
        metric_to_count(i,j,:,:) = tmp(1,2);
    end
end
figure
boxplot(metric_to_count,'Labels',{'Clean','Noisy','ICA','wICA','sICA','MNA'})
%% Compare first level source metrics from second batch
clear metric_to_count
level = 'source';
cleaning = {'clean'; 'noisy'; 'ICA'; 'wICA'; 'ICA_strips'; 'ICA_A_rebuilt'};
metric = 'RCGCIM_missclassifications_act';
for j = 1:6
    conf_all(:,:,j) = [0 0; 0 0]
end
for i = 1:100
    filename = ['Pink_data_analysis2_' num2str(i,'%02d')]
    A3 = load(filename).A3;
    for j = 1:length(cleaning)
        metrics_from = [level '_metrics_' cleaning(j)];
        metrics_from = string([metrics_from{1} metrics_from{2} metrics_from{3} '_ext']);
        %tmp = A3.(metrics_from).(metric);
        %conf_all(:,:,j) = conf_all(:,:,j) + tmp;
        metric_to_count(i,j,:,:) = A3.(metrics_from).(metric);
    end
end
%figure
%boxplot(metric_to_count,'Labels',{'Clean','Noisy','ICA','wICA','sICA','MNA'})
%% compare second level source metrics from second batch
clear metric_to_count
level = 'source';
%remove/add 'clean' to switch between 7-ROI and 68-ROI level metrics
cleaning = {'clean'; 'noisy'; 'ICA'; 'wICA'; 'ICA_strips'; 'ICA_A_rebuilt'};
%switch between RCGCIM - CGCIM and _act - _calm - 
metric = 'RCGCIM_confmat_act_metrics';
%switch between sensitivity - specificity - accuracy
metric2 = 'sensitivity';
for i = 1:100
    filename = ['Pink_data_analysis2_' num2str(i,'%02d')]
    A3 = load(filename).A3;
    for j = 1:length(cleaning)
        metrics_from = [level '_metrics_' cleaning(j)];
        %remove/add _ext from line below to switch between 7-ROI and 68-ROI level metrics.
        metrics_from = string([metrics_from{1} metrics_from{2} metrics_from{3} '_ext']);
        metric_to_count(i,j,:,:) = A3.(metrics_from).(metric).(metric2);
    end
end
%% Compare first level channel metrics
clear metric_to_count
level = 'channel'
cleaning = {'clean'; 'noisy'; 'ICA'; 'wICA'; 'ICA_strips'; 'ICA_A_rebuilt'};
metric = 'channel_mean_sq_diff_act';
error_all = zeros(64,64,6);
strip.a = [1 33 34 2 3 37 36 35]; %front of the head strip
strip.b = [4:11 38:47];
strip.c = [12:19 32 48:56];
strip.d = [20:26 30 31 57:63];
strip.e = [27 28 29 64]; %back of the head strip
names = ['a' 'b' 'c' 'd' 'e'];

for i = 1:100
    filename = ['Pink_data_analysis2_' num2str(i,'%02d')]
    A3 = load(filename).A3;
    for j = 1:length(cleaning)
        metrics_from = [level '_metrics_' cleaning(j)];
        metrics_from = string([metrics_from{1} metrics_from{2} metrics_from{3}]);
        metric_to_count(i,j,:,:) = A3.(metrics_from).(metric);
        %tmp = A3.(metrics_from).(metric)
        %error_all(:,:,j) = error_all(:,:,j) + tmp;
    end
end
metric_mean = squeeze(mean(metric_to_count,1));
for i=1:5
    current_strip = names(i);
    current_strip = strip.(current_strip);
    strip_metric_mean(:,i) = mean(metric_mean(:,current_strip),2);

end
% for j = 1:5
%     a1 = reshape(squeeze(error_all(:,:,j)),[64*64,1]);
%     [a(:,j) I(:,j)] = maxk(a1,50);
% end
% I2 = floor(I/64)+1;
% I1 = mod(I,64)
% most_common_errors= zeros(64,64,5);
% for i = 1:5
%     for j = 1:50
%      most_common_errors(I1(j,i),I2(j,i),i) = 1;
%     end
% end
%% compare second level channel metrics
clear metric_to_count
level = 'channel';
cleaning = {'clean';'noisy'; 'ICA'; 'wICA'; 'ICA_strips'; 'ICA_A_rebuilt'};
%switch between and _act - _calm - 
metric = 'RCGCIM_confmat_act_metrics';
%switch between sensitivity - specificity - accuracy
metric2 = 'sensitivity';
for j = 1:6
    conf_all(:,:,j) = [0 0; 0 0]
end
for i = 1:100
    filename = ['Pink_data_analysis2_' num2str(i,'%02d')]
    A3 = load(filename).A3;
    for j = 1:length(cleaning)
        metrics_from = [level '_metrics_' cleaning(j)];
        %remove/add _ext from line below to switch between 7-ROI and 68-ROI level metrics.
        metrics_from = string([metrics_from{1} metrics_from{2} metrics_from{3}]);
        metric_to_count(i,j,:,:) = A3.(metrics_from).(metric).(metric2);
        %tmp = A3.(metrics_from).(metric);
        %conf_all(:,:,j) = conf_all(:,:,j) + tmp;
    end
end
%%
%% Compare localization error
clear metric_to_count
level = 'channel'
error_all = zeros(64,64,6);
strip.a = [1 33 34 2 3 37 36 35]; %front of the head strip
strip.b = [4:11 38:47];
strip.c = [12:19 32 48:56];
strip.d = [20:26 30 31 57:63];
strip.e = [27 28 29 64]; %back of the head strip
names = ['a' 'b' 'c' 'd' 'e'];
cleaning = {'ICA'; 'ICA_strips'; 'ICA_A_rebuilt'};
metric = 'localization_error2';
for i = 1:100
    filename = ['Pink_data_analysis2_' num2str(i,'%02d')]
    A3 = load(filename).A3;
    for j = 1:length(cleaning)
        metrics_from = [level '_metrics_' cleaning(j)];
        metrics_from = string([metrics_from{1} metrics_from{2} metrics_from{3}]);
        metric_to_count(i,j,:,:) = A3.(metrics_from).(metric);
    end
end
metric_to_count = squeeze(metric_to_count);
metric_to_count_mean = squeeze(mean(abs(metric_to_count),1));
metric_to_count_mean2 = squeeze(mean(metric_to_count,1));


for i=1:5
    current_strip = names(i);
    current_strip = strip.(current_strip);
    strip_metric_mean(:,i) = mean(metric_to_count_mean(:,current_strip),2);

end
figure
boxplot(metric_to_count_mean','Labels',{'ICA', 'sICA', 'NMA'})   
%% Adjacency matrix for visualization with BrainNet
active_regions = [1 2 5 6 7 9 33];
ROI_mapping = [1:2:9 13:2:17 21:2:31 35 33 37:2:63 11 65 67 19 2:2:10 14:2:18 22:2:32 36 34 38:2:64 12 66 68 20];
Adj_mat_ext = A3.source_metrics_clean_ext.RCGCIM_network;
%Adj_mat_ext = extend_mat(Adj_mat_ext,active_regions,68);
for i=1:68
    for j = 1:68
        Adj_mat_ext2(i,j) = Adj_mat_ext(ROI_mapping(i),ROI_mapping(j));
    end
end
writematrix(Adj_mat_ext2,'tabledata2.txt','Delimiter','\t')
%% plot electrode networks
aij64 = A3.channel_metrics_clean.RCGCIM_network;
figure(4);
f_PlotEEG_BrainNetwork(64)

nch = 64; %take care of this variable (nch must be according to matrix size you want to plot it)
p = 0.03;   %proportion of weigthed links to keep for.
aij = threshold_proportional(aij64, p); %thresholding networks due to proportion p
ijw = adj2edgeL(triu(aij));             %passing from matrix form to edge list form
n_features = sum(aij, 2);       % in this case the feature is the Strenght
cbtype = 'ncb';                 % colorbar for weigth
figure(5);
f_PlotEEG_BrainNetwork(nch, ijw, 'w_intact');
%% fix electrode
elec_adj = A3.channel_metrics_clean.RCGCIM_network;
figure()
h = plotnetworktitle(elec_adj)
%% make electrode graph plot
info = most_common_errors(:,:,1);
figure
G = digraph(info);
p = plot(G,'layout','circle')
p.XData = [EEG.chanlocs(:).X];
p.YData = [EEG.chanlocs(:).Y];
p.ZData = [EEG.chanlocs(:).Z];
%% change CGCIM threshold


clear all
[~,Adj_mat] = VAR1RingStructure(5000,7); %generate simulated data
Adj_mat = Adj_mat > 0.05;
Adj_mat = Adj_mat';

level = 'source'
cleaning = {'clean'; 'noisy'; 'ICA'; 'wICA'; 'ICA_strips'; 'ICA_A_rebuilt'};
for i = 2:100
    filename = ['Pink_data_analysis2_' num2str(i,'%02d')]
    A3 = load(filename).A1;
    for j = 1:length(cleaning)
        metrics_from = [level '_metrics_' cleaning(j)];
        metrics_from = string([metrics_from{1} metrics_from{2} metrics_from{3}]);
        A3.(metrics_from).CGCIM_network = A3.(metrics_from).CGCIM > 0.05;
        A3.(metrics_from).CGCIM_network_calm = A3.(metrics_from).CGCIM_calm > 0.05;
        A3.(metrics_from).CGCIM_network_act = A3.(metrics_from).CGCIM_act > 0.05;
        [A3.(metrics_from).CGCIM_missclassifications A3.(metrics_from).CGCIM_confmat A3.(metrics_from).CGCIM_confmat_metrics A3.(metrics_from).CGCIM_error_network] = compute_network_variance(Adj_mat ,A3.(metrics_from).CGCIM_network);
        [A3.(metrics_from).CGCIM_missclassifications_calm A3.(metrics_from).CGCIM_confmat_calm A3.(metrics_from).CGCIM_confmat_calm_metrics A3.(metrics_from).CGCIM_error_network_calm] = compute_network_variance(Adj_mat,A3.(metrics_from).CGCIM_network_calm);
        [A3.(metrics_from).CGCIM_missclassifications_act A3.(metrics_from).CGCIM_confmat_act A3.(metrics_from).CGCIM_confmat_act_metrics A3.(metrics_from).CGCIM_error_network_act] = compute_network_variance(Adj_mat,A3.(metrics_from).CGCIM_network_act);

    end
    save(filename,'A3')
end
