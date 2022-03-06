function metrics = compute_source_metrics_ext(source_data,Fs,calm_start,calm_end,act_start,act_end)
alpha = 0.1; %threshold for CGCIM
beta = 0.25; %threshold for RCGCIM


pmax = 1;
maketest = 1;
[metrics.RCGCIM,metrics.pRCGCIM]=mBTSCGCImatrix(source_data(:,5000:10000)',pmax,maketest);
[metrics.RCGCIM_calm,metrics.pRCGCIM_calm]=mBTSCGCImatrix(source_data(:,calm_start:calm_end)',pmax,maketest);
[metrics.RCGCIM_act,metrics.pRCGCIM_act]=mBTSCGCImatrix(source_data(:,act_start:act_end)',pmax,maketest);
metrics.RCGCIM_network = metrics.RCGCIM > beta;
metrics.RCGCIM_network_calm = metrics.RCGCIM_calm > beta;
metrics.RCGCIM_network_act = metrics.RCGCIM_act > beta;



% initialiase
% metrics.CGCIM_missclassifications = -1;
% metrics.CGCIM_missclassifications_calm = -1;
% metrics.CGCIM_missclassifications_act = -1;
% metrics.CGCIM_confmat = zeros(2);
% metrics.CGCIM_confmat_act = zeros(2);
% metrics.CGCIM_confmat_calm = zeros(2);

metrics.RCGCIM_missclassifications = -1;
metrics.RCGCIM_missclassifications_calm = -1;
metrics.RCGCIM_missclassifications_act = -1;
metrics.RCGCIM_confmat = zeros(2);
metrics.RCGCIM_confmat_calm = zeros(2);
metrics.RCGCIM_confmat_act = zeros(2);

% tic
% %PDF - DTF
% L = size(source_data,2);    % Number of samples
% CH = size(source_data,1);   % Number of channels
% c21 = zeros(1,L); % First time-varying parameter
% c23 = zeros(1,L); % Second time-varying parameter
% c12 = zeros(1,L); % Third time-varying parameter
% N_freq = 50;      % Number of frequency bins
%                   % It doesn't matter for the simulated model what the value of Fs is.
% Fmax = Fs/2;      % Maximum frequency limit in the PDC and DTF 
% [w, A_TI, C, sbc, fpe, th] = arfit(source_data', 1, 20, 'sbc'); % ---> ARFIT toolbox
% [tmp,p_opt] = min(sbc); % Optimum order for the MVAR model using the SBC approach
% 
% UC = .0001; % Update Coefficient for the adaptive AR modelling
% Down_Sample = 10; % Downsampling rate 
% %%%%%% AAR estimation from the BioSig toolbox
% Mode2 = 2; % Update mode of the process noise covariance matrix Q
% source_data_act = source_data(:,act_start:act_end);
% [A_TV,e,Kalman,C_TV] = mvaar(source_data_act',p_opt,UC,Mode2); % Multivariate Adaptive AR estimation using Kalman filter --> BioSig toolbox
% A_TV_reshape = reshape(A_TV', CH*p_opt, CH, size(A_TV,1)); % (CH*p x CH x T)
% for i = 1 : size(A_TV_reshape,3)
%     A_TV3(:,:,i) = A_TV_reshape(:,:,i)';
% end
% A_TV3 = A_TV3(:,:,1:Down_Sample:end); % Down sampling in the AR coefficients
% T = size(A_TV3,3); % Number of time points after down sampling
% [PDC_TV_act_all_freq, DTF_TV_act_all_freq] = PDC_DTF_matrix(A_TV3,p_opt,Fs,Fmax,N_freq); % Time-varying PDC and DTF
% PDC_TV_act_all_freq = mean(PDC_TV_act_all_freq,4);
% DTF_TV_act_all_freq = mean(DTF_TV_act_all_freq,4);
% 
% %frequency band delta
% PDC_TV_act.D = PDC_TV_act_all_freq(:,:,1:2);
% DTF_TV_act.D = DTF_TV_act_all_freq(:,:,1:2);
% PDC_TV_act.D = mean(PDC_TV_act.D,3);
% DTF_TV_act.D = mean(DTF_TV_act.D,3);
% metrics.PDC_network_act.D = PDC_TV_act.D' > beta;
% metrics.DTF_network_act.D = DTF_TV_act.D' > beta;
% metrics.PDC_network_act.D = metrics.PDC_network_act.D - diag(diag(metrics.PDC_network_act.D));
% metrics.DTF_network_act.D = metrics.DTF_network_act.D - diag(diag(metrics.DTF_network_act.D));
% metrics.PDC_network_act.D = logical(metrics.PDC_network_act.D);
% metrics.DTF_network_act.D = logical(metrics.DTF_network_act.D);
% 
% %frequency band theta
% PDC_TV_act.Th = PDC_TV_act_all_freq(:,:,2:4);
% DTF_TV_act.Th = DTF_TV_act_all_freq(:,:,2:4);
% PDC_TV_act.Th = mean(PDC_TV_act.Th,3);
% DTF_TV_act.Th = mean(DTF_TV_act.Th,3);
% metrics.PDC_network_act.Th = PDC_TV_act.Th' > beta;
% metrics.DTF_network_act.Th = DTF_TV_act.Th' > beta;
% metrics.PDC_network_act.Th = metrics.PDC_network_act.Th - diag(diag(metrics.PDC_network_act.Th));
% metrics.DTF_network_act.Th = metrics.DTF_network_act.Th - diag(diag(metrics.DTF_network_act.Th));
% metrics.PDC_network_act.Th = logical(metrics.PDC_network_act.Th);
% metrics.DTF_network_act.Th = logical(metrics.DTF_network_act.Th);
% 
% %frequency band alpha
% PDC_TV_act.A = PDC_TV_act_all_freq(:,:,3:5);
% DTF_TV_act.A = DTF_TV_act_all_freq(:,:,3:5);
% PDC_TV_act.A = mean(PDC_TV_act.A,3);
% DTF_TV_act.A = mean(DTF_TV_act.A,3);
% metrics.PDC_network_act.A = PDC_TV_act.A' > beta;
% metrics.DTF_network_act.A = DTF_TV_act.A' > beta;
% metrics.PDC_network_act.A = metrics.PDC_network_act.A - diag(diag(metrics.PDC_network_act.A));
% metrics.DTF_network_act.A = metrics.DTF_network_act.A - diag(diag(metrics.DTF_network_act.A));
% metrics.PDC_network_act.A = logical(metrics.PDC_network_act.A);
% metrics.DTF_network_act.A = logical(metrics.DTF_network_act.A);
% 
% %frequency band beta
% PDC_TV_act.B = PDC_TV_act_all_freq(:,:,5:15);
% DTF_TV_act.B = DTF_TV_act_all_freq(:,:,5:15);
% PDC_TV_act.B = mean(PDC_TV_act.B,3);
% DTF_TV_act.B = mean(DTF_TV_act.B,3);
% metrics.PDC_network_act.B = PDC_TV_act.B' > beta;
% metrics.DTF_network_act.B = DTF_TV_act.B' > beta;
% metrics.PDC_network_act.B = metrics.PDC_network_act.B - diag(diag(metrics.PDC_network_act.B));
% metrics.DTF_network_act.B = metrics.DTF_network_act.B - diag(diag(metrics.DTF_network_act.B));
% metrics.PDC_network_act.B = logical(metrics.PDC_network_act.B);
% metrics.DTF_network_act.B = logical(metrics.DTF_network_act.B);
% 
% %all low frequencies
% PDC_TV_act.low = PDC_TV_act_all_freq(:,:,2:18);
% DTF_TV_act.low = DTF_TV_act_all_freq(:,:,2:18);
% PDC_TV_act.low = mean(PDC_TV_act.low,3);
% DTF_TV_act.low = mean(DTF_TV_act.low,3);
% metrics.PDC_network_act.low = PDC_TV_act.low' > beta;
% metrics.DTF_network_act.low = DTF_TV_act.low' > beta;
% metrics.PDC_network_act.low = metrics.PDC_network_act.low - diag(diag(metrics.PDC_network_act.low));
% metrics.DTF_network_act.low = metrics.DTF_network_act.low - diag(diag(metrics.DTF_network_act.low));
% metrics.PDC_network_act.low = logical(metrics.PDC_network_act.low);
% metrics.DTF_network_act.low = logical(metrics.DTF_network_act.low);
% 
% 
% metrics.DTF_missclassifications_act.D = -1;
% metrics.DTF_confmat_act.D = zeros(2);
% metrics.DTF_missclassifications_act.Th = -1;
% metrics.DTF_confmat_act.Th = zeros(2);
% metrics.DTF_missclassifications_act.A = -1;
% metrics.DTF_confmat_act.A = zeros(2);
% metrics.DTF_missclassifications_act.B = -1;
% metrics.DTF_confmat_act.B = zeros(2);
% metrics.DTF_missclassifications_act.low = -1;
% metrics.DTF_confmat_act.low = zeros(2);
% metrics.PDC_missclassifications_act.D = -1;
% metrics.PDC_confmat_act.D = zeros(2);
% metrics.PDC_missclassifications_act.Th = -1;
% metrics.PDC_confmat_act.Th = zeros(2);
% metrics.PDC_missclassifications_act.A = -1;
% metrics.PDC_confmat_act.A = zeros(2);
% metrics.PDC_missclassifications_act.B = -1;
% metrics.PDC_confmat_act.B = zeros(2);
% toc