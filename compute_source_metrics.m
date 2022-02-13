function metrics = compute_source_metrics(source_data,real_source_data,Fs,calm_start,calm_end,act_start,act_end)
alpha = 0.1; %threshold for CGCIM
beta = 0.15; %threshold for DTF / PDC
coef = zeros(size(source_data,1),1);
for i=1:size(source_data,1)
        x1 = source_data(i,:);
        x2 = real_source_data(i,:);
        tmpcoef = corrcoef(x1,x2);
        coef(i) = tmpcoef(1,2);
end
metrics.source_corr = coef;
[CGCIM metrics.pCGCIM] = CGCinall(source_data',2,1);
[CGCIM_act metrics.pCGCIM_act] = CGCinall(source_data(:,act_start:act_end)',2,1);
[CGCIM_calm metrics.pCGCIM_calm] = CGCinall(source_data(:,calm_start:calm_end)',2,1);
metrics.CGCIM_network = CGCIM > alpha;
metrics.CGCIM_network_act = CGCIM_act > alpha;
metrics.CGCIM_network_calm = CGCIM_calm > alpha;
metrics.CGCIM = CGCIM;
metrics.CGCIM_calm = CGCIM_calm;
metrics.CGCIM_act = CGCIM_act;


%PDF - DTF
L = size(source_data,2);    % Number of samples
CH = size(source_data,1);   % Number of channels
c21 = zeros(1,L); % First time-varying parameter
c23 = zeros(1,L); % Second time-varying parameter
c12 = zeros(1,L); % Third time-varying parameter
N_freq = 50;      % Number of frequency bins
                  % It doesn't matter for the simulated model what the value of Fs is.
Fmax = Fs/2;      % Maximum frequency limit in the PDC and DTF 
[w, A_TI, C, sbc, fpe, th] = arfit(source_data', 1, 20, 'sbc'); % ---> ARFIT toolbox
[tmp,p_opt] = min(sbc); % Optimum order for the MVAR model using the SBC approach

UC = .0001; % Update Coefficient for the adaptive AR modelling
Down_Sample = 10; % Downsampling rate 
%%%%%% AAR estimation from the BioSig toolbox
Mode2 = 2; % Update mode of the process noise covariance matrix Q
[A_TV,e,Kalman,C_TV] = mvaar(source_data',p_opt,UC,Mode2); % Multivariate Adaptive AR estimation using Kalman filter --> BioSig toolbox
A_TV_reshape = reshape(A_TV', CH*p_opt, CH, size(A_TV,1)); % (CH*p x CH x T)
for i = 1 : size(A_TV_reshape,3)
    A_TV3(:,:,i) = A_TV_reshape(:,:,i)';
end
A_TV3 = A_TV3(:,:,1:Down_Sample:end); % Down sampling in the AR coefficients
T = size(A_TV3,3); % Number of time points after down sampling
[PDC_TV_all_freq, DTF_TV_all_freq] = PDC_DTF_matrix(A_TV3,p_opt,Fs,Fmax,N_freq); % Time-varying PDC and DTF
PDC_TV_all_freq = mean(PDC_TV_all_freq,4);
DTF_TV_all_freq = mean(DTF_TV_all_freq,4);

%frequency band delta
PDC_TV.D = PDC_TV_all_freq(:,:,1:2);
DTF_TV.D = DTF_TV_all_freq(:,:,1:2);
PDC_TV.D = mean(PDC_TV.D,3);
DTF_TV.D = mean(DTF_TV.D,3);
metrics.PDC_network.D = PDC_TV.D' > beta;
metrics.DTF_network.D = DTF_TV.D' > beta;
metrics.PDC_network.D = metrics.PDC_network.D - diag(diag(metrics.PDC_network.D));
metrics.DTF_network.D = metrics.DTF_network.D - diag(diag(metrics.DTF_network.D));
metrics.PDC_network.D = logical(metrics.PDC_network.D);
metrics.DTF_network.A = logical(metrics.DTF_network.D);

%frequency band theta
PDC_TV.Th = PDC_TV_all_freq(:,:,2:4);
DTF_TV.Th = DTF_TV_all_freq(:,:,2:4);
PDC_TV.Th = mean(PDC_TV.Th,3);
DTF_TV.Th = mean(DTF_TV.Th,3);
metrics.PDC_network.Th = PDC_TV.Th' > beta;
metrics.DTF_network.Th = DTF_TV.Th' > beta;
metrics.PDC_network.Th = metrics.PDC_network.Th - diag(diag(metrics.PDC_network.Th));
metrics.DTF_network.Th = metrics.DTF_network.Th - diag(diag(metrics.DTF_network.Th));
metrics.PDC_network.Th = logical(metrics.PDC_network.Th);
metrics.DTF_network.Th = logical(metrics.DTF_network.Th);

%frequency band alpha
PDC_TV.A = PDC_TV_all_freq(:,:,3:5);
DTF_TV.A = DTF_TV_all_freq(:,:,3:5);
PDC_TV.A = mean(PDC_TV.A,3);
DTF_TV.A = mean(DTF_TV.A,3);
metrics.PDC_network.A = PDC_TV.A' > beta;
metrics.DTF_network.A = DTF_TV.A' > beta;
metrics.PDC_network.A = metrics.PDC_network.A - diag(diag(metrics.PDC_network.A));
metrics.DTF_network.A = metrics.DTF_network.A - diag(diag(metrics.DTF_network.A));
metrics.PDC_network.A = logical(metrics.PDC_network.A);
metrics.DTF_network.A = logical(metrics.DTF_network.A);

%frequency band beta
PDC_TV.B = PDC_TV_all_freq(:,:,5:15);
DTF_TV.B = DTF_TV_all_freq(:,:,5:15);
PDC_TV.B = mean(PDC_TV.B,3);
DTF_TV.B = mean(DTF_TV.B,3);
metrics.PDC_network.B = PDC_TV.B' > beta;
metrics.DTF_network.B = DTF_TV.B' > beta;
metrics.PDC_network.B = metrics.PDC_network.B - diag(diag(metrics.PDC_network.B));
metrics.DTF_network.B = metrics.DTF_network.B - diag(diag(metrics.DTF_network.B));
metrics.PDC_network.B = logical(metrics.PDC_network.B);
metrics.DTF_network.B = logical(metrics.DTF_network.B);

%all low frequencies
PDC_TV.low = PDC_TV_all_freq(:,:,2:18);
DTF_TV.low = DTF_TV_all_freq(:,:,2:18);
PDC_TV.low = mean(PDC_TV.low,3);
DTF_TV.low = mean(DTF_TV.low,3);
metrics.PDC_network.low = PDC_TV.low' > beta;
metrics.DTF_network.low = DTF_TV.low' > beta;
metrics.PDC_network.low = metrics.PDC_network.low - diag(diag(metrics.PDC_network.low));
metrics.DTF_network.low = metrics.DTF_network.low - diag(diag(metrics.DTF_network.low));
metrics.PDC_network.low = logical(metrics.PDC_network.low);
metrics.DTF_network.low = logical(metrics.DTF_network.low);

source_data_calm = source_data(:,calm_start:calm_end);
[A_TV,e,Kalman,C_TV] = mvaar(source_data_calm',p_opt,UC,Mode2); % Multivariate Adaptive AR estimation using Kalman filter --> BioSig toolbox
A_TV_reshape = reshape(A_TV', CH*p_opt, CH, size(A_TV,1)); % (CH*p x CH x T)
for i = 1 : size(A_TV_reshape,3)
    A_TV3(:,:,i) = A_TV_reshape(:,:,i)';
end
A_TV3 = A_TV3(:,:,1:Down_Sample:end); % Down sampling in the AR coefficients
T = size(A_TV3,3); % Number of time points after down sampling
[PDC_TV_calm_all_freq, DTF_TV_calm_all_freq] = PDC_DTF_matrix(A_TV3,p_opt,Fs,Fmax,N_freq); % Time-varying PDC and DTF
PDC_TV_calm_all_freq = mean(PDC_TV_calm_all_freq,4);
DTF_TV_calm_all_freq = mean(DTF_TV_calm_all_freq,4);

%frequency band delta
PDC_TV_calm.D = PDC_TV_calm_all_freq(:,:,1:2);
DTF_TV_calm.D = DTF_TV_calm_all_freq(:,:,1:2);
PDC_TV_calm.D = mean(PDC_TV_calm.D,3);
DTF_TV_calm.D = mean(DTF_TV_calm.D,3);
metrics.PDC_network_calm.D = PDC_TV_calm.D' > beta;
metrics.DTF_network_calm.D = DTF_TV_calm.D' > beta;
metrics.PDC_network_calm.D = metrics.PDC_network_calm.D - diag(diag(metrics.PDC_network_calm.D));
metrics.DTF_network_calm.D = metrics.DTF_network_calm.D - diag(diag(metrics.DTF_network_calm.D));
metrics.PDC_network_calm.D = logical(metrics.PDC_network_calm.D);
metrics.DTF_network_calm.D = logical(metrics.DTF_network_calm.D);

%frequency band theta
PDC_TV_calm.Th = PDC_TV_calm_all_freq(:,:,2:4);
DTF_TV_calm.Th = DTF_TV_calm_all_freq(:,:,2:4);
PDC_TV_calm.Th = mean(PDC_TV_calm.Th,3);
DTF_TV_calm.Th = mean(DTF_TV_calm.Th,3);
metrics.PDC_network_calm.Th = PDC_TV_calm.Th' > beta;
metrics.DTF_network_calm.Th = DTF_TV_calm.Th' > beta;
metrics.PDC_network_calm.Th = metrics.PDC_network_calm.Th - diag(diag(metrics.PDC_network_calm.Th));
metrics.DTF_network_calm.Th = metrics.DTF_network_calm.Th - diag(diag(metrics.DTF_network_calm.Th));
metrics.PDC_network_calm.Th = logical(metrics.PDC_network_calm.Th);
metrics.DTF_network_calm.Th = logical(metrics.DTF_network_calm.Th);

%frequency band alpha
PDC_TV_calm.A = PDC_TV_calm_all_freq(:,:,3:5);
DTF_TV_calm.A = DTF_TV_calm_all_freq(:,:,3:5);
PDC_TV_calm.A = mean(PDC_TV_calm.A,3);
DTF_TV_calm.A = mean(DTF_TV_calm.A,3);
metrics.PDC_network_calm.A = PDC_TV_calm.A' > beta;
metrics.DTF_network_calm.A = DTF_TV_calm.A' > beta;
metrics.PDC_network_calm.A = metrics.PDC_network_calm.A - diag(diag(metrics.PDC_network_calm.A));
metrics.DTF_network_calm.A = metrics.DTF_network_calm.A - diag(diag(metrics.DTF_network_calm.A));
metrics.PDC_network_calm.A = logical(metrics.PDC_network_calm.A);
metrics.DTF_network_calm.A = logical(metrics.DTF_network_calm.A);

%frequency band beta
PDC_TV_calm.B = PDC_TV_calm_all_freq(:,:,5:15);
DTF_TV_calm.B = DTF_TV_calm_all_freq(:,:,5:15);
PDC_TV_calm.B = mean(PDC_TV_calm.B,3);
DTF_TV_calm.B = mean(DTF_TV_calm.B,3);
metrics.PDC_network_calm.B = PDC_TV_calm.B' > beta;
metrics.DTF_network_calm.B = DTF_TV_calm.B' > beta;
metrics.PDC_network_calm.B = metrics.PDC_network_calm.B - diag(diag(metrics.PDC_network_calm.B));
metrics.DTF_network_calm.B = metrics.DTF_network_calm.B - diag(diag(metrics.DTF_network_calm.B));
metrics.PDC_network_calm.B = logical(metrics.PDC_network_calm.B);
metrics.DTF_network_calm.B = logical(metrics.DTF_network_calm.B);

%all low frequencies
PDC_TV_calm.low = PDC_TV_calm_all_freq(:,:,2:18);
DTF_TV_calm.low = DTF_TV_calm_all_freq(:,:,2:18);
PDC_TV_calm.low = mean(PDC_TV_calm.low,3);
DTF_TV_calm.low = mean(DTF_TV_calm.low,3);
metrics.PDC_network_calm.low = PDC_TV_calm.low' > beta;
metrics.DTF_network_calm.low = DTF_TV_calm.low' > beta;
metrics.PDC_network_calm.low = metrics.PDC_network_calm.low - diag(diag(metrics.PDC_network_calm.low));
metrics.DTF_network_calm.low = metrics.DTF_network_calm.low - diag(diag(metrics.DTF_network_calm.low));
metrics.PDC_network_calm.low = logical(metrics.PDC_network_calm.low);
metrics.DTF_network_calm.low = logical(metrics.DTF_network_calm.low);

source_data_act = source_data(:,act_start:act_end);
[A_TV,e,Kalman,C_TV] = mvaar(source_data_act',p_opt,UC,Mode2); % Multivariate Adaptive AR estimation using Kalman filter --> BioSig toolbox
A_TV_reshape = reshape(A_TV', CH*p_opt, CH, size(A_TV,1)); % (CH*p x CH x T)
for i = 1 : size(A_TV_reshape,3)
    A_TV3(:,:,i) = A_TV_reshape(:,:,i)';
end
A_TV3 = A_TV3(:,:,1:Down_Sample:end); % Down sampling in the AR coefficients
T = size(A_TV3,3); % Number of time points after down sampling
[PDC_TV_act_all_freq, DTF_TV_act_all_freq] = PDC_DTF_matrix(A_TV3,p_opt,Fs,Fmax,N_freq); % Time-varying PDC and DTF
PDC_TV_act_all_freq = mean(PDC_TV_act_all_freq,4);
DTF_TV_act_all_freq = mean(DTF_TV_act_all_freq,4);

%frequency band delta
PDC_TV_act.D = PDC_TV_act_all_freq(:,:,1:2);
DTF_TV_act.D = DTF_TV_act_all_freq(:,:,1:2);
PDC_TV_act.D = mean(PDC_TV_act.D,3);
DTF_TV_act.D = mean(DTF_TV_act.D,3);
metrics.PDC_network_act.D = PDC_TV_act.D' > beta;
metrics.DTF_network_act.D = DTF_TV_act.D' > beta;
metrics.PDC_network_act.D = metrics.PDC_network_act.D - diag(diag(metrics.PDC_network_act.D));
metrics.DTF_network_act.D = metrics.DTF_network_act.D - diag(diag(metrics.DTF_network_act.D));
metrics.PDC_network_act.D = logical(metrics.PDC_network_act.D);
metrics.DTF_network_act.D = logical(metrics.DTF_network_act.D);

%frequency band theta
PDC_TV_act.Th = PDC_TV_act_all_freq(:,:,2:4);
DTF_TV_act.Th = DTF_TV_act_all_freq(:,:,2:4);
PDC_TV_act.Th = mean(PDC_TV_act.Th,3);
DTF_TV_act.Th = mean(DTF_TV_act.Th,3);
metrics.PDC_network_act.Th = PDC_TV_act.Th' > beta;
metrics.DTF_network_act.Th = DTF_TV_act.Th' > beta;
metrics.PDC_network_act.Th = metrics.PDC_network_act.Th - diag(diag(metrics.PDC_network_act.Th));
metrics.DTF_network_act.Th = metrics.DTF_network_act.Th - diag(diag(metrics.DTF_network_act.Th));
metrics.PDC_network_act.Th = logical(metrics.PDC_network_act.Th);
metrics.DTF_network_act.Th = logical(metrics.DTF_network_act.Th);

%frequency band alpha
PDC_TV_act.A = PDC_TV_act_all_freq(:,:,3:5);
DTF_TV_act.A = DTF_TV_act_all_freq(:,:,3:5);
PDC_TV_act.A = mean(PDC_TV_act.A,3);
DTF_TV_act.A = mean(DTF_TV_act.A,3);
metrics.PDC_network_act.A = PDC_TV_act.A' > beta;
metrics.DTF_network_act.A = DTF_TV_act.A' > beta;
metrics.PDC_network_act.A = metrics.PDC_network_act.A - diag(diag(metrics.PDC_network_act.A));
metrics.DTF_network_act.A = metrics.DTF_network_act.A - diag(diag(metrics.DTF_network_act.A));
metrics.PDC_network_act.A = logical(metrics.PDC_network_act.A);
metrics.DTF_network_act.A = logical(metrics.DTF_network_act.A);

%frequency band beta
PDC_TV_act.B = PDC_TV_act_all_freq(:,:,5:15);
DTF_TV_act.B = DTF_TV_act_all_freq(:,:,5:15);
PDC_TV_act.B = mean(PDC_TV_act.B,3);
DTF_TV_act.B = mean(DTF_TV_act.B,3);
metrics.PDC_network_act.B = PDC_TV_act.B' > beta;
metrics.DTF_network_act.B = DTF_TV_act.B' > beta;
metrics.PDC_network_act.B = metrics.PDC_network_act.B - diag(diag(metrics.PDC_network_act.B));
metrics.DTF_network_act.B = metrics.DTF_network_act.B - diag(diag(metrics.DTF_network_act.B));
metrics.PDC_network_act.B = logical(metrics.PDC_network_act.B);
metrics.DTF_network_act.B = logical(metrics.DTF_network_act.B);

%all low frequencies
PDC_TV_act.low = PDC_TV_act_all_freq(:,:,2:18);
DTF_TV_act.low = DTF_TV_act_all_freq(:,:,2:18);
PDC_TV_act.low = mean(PDC_TV_act.low,3);
DTF_TV_act.low = mean(DTF_TV_act.low,3);
metrics.PDC_network_act.low = PDC_TV_act.low' > beta;
metrics.DTF_network_act.low = DTF_TV_act.low' > beta;
metrics.PDC_network_act.low = metrics.PDC_network_act.low - diag(diag(metrics.PDC_network_act.low));
metrics.DTF_network_act.low = metrics.DTF_network_act.low - diag(diag(metrics.DTF_network_act.low));
metrics.PDC_network_act.low = logical(metrics.PDC_network_act.low);
metrics.DTF_network_act.low = logical(metrics.DTF_network_act.low);

%initialiase
metrics.CGCIM_missclassifications = -1;
metrics.CGCIM_missclassifications_calm = -1;
metrics.CGCIM_missclassifications_act = -1;
metrics.CGCIM_confmat = zeros(2);
metrics.CGCIM_confmat_act = zeros(2);
metrics.CGCIM_confmat_calm = zeros(2);

metrics.DTF_missclassifications.D = -1;
metrics.DTF_missclassifications_calm.D = -1;
metrics.DTF_missclassifications_act.D = -1;
metrics.DTF_confmat.D = zeros(2);
metrics.DTF_confmat_act.D = zeros(2);
metrics.DTF_confmat_calm.D = zeros(2);
metrics.DTF_missclassifications.Th = -1;
metrics.DTF_missclassifications_calm.Th = -1;
metrics.DTF_missclassifications_act.Th = -1;
metrics.DTF_confmat.Th = zeros(2);
metrics.DTF_confmat_act.Th = zeros(2);
metrics.DTF_confmat_calm.Th = zeros(2);
metrics.DTF_missclassifications.A = -1;
metrics.DTF_missclassifications_calm.A = -1;
metrics.DTF_missclassifications_act.A = -1;
metrics.DTF_confmat.A = zeros(2);
metrics.DTF_confmat_act.A = zeros(2);
metrics.DTF_confmat_calm.A = zeros(2);
metrics.DTF_missclassifications.B = -1;
metrics.DTF_missclassifications_calm.B = -1;
metrics.DTF_missclassifications_act.B = -1;
metrics.DTF_confmat.B = zeros(2);
metrics.DTF_confmat_act.B = zeros(2);
metrics.DTF_confmat_calm.B = zeros(2);
metrics.DTF_missclassifications.low = -1;
metrics.DTF_missclassifications_calm.low = -1;
metrics.DTF_missclassifications_act.low = -1;
metrics.DTF_confmat.low = zeros(2);
metrics.DTF_confmat_act.low = zeros(2);
metrics.DTF_confmat_calm.low = zeros(2);

metrics.PDC_missclassifications.D = -1;
metrics.PDC_missclassifications_calm.D = -1;
metrics.PDC_missclassifications_act.D = -1;
metrics.PDC_confmat.D = zeros(2);
metrics.PDC_confmat_act.D = zeros(2);
metrics.PDC_confmat_calm.D = zeros(2);
metrics.PDC_missclassifications.Th = -1;
metrics.PDC_missclassifications_calm.Th = -1;
metrics.PDC_missclassifications_act.Th = -1;
metrics.PDC_confmat.Th = zeros(2);
metrics.PDC_confmat_act.Th = zeros(2);
metrics.PDC_confmat_calm.Th = zeros(2);
metrics.PDC_missclassifications.A = -1;
metrics.PDC_missclassifications_calm.A = -1;
metrics.PDC_missclassifications_act.A = -1;
metrics.PDC_confmat.A = zeros(2);
metrics.PDC_confmat_act.A = zeros(2);
metrics.PDC_confmat_calm.A = zeros(2);
metrics.PDC_missclassifications.B = -1;
metrics.PDC_missclassifications_calm.B = -1;
metrics.PDC_missclassifications_act.B = -1;
metrics.PDC_confmat.B = zeros(2);
metrics.PDC_confmat_act.B = zeros(2);
metrics.PDC_confmat_calm.B = zeros(2);