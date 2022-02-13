function comparison = compare_channel_metrics(metrics_1,metrics_2,act_start,act_end)
%values higher than 0 = > metrics 1 is better than metrics 2
comparison.error = metrics_2.channel_mean_sq_diff - metrics_1.channel_mean_sq_diff;
comparison.error_calm = metrics_2.channel_mean_sq_diff_calm - metrics_1.channel_mean_sq_diff_calm;
comparison.error_act = metrics_2.channel_mean_sq_diff_act - metrics_1.channel_mean_sq_diff_act;

%values higher than 1 = > metrics 1 is better than metrics 2
relative_error = metrics_2.channel_mean_sq_diff./metrics_1.channel_mean_sq_diff;
comparison.relative_error = log(relative_error);
relative_error = metrics_2.channel_mean_sq_diff_calm./metrics_1.channel_mean_sq_diff_calm ;
comparison.relative_error_calm = log(relative_error);
relative_error = metrics_2.channel_mean_sq_diff_act./metrics_1.channel_mean_sq_diff_act ;
comparison.relative_error_act = log(relative_error);

act_error1 = mean(abs(metrics_1.channel_diff(:,act_start:act_end)),2);
act_error2 = mean(abs(metrics_2.channel_diff(:,act_start:act_end)),2);
comparison.error_pct =100*( act_error1./act_error2);

% A = metrics_1.coh_error_network;
% B = 2*metrics_2.coh_error_network;
% comparison.coh_errors = A + B;
% 
% A = metrics_1.coh_error_network_calm;
% B = 2*metrics_2.coh_error_network_calm;
% comparison.coh_errors_calm = A + B;
% 
% A = metrics_1.coh_error_network_act;
% B = 2*metrics_2.coh_error_network_act;
% comparison.coh_errors_act = A + B;