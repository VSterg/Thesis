function metrics = compute_conf_mat_metrics(conf_mat)
metrics.accuracy = sum(diag(conf_mat))/sum(sum(conf_mat));
metrics.specificity = conf_mat(1,1)/(conf_mat(1,1)+conf_mat(1,2));
metrics.sensitivity = conf_mat(2,2)/(conf_mat(2,2)+conf_mat(2,1));