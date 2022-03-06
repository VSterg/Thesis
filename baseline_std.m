function baseline_std = baseline_std(timeseries,base)
N = size(timeseries,1);
timeseries = timeseries - base;
S = sum(timeseries.*timeseries);
S = S/N;
baseline_std = sqrt(S);