function [noiseidx maxcov_idx] = find_eyeblink_comp(eyeblink,ICA_components)
noiseidx = [];
maxcov_idx = -1;
maxcov = 0;
for i=1:size(ICA_components,1)
    covariance = cov(ICA_components(i,:),eyeblink);
    if abs(covariance(1,2)) > 5;
        noiseidx = [noiseidx i];
    end
    if maxcov < abs(covariance(1,2))
            maxcov = abs(covariance(1,2));
            maxcov_idx = i;
    end
end