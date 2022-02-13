function xM = extend_ts(xM,active_regions,desired_size)
tmp_ts = zeros(size(xM,1),desired_size);
for i=1:length(active_regions)
    tmp_ts(:,active_regions(i)) = xM(:,i);
end
xM = tmp_ts;
