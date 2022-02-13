function mat = extend_mat(Adj_mat,active_regions,desired_size)
mat = zeros(desired_size);
for i=1:length(active_regions)
    for j=1:length(active_regions)
        mat(active_regions(i),active_regions(j)) = Adj_mat(i,j);
    end
end

