function z_adj = apply_adjusted_z(path_table, z_disp_matrix)
%--------------------------------------------------------------------------
% Calculated adjusted z positions from adjustment matrix.
%--------------------------------------------------------------------------

% Start from lowest z position 
c = min(path_table.z)-1;
path_table.z_adj = path_table.z-c;

% Apply z adjustments from adjustment matrix
[nrows, ncols] = size(z_disp_matrix);
for i = 1:nrows
    for j = 1:ncols
        idx = path_table.y == i & path_table.x == j;
        path_table(idx,:).z_adj = path_table(idx,:).z_adj - z_disp_matrix(i,j);
    end
end

% Remove tiles where not all z positions are present
v = sum(path_table.z_adj==path_table.z_adj');
path_table(v ~= max(v),:) = [];
path_table.z_adj = path_table.z_adj-min(path_table.z_adj)+1;
z_adj = [path_table(:,1) path_table(:,end)];

end