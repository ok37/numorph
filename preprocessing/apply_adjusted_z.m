function z_adj = apply_adjusted_z(path_table, z_disp_matrix)
    %%% Add adjusted z values based on overlapping results
    c = min(path_table.z)-1;
    path_table.z_adj = path_table.z-c;

    for i = 1:height(path_table)
        z_adj = z_disp_matrix(path_table(i,:).y,path_table(i,:).x);
        path_table.z_adj(i) = path_table.z_adj(i) - z_adj;
    end

    % Remove tiles where not all z positions are present
    v = sum(path_table.z_adj==path_table.z_adj');
    path_table(v ~= max(v),:) = [];
    path_table.z_adj = path_table.z_adj-min(path_table.z_adj)+1;
    z_adj = [path_table(:,1) path_table(:,end)];
end