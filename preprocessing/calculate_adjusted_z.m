function z_adj = calculate_adjusted_z(path_table,nrows,ncols,markers,overlap,z_positions,z_window,lowerThresh,output_directory)
%Calculate adjusted z positions for iterative 2D image stitching based on
%cross-correlation of overlapping horizontal and vertial regions in
%adjacent tiles

%Initialize Stitching Displacement Matrices
h_disp_matrix = zeros(nrows, ncols);
v_disp_matrix = zeros(nrows, ncols);
q_h_disp_matrix = zeros(size(h_disp_matrix));
q_v_disp_matrix = zeros(size(v_disp_matrix));
h_low_flag = zeros(size(h_disp_matrix));
v_low_flag = zeros(size(v_disp_matrix));

tiles_x = unique(path_table.x)';
tiles_y = unique(path_table.y)';

% Perform overlap measurements only in first (reference) channel
path_table_ref = path_table(path_table.markers == markers(1),:);

%%% Find best pairwise positions of overlapping regions for horizontal tiles
fprintf(strcat(char(datetime('now')),'\t Performing horizontal pairwise z alignment \n'))
direction = 1;
tform_h_flag = zeros([z_positions,1]); tform_h_x = tform_h_flag; tform_h_y = tform_h_flag;
for i = 1:nrows
    for j = 2:ncols
        fprintf(strcat(char(datetime('now')),'\t Aligning tiles 0%dx0%d \n'),[i,j])
        path_ref = path_table_ref(path_table_ref.x==tiles_x(j-1) & path_table_ref.y==tiles_y(i),:);
        path_mov = path_table_ref(path_table_ref.x==tiles_x(j) & path_table_ref.y==tiles_y(i),:);
        [z_displacement, q, low_flag] = zAlign_stitch(path_mov,path_ref,overlap,z_positions,z_window,direction,lowerThresh(1));

        %Store results
        h_disp_matrix(i,j) = z_displacement;
        q_h_disp_matrix(i,j) = q;
        h_low_flag(i,j) = low_flag;
    end
end

if any(h_low_flag(:)==1)
    sub_displacement = mode(h_disp_matrix(:));
    h_disp_matrix(isnan(h_disp_matrix)) = sub_displacement;
end

fprintf(strcat(char(datetime('now')),'\t Final horizontal displacement matrix:\t',mat2str(h_disp_matrix),'\n'))

%%% Find best pairwise positions of overlapping regions for vertical tiles
fprintf(strcat(char(datetime('now')),'\t Performing vertical pairwise z alignment \n'))
direction = 2;
tform_v_flag = zeros([z_positions,1]); tform_v_x = tform_v_flag; tform_v_y = tform_v_flag;
for i = 1:ncols
    for j = 2:nrows
        fprintf(strcat(char(datetime('now')),'\t Aligning tiles 0%dx0%d \n'),[j,i])
        path_ref = path_table_ref(path_table_ref.x==tiles_x(i) & path_table_ref.y==tiles_y(j-1),:);
        path_mov = path_table_ref(path_table_ref.x==tiles_x(i) & path_table_ref.y==tiles_y(j),:);
        [z_displacement, q, low_flag] = zAlign_stitch(path_mov,path_ref,overlap,z_positions,z_window,direction,lowerThresh(1));

        %Store results
        v_disp_matrix(j,i) = z_displacement;
        q_v_disp_matrix(j,i) = q;
        v_low_flag(j,i) = low_flag;
    end
end

if any(v_low_flag(:)==1)
    sub_displacement = mode(v_disp_matrix(:));
    v_disp_matrix(isnan(v_disp_matrix)) = sub_displacement;
end
    fprintf(strcat(char(datetime('now')),'\t Final vertical displacement matrix:\t',mat2str(v_disp_matrix),'\n'))
%[0 -2 -2 -1;0 -1 -1 -1;0 -2 -2 -1;0 -2 -2 -2]
    %%% Use minimum spanning tree to get final z displacement matrix
    z_disp_matrix = min_span_tree_2(v_disp_matrix, h_disp_matrix, q_v_disp_matrix, q_h_disp_matrix);
    fprintf(strcat(char(datetime('now')),'\t Final displacement matrix:\t',mat2str(z_disp_matrix),'\n'))
    save(fullfile(char(output_directory),'variables','z_disp_matrix'),'z_disp_matrix')
    
    z_adj = apply_adjusted_z(path_table, z_disp_matrix);
    save(fullfile(char(output_directory),'variables','adjusted_z.mat'), 'z_adj') 
end
