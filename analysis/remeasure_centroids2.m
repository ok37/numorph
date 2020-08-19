function [centroids_new, centroids_back] = remeasure_centroids2(centroids,path_table_stitched,config)
%--------------------------------------------------------------------------
% Remeasure centroid intensities from centroid list.
%--------------------------------------------------------------------------
res = config.resolution; % Voxel resolution
method = 'max';        % How to measure centroids. mean or max
s = 1;  % Sampling neighbors. If set to 1, this will sample the centroid + 8 neighbors
sub_back = 'true';  % Option to subtract background
measure_markers = [2,3];  % Which markers to measure

markers = config.markers(measure_markers);

tempI = loadtiff(path_table_stitched.file{1});
[nrows, ncols] = size(tempI);

% +1 to go from python indexing
z_pos = unique(centroids(:,3))+1;

new_vals1 = zeros(length(markers),size(centroids,1));
new_vals2 = zeros(length(markers),size(centroids,1));


for i = 600%1:length(z_pos)
    fprintf('%s\t Remeasuring image %d \n',datetime('now'),z_pos(i));
    
    new_vals = zeros(length(markers),size(centroids,1));
    new_back_vals = zeros(length(markers),size(centroids,1));
    
    idx = centroids(:,3)==z_pos(i);
    
    % Get linear indexes for cooridnates
    y = centroids(idx,1)+1;
    x = centroids(idx,2)+1;
    idxs = get_mesh_indexes(y,x,s,nrows,ncols);
    
    % Read image and measure intensity
    for j = 1:length(markers)
        img_file = path_table_stitched(path_table_stitched.z == z_pos(i) &...
            path_table_stitched.markers == markers(j),:);  
    
        I = loadtiff(img_file.file{1});
        I = imadjust(I);

    
        % Get values at indexes
        if isequal(method,'mean')
            vals = cellfun(@(k) mean(I(k)),idxs);
        else
            vals = cellfun(@(k) max(I(k)),idxs);
        end
        
        % Save values
        new_vals(j,idx) = vals';
        
        % Optional: measure background
        if isequal(sub_back,'true')
            %[~, B] = smooth_background_subtraction(I,'false',80);
            %B = single(imgaussfilt(I,30));
                        
            idxs = get_mesh_indexes(y,x,30,nrows,ncols);
            B = single(I);
            mean_back = cellfun(@(k) mean(B(k)),idxs);
            std_back = cellfun(@(k) std(B(k)),idxs);
            
            vals_back = (single(vals)-mean_back)./std_back;
            
            % Save values
            new_back_vals(j,idx) = vals_back';
        end
    end
    
    % Add to final matrix
    new_vals1 = new_vals1 + new_vals;
    if isequal(sub_back,'true')
        new_vals2 = new_vals2 + new_back_vals;
    end
end

% Update measurements
centroids_new = centroids;
col_idx = 4+measure_markers;
centroids_new(:,col_idx) = new_vals1';

if isequal(sub_back,'true')
    centroids_back = new_vals2';
    save_path = fullfile(config.output_directory, 'background');
    if exist(save_path,'dir') ~= 7
        mkdir(save_path);
    end
    writematrix(centroids_back,fullfile(save_path,'back_values.csv'))
else
    centroids_back = [];
end

end


function idxs = get_mesh_indexes(y,x,s,nrows,ncols)

y1 = arrayfun(@(a,b) meshgrid(a-s:a+s,b-s:b+s),y,x,'UniformOutput',false);
x1 = arrayfun(@(a,b) meshgrid(a-s:a+s,b-s:b+s)',x,y,'UniformOutput',false);

% Remove coordinates outside of bounds of images
for j = 1:length(x1)
rm_idx = y1{j} < 1 | y1{j} > nrows;
rm_idx = rm_idx |  x1{j} < 1 | x1{j} > ncols;
y1{j} = y1{j}(~rm_idx);
x1{j} = x1{j}(~rm_idx);
end

% Linearize
idxs = cellfun(@(a,b) sub2ind([nrows,ncols],a,b),y1,x1,...
    'UniformOutput',false);

end