function A = read_stitching_grid(img_grid,stitch_channels,markers,adj_params,alignment_table)
%--------------------------------------------------------------------------
% Read cell grid of images and apply adjustments prior to stitching. Return
% as single.
%--------------------------------------------------------------------------

[nrows, ncols,nchannels] = size(img_grid);
A = cell(nrows,ncols,nchannels);
A{1} = imread(img_grid{1});
[img_height,img_width] = size(A{1});
ref_fixed = imref2d([img_height img_width]);

% Read grid
for i = 1:nrows
    for j = 1:ncols
        for k = 1:length(stitch_channels)
            c_idx = stitch_channels(k);            
            
            % Read image
            A{i,j,k} = imread(img_grid{i,j,k});
                    
            % Crop or pad images
            if ~isequal(size(A{i,j,k}),[img_height,img_width])
                A{i,j,k} = crop_to_ref(A{1},A{i,j,k});
            end
            
            % Apply intensity adjustments
            if ~isempty(adj_params)
               % Crop laser width adjustment if necessary
                if length(adj_params{c_idx}.y_adj) ~= length(adj_params{1}.y_adj)
                    adj_params{c_idx}.y_adj = crop_to_ref(adj_params{1}.y_adj,adj_params{c_idx}.y_adj);
                    A{i,j,k} = crop_to_ref(A{1,1,1},A{i,j,k});
                end
                % Apply adjustments
                A{i,j,k} = apply_intensity_adjustment(A{i,j,k},adj_params{c_idx},...
                    'r',i,'c',j);
            end

            % Apply alignment transforms
            if ~isempty(alignment_table) && c_idx>1                
                % Get row index of image
                r_idx = find(cellfun(@(s) isequal(s,img_grid{i,j,k}),alignment_table{i,j}{:,2}));
                
                if ~isempty(r_idx)
                    % Get x,y translations and shift
                    x = alignment_table{i,j}{r_idx,"X_Shift_"+markers(c_idx)};
                    y = alignment_table{i,j}{r_idx,"Y_Shift_"+markers(c_idx)};
                
                    tform = affine2d([1 0 0; 0 1 0; x y 1]);
                    A{i,j,k} = imwarp(A{i,j,k}, tform,'OutputView',ref_fixed,'FillValues',0); 
                end
            end
            
        end
    end
end
            
%Convert images to single
A = cellfun(@(s) single(s),A,'UniformOutput',false);
            
end