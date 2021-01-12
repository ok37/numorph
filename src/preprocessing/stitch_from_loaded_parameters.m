function stitch_from_loaded_parameters(path_table, h_stitch_tforms, v_stitch_tforms, config)
% Stitch using previously calculated translations

% Stitch images using previously calculated parameters
fprintf('%s\t Begin stitching \n',datetime('now'))

% Create directory for stitched images
if ~exist(fullfile(config.output_directory,'stitched'),'dir')
  mkdir(fullfile(config.output_directory,'stitched'))
end

% Check if only certain channels to be stitched
if isempty(config.stitch_sub_channel)
    config.stitch_sub_channel = 1:length(config.markers);
end

% Generate image name grid
img_name_grid = cell(max(path_table.y),max(path_table.x),...
    length(config.stitch_sub_channel), max(path_table.z_adj));
path_table = sortrows(path_table,["z_adj","channel_num","x","y"],'ascend');

% Arrange images into correct positions
try
    img_name_grid = reshape(path_table.file,size(img_name_grid));
catch ME
    disp(ME.message)
    if isequal(ME.identifier,'MATLAB:getReshapeDims:notSameNumel')
        error('%s\t Inconsistent image file information.Recalculate adjusted z and/or check configuration \n',char(datetime('now')))
    end
end

% Check if only certain sub-section is to be stitched
if ~isempty(config.stitch_sub_stack)
    z_range = config.stitch_sub_stack;
else
    z_range = 1:size(img_name_grid,4);
end

% Check if intensity adjustments are specified
if isequal(config.adjust_intensity,"true")
    if ~isfield(config,'adj_params')
        error('%s\t Intensity adjustments requested but could not locate adjustment parameters.',...
            datetime('now'))
    end
    adj_params = cell(1,length(config.markers));
    for i = 1:length(config.markers)
        adj_params{i} = config.adj_params.(config.markers(i));
    end
    fprintf('%s\t Applying some intensity adjustments \n',datetime('now'))
else
    fprintf('%s\t Stitching without intensity adjustments \n',datetime('now'))
    adj_params = [];
end

% Check for alignment table
if isequal(config.load_alignment_params,"true") && ~isempty(config.alignment_table)
    fprintf('%s\t Applying channel alignment parameters during stitching \n',datetime('now'))
end

%Start parallel pool
try
    p = gcp('nocreate');
catch
    warning("Could not load Parallel Computing Toolbox. It's recommended "+...
        "that this toolbox is installed to speed up stitching and subsequent "+...
        "analysis.")
    p = 1;
end

if isempty(p) && length(z_range)>1
    parpool
elseif isfield(p,'NumWorkers') && p.NumWorkers == 2
    delete(p)
    parpool
end

% Begin stitching
if ~isempty(p)
    parfor i = z_range
        % Print image being stitched
        fprintf('%s\t Stitching image %d \n',datetime('now'),i);
        stitch_worker_loaded(img_name_grid(:,:,:,i),h_stitch_tforms(:,i),...
            v_stitch_tforms(:,i),config,i,adj_params);
    end
else
    for i = z_range
        % Print image being stitched
        fprintf('%s\t Stitching image %d \n',datetime('now'),i);
        stitch_worker_loaded(img_name_grid(:,:,:,i),h_stitch_tforms(:,i),...
            v_stitch_tforms(:,i),config,i,adj_params);
    end
end
end

function stitch_worker_loaded(img_grid,h_tforms,v_tforms, config, z_idx, adj_params)
% Stitch image grid

% Border padding amount along the edge to compensate for empty space left
% after alignemnt
border_pad = config.border_pad;

% Image grid info
[nrows,ncols,nchannels] = size(img_grid);

% Read images, adjust intensities, apply translations for multichannel
A = read_stitching_grid(img_grid,config.stitch_sub_channel,config.markers,...
    adj_params,config.alignment_table);
[img_height,img_width] = size(A{1});

% Calculate overlaps in pixels
h_overlap = ceil(img_width*config.overlap);
v_overlap = ceil(img_height*config.overlap);
x0 = ceil(h_overlap/2);
y0 = ceil(v_overlap/2);

overlap_h_min = 1:h_overlap;
overlap_v_min = 1:v_overlap;

% Sizes of optimal, fully stitched image
full_width = img_width*ncols-h_overlap*(ncols-1);
full_height = img_height*nrows-v_overlap*(nrows-1);

%Save first column before horizontal stitching
B = A(:,1,:);

% Apply Horizontal Translations
a = 1;
for i = 1:nrows
    if ncols == 1
        continue
    end
    for j = 1:ncols-1
        % Load stitch parameters
        final_tform = affine2d([1 0 0; 0 1 0; h_tforms(a) h_tforms(a+1) 1]);

        % Adjust horizontally overlapped pixels based on translations
        ref_fixed = imref2d([img_height img_width+floor(final_tform.T(3))]);
        
        % Create adjusted blending weights
        x1 = x0 - final_tform.T(3);
        if any(arrayfun(@(s) isequal(s,"sigmoid"),config.blending_method))
            w_h_adj = 1-1./(1 + exp(config.sd*(overlap_h_min-x1)));
        else
            w_h_adj = overlap_h_min/h_overlap;
        end
        
        % Clip ends based on horizontal translation where image intensity
        % is 0
        if final_tform.T(3) > 0
            adj_left = max(border_pad,ceil(final_tform.T(3)));
            w_h_adj(1:adj_left) = 0;
        else
            adj_right = max(border_pad,ceil(abs(final_tform.T(3))));
            w_h_adj(1:adj_right) = 0;
        end
        
        % Rescale non-cropped areas
        min_w_h_adj = min(w_h_adj(w_h_adj>0));
        w_h_adj(w_h_adj>=min_w_h_adj) = (w_h_adj(w_h_adj>=(min_w_h_adj)) - min_w_h_adj)./(1-min_w_h_adj);

        % Transform and merge images (faster to for loop on each channel)
        for k = 1:nchannels
            reg_img = imwarp(A{i,j+1,k},final_tform,'OutputView',ref_fixed,'FillValues',0);

            %Adjust intensity again?
            if isequal(config.adjust_tile_position,"true")
                overlap_h_max = size(B{i,1},2)-length(overlap_h_min)+1:size(B{i,1},2);

                adj_factor = prctile(B{i,k}(:,overlap_h_max),75,'all')/prctile(reg_img(:,overlap_h_min),75,'all');
                reg_img = reg_img * adj_factor;
            end

            B{i,k} = blend_images(reg_img,B{i,k},false,config.blending_method(k),w_h_adj);  
        end
        a = a+2;
    end
end

%Crop horizontally stitched images to minimum width
min_width = min(cellfun(@(s) size(s,2),B(:,1)));
B = cellfun(@(s) s(:,1:min_width), B,'UniformOutput',false);
I = B(1,:);

% Apply Vertical Translations
b = 1;
for i = 1:length(B)-1
    % Load stitch parameters
    final_tform = affine2d([1 0 0; 0 1 0; v_tforms(b) v_tforms(b+1) 1]);
    
    % Adjust vertically overlapped pixels based on translations
    ref_fixed = imref2d([img_height+floor(final_tform.T(6)) size(I{1},2)]);
    
    % Create adjusted blending weights
    y1 = y0 - final_tform.T(6);
    if isequal(config.blending_method(k),"sigmoid")
        w_v_adj = 1-1./(1 + exp(config.sd*(overlap_v_min-y1)))';
    else
        w_v_adj = (overlap_v_min/v_overlap)';
    end
    
    % Clip ends based on vertical translation where image intensity
    % is 0
    if final_tform.T(6) > 0
        adj_top = max(border_pad,ceil(final_tform.T(6)));
        w_v_adj(1:adj_top) = 0;
    else
        adj_bottom = max(border_pad,ceil(abs(final_tform.T(6))));
        w_v_adj(1:adj_bottom) = 0;
    end
    
    % Rescale non-cropped areas
    min_w_v_adj = min(w_v_adj(w_v_adj>0));
    w_v_adj(w_v_adj>=min_w_v_adj) = (w_v_adj(w_v_adj>=(min_w_v_adj)) - min_w_v_adj)./(1-min_w_v_adj);
        
    for k = 1:nchannels
        reg_img = imwarp(B{i+1,k},final_tform,'OutputView',ref_fixed,'FillValues',0,'SmoothEdges',true);
        %Adjust intensity again?
        if isequal(config.adjust_tile_position,"true")
            overlap_v_max = size(I{1},1)-v_overlap+1:size(I{1},1);

            adj_factor = prctile(I{k}(overlap_v_max,:),75,'all')/prctile(reg_img(overlap_v_min,:),75,'all');
            reg_img = reg_img * adj_factor;
        end
        I{k} = blend_images(reg_img,I{k},false,config.blending_method(k),w_v_adj); 
    end
    b = b+2;
end

% Postprocess the image with various filters, background subtraction, etc.
c_idx = config.stitch_sub_channel;
for i = 1:length(c_idx)
    I{i} = postprocess_image(config, I{i}, c_idx(i));
end

%Crop or pad images based on ideal size
I = cellfun(@(s) crop_to_ref(zeros(full_height,full_width),s),I,'UniformOutput',false);

%Save images as individual channels (will be large)
for i = 1:length(c_idx)
    img_name = sprintf('%s_%s_C%d_%s_stitched.tif',...
        config.sample_name,num2str(z_idx,'%04.f'),c_idx(i),config.markers(c_idx(i)));
    img_path = fullfile(char(config.output_directory),'stitched',img_name);
    imwrite(uint16(I{i}),img_path)
end

end
