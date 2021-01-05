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
if isequal(config.adjust_intensity,"true") && ~isempty(config.adj_params)
    fprintf('%s\t Applying some intensity adjustments \n',datetime('now'))
elseif isequal(config.adjust_intensity,"true") && isempty(config.adj_params)
    error('%s\t Intensity adjustments requested but could not locate adjustment parameters.',...
        datetime('now'))
elseif isequal(config.adjust_intensity,"false")
    fprintf('%s\t Stitching without intensity adjustments \n',datetime('now'))
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
end

% Begin stitching
if ~isempty(p)
    parfor i = z_range
        % Print image being stitched
        fprintf('%s\t Stitching image %d \n',datetime('now'),i);
        stitch_worker_loaded(img_name_grid(:,:,:,i),h_stitch_tforms(:,i),...
            v_stitch_tforms(:,i),config,i);
    end
else
    for i = z_range
        % Print image being stitched
        fprintf('%s\t Stitching image %d \n',datetime('now'),i);
        stitch_worker_loaded(img_name_grid(:,:,:,i),h_stitch_tforms(:,i),...
            v_stitch_tforms(:,i),config,i);
    end
end
end

function stitch_worker_loaded(img_grid,h_tforms,v_tforms, config, z_idx)

% Border padding amount along the edge to compensate for empty space left
% after alignemnt
border_pad = config.border_pad;

%Image grid info
[nrows,ncols,nchannels] = size(img_grid);

%Read images, adjust intensities, apply translations for multichannel
A = cell(size(img_grid));
A{1} = imread(img_grid{1});
[img_height,img_width,~] = size(A{1});
ref_fixed = imref2d([img_height img_width]);
for k = 1:nchannels 
    for i = 1:nrows
        for j = 1:ncols
            c_idx = config.stitch_sub_channel(k);
            
            % Read image
            A{i,j,k} = imread(img_grid{i,j,k});
                    
            % Crop or pad images
            if ~isequal(size(A{i,j,k}),[img_height,img_width])
                A{i,j,k} = crop_to_ref(A{1},A{i,j,k});
            end
            
            % Apply intensity adjustments
            if isequal(config.adjust_intensity,"true")
                adj_params.r = i;
                adj_params.c = j;
                adj_params.idx = c_idx;
                
               % Crop laser width adjustment if necessary
                if length(config.adj_params.y_adj{k}) ~= length(config.adj_params.y_adj{1})
                    config.adj_params.y_adj{k} = crop_to_ref(config.adj_params.y_adj{1},...
                        config.adj_params.y_adj{k});
                    A{i,j,k} = crop_to_ref(A{1,1,1},A{i,j,k});
                end
                A{i,j,k} = apply_intensity_adjustment(A{i,j,k},config.adj_params,...
                    'r',i,'c',j,'idx',c_idx);
            end

            %Apply alignment transforms
            if ~isempty(config.alignment_table) && k>1                
                alignment_table_sub = config.alignment_table{i,j}...
                    (config.alignment_table{i,j}.Reference_Z == z_idx,:);
                channel_idx = nchannels+1+(k-2)*3;
                x = table2array(alignment_table_sub(:,channel_idx));
                y = table2array(alignment_table_sub(:,channel_idx+1));
                
                tform = affine2d([1 0 0; 0 1 0; x y 1]);
                A{i,j,k} = imwarp(A{i,j,k}, tform,'OutputView',ref_fixed,'FillValues',0); 
            end
        end
    end
end

%Convert images to single
A = cellfun(@(s) single(s),A,'UniformOutput',false);

%Calculate overlaps in pixels
v_overlap = round(img_height*config.overlap);
overlap_v_min = 1:v_overlap;
h_overlap = round(img_width*config.overlap);
overlap_h_min = 1:h_overlap;

%Sizes of optimal, fully stitched image
full_width = img_width*ncols-h_overlap*(ncols-1);
full_height = img_height*nrows-v_overlap*(nrows-1);

%Generate pixel merge weights using sigmoid function
%Horizontal
w_h = linspace(-config.sd,config.sd,h_overlap-border_pad*2);
w_h = 1./(1 + exp(-(w_h)));
w_h = (w_h - min(w_h))/(max(w_h)-min(w_h));
w_h = cat(2,zeros(1,border_pad),w_h,ones(1,border_pad));

%Save first column before horizontal stitching
B = A(:,1,:);

% Apply Horizontal Translations
a = 1;
for i = 1:nrows
for j = 1:ncols-1
    % Update overlap region of the left image 
    overlap_h_max = size(B{i,1},2)-h_overlap+1:size(B{i,1},2);
    
    % Load stitch parameters
    final_tform = affine2d([1 0 0; 0 1 0; h_tforms(a) h_tforms(a+1) 1]);

    % Adjust horizontally overlapped pixels based on translations
    ref_fixed2 = imref2d([img_height img_width+ceil(final_tform.T(3))]);

    % Transform and merge images (faster to for loop on each channel)
    for k = 1:nchannels
        reg_img = imwarp(A{i,j+1,k},final_tform,'OutputView',ref_fixed2,'FillValues',0);
        
        %Adjust intensity again?
        %adj_factor = prctile(B{i,k}(:,overlap_h_max),75)/prctile(reg_img(:,overlap_h_min),75);
        %reg_img = reg_img * adj_factor;
        
        B{i,k} = blend_images(reg_img,B{i,k},w_h,overlap_h_min,overlap_h_max,...
            config.blending_method(k),'horizontal');
    end
    a = a+2;
end
end

%Vertical
w_v = linspace(-config.sd,config.sd,v_overlap-border_pad*2);
w_v = 1./(1 + exp(-(w_v)))';
w_v = (w_v - min(w_v))/(max(w_v)-min(w_v));
w_v = cat(1,zeros(border_pad,1),w_v,ones(border_pad,1));

%Crop horizontally stitched images to minimum width
min_width = min(cellfun(@(s) size(s,2),B(:,1)));
B = cellfun(@(s) s(:,1:min_width), B,'UniformOutput',false);
I = B(1,:);

% Apply Vertical Translations
b = 1;
for i = 1:length(B)-1
    % Update overlap region of the top image 
    overlap_v_max = size(I{1},1)-v_overlap+1:size(I{1},1);
    
    % Load stitch parameters
    final_tform = affine2d([1 0 0; 0 1 0; v_tforms(b) v_tforms(b+1) 1]);
    
    % Adjust vertically overlapped pixels based on translations
    ref_fixed2 = imref2d([img_height+ceil(final_tform.T(6)) size(I{1},2)]);
        
    % Transform image
    for k = 1:nchannels
        reg_img = imwarp(B{i+1,k},final_tform,'OutputView',ref_fixed2,'FillValues',0);
                        
        %Adjust intensity again?
        adj_factor = prctile(I{k}(overlap_v_max,:),75)/prctile(reg_img(overlap_v_min,:),75);
        reg_img = reg_img * adj_factor;
                        
        I{k} = blend_images(reg_img,I{k},w_v,overlap_v_min,overlap_v_max,...
            config.blending_method(k),'vertical'); 
    end
    b = b+2;
end

% Postprocess the image with various filters, background subtraction, etc.
for i = 1:nchannels
    c_idx = config.stitch_sub_channel(i);
    I{c_idx} = postprocess_image(config, I{c_idx}, c_idx);
end

%Crop or pad images based on ideal size
I = cellfun(@(s) crop_to_ref(zeros(full_height,full_width),s),I,'UniformOutput',false);

%Save images as individual channels (will be large)
for i = 1:nchannels
    img_name = sprintf('%s_%s_C%d_%s_stitched.tif',config.sample_name,num2str(z_idx,'%04.f'),i,config.markers(i));
    img_path = fullfile(char(config.output_directory),'stitched',img_name);
    imwrite(uint16(I{i}),img_path)
end

end

function ref_img = blend_images(mov_img,ref_img,w,overlap_min,overlap_max,...
    blending_method,direction)


switch blending_method
    case 'sigmoid'
        if isequal(direction,'horizontal')
            %Perform non-linear weight
            inv_w = abs(1-w);
            ref_img(:,overlap_max) = ref_img(:,overlap_max).*inv_w +...
                mov_img(:,overlap_min).*w;
            mov_img(:,overlap_min) = [];

            %Concatanate images
            ref_img = horzcat(ref_img,mov_img);
        else
            %Perform non-linear weight
            inv_w = abs(1-w);
            ref_img(overlap_max,:) = ref_img(overlap_max,:).*inv_w +...
                mov_img(overlap_min,:).*w;
            mov_img(overlap_min,:) = [];
        
            %Concatanate images
            ref_img = vertcat(ref_img,mov_img);
        end
    case 'max'
        if isequal(direction,'horizontal')
            %Take max in the overlapping region
            ref_img(:,overlap_max) = max(ref_img(:,overlap_max),mov_img(:,overlap_min));
            mov_img(:,overlap_min) = [];
            
            %Concatanate images
            ref_img = horzcat(ref_img,mov_img);
            
        else
            %Take max in the overlapping region
            ref_img(overlap_max,:) = max(ref_img(overlap_max,:),mov_img(overlap_min,:));
            mov_img(overlap_min,:) = [];
                     
            %Concatanate images
            ref_img = vertcat(ref_img,mov_img);
        end    
end   
end