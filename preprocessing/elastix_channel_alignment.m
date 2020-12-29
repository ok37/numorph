function out = elastix_channel_alignment(config, path_table, warp_images, chunk_pad)
%--------------------------------------------------------------------------
% Aligning multiple channels to a reference nuclei channel using a coarse
% non-linear B-Spline registration. This uses a modified version of
% melastix, a MATLAB wrapper for the elastix registration toolbox. Melastix
% will need to be in your MATLAB path and Elastix libraries will need to be
% properly setup. In its current implementation, this function contains many
% parameters that were empirically adjusted to account for subtle
% deformations. Further fine-tuning is likley necessary if the expected
% deformations are large. Also note, this has relatively high memory
% requirements if the image stacks are large. Outputs are 2D series of
% aligned channel images + reference channel images in the 'aligned'
% directory
%--------------------------------------------------------------------------

if nargin<4 && ~isfield(config,'chunk_pad')
    chunk_pad = 30; %Important parameter. If misalignment is very large, you will have to increase to prevent zeroed areas from occuring in the output chunk edges
else
    chunk_pad = config.chunk_pad;
end

% Parameters not found in NMp_template. Defaults
img_gamma_adj = [1,1,1]; % Apply gamma during intensity adjustment. Decrease for channels with low background
smooth_blobs = "false"; % Smooth blobs using gaussian filter

% Unpack image information variables
sample_name = config.sample_name;
markers = config.markers;
output_directory = config.output_directory;
home_path = config.home_path;

% Registration-specific parameters
max_chunk_size = config.max_chunk_size; 
mask_int_threshold = config.mask_int_threshold;
s = config.resample_s;
histogram_bins = config.hist_match; % Match histogram bins to reference channel? If so, specify number of bins. Otherwise leave at 0. This can be useful for low contrast images

if isequal(config.pre_align,"true") && isfield(config,'alignment_table')
    alignment_table = config.alignment_table;
else
    alignment_table = [];
end

% Unpack intensity adjustment variables
lowerThresh = config.lowerThresh;
upperThresh = config.upperThresh;

% Check if equal resolution
equal_res = cellfun(@(s) all(config.resolution{1}(1:2) == s(1:2)),config.resolution);
if ~all(equal_res)
    res_adj = cellfun(@(s) config.resolution{1}(1:2)./s(1:2),config.resolution,'UniformOutput',false);
else
    res_adj = repmat({ones(1,2)},1,length(markers));
end
    
% Get x,y tile positions from path_table
x = unique(path_table.x);
y = unique(path_table.y);
assert(length(x) == 1, "More than 1 tile column detected")
assert(length(y) == 1, "More than 1 tile row detected")

% Paths to elastix parameter files
parameter_path = fullfile(home_path,'supplementary_data','elastix_parameter_files','channel_alignment');
transform_path = cell(length(markers)-1,3); outputDir = cell(1,length(markers)-1);
for i = 1:length(markers)-1
    parameter_dir = dir(parameter_path);
    idx = find({parameter_dir.name} == config.param_folder(i), 1);
    if ~isempty(idx)
        % Detect which transform is in the parameter file
        parameterSub = dir(fullfile(parameter_path,config.param_folder(i)));
        parameterSub = parameterSub(arrayfun(@(s) endsWith(s.name,'.txt'),parameterSub));
        for j = 1:length(parameterSub)
            file_path = fullfile(parameterSub(1).folder,parameterSub(j).name);
            text = textread(file_path,'%s','delimiter','\n');
            n = find(cellfun(@(s) contains(s,'(Transform '),text));
            if contains(text(n),'Translation')
                transform_path{i,1} = file_path;
            elseif contains(text(n),'Affine') || contains(text(n),'Euler')
                transform_path{i,2} = file_path;
            elseif contains(text(n),'BSpline')
                transform_path{i,3} = file_path;
            end
        end
    else
        error("Could not locate elastix parameter folder %s "+...
            "in %s",config.param_folder(i),parameter_path)
    end
    % Create empty tmp directory for saving images/transforms
    outputDir{i} = string(fullfile(config.output_directory,sprintf('tmp%d',i)));
    if ~exist(outputDir{i},'dir')
        mkdir(outputDir{i})
    else
        rmdir(outputDir{i},'s')
        mkdir(outputDir{i})
    end
end

%% Check if previous registration parameters exist
if isequal(config.load_alignment_params,"true") || isequal(config.load_alignment_params,"update")
   fprintf("%s\t Loading previous alignment parameters \n",datetime('now'))
   varDir = fullfile(output_directory,'variables','alignment_params.mat');   
   m = matfile(varDir);
    if isempty(m.alignment_params(x,y))
        warning("Could not locate elastix parameters in variables folder. Running alignment from scratch ")
        using_loaded_parameters = 'false';
    else
        out = m.alignment_params(x,y);
        out = out{1};
        using_loaded_parameters = 'true';
    end
else
    using_loaded_parameters = 'false';
end

%% Subset chunk ranges or channels
if ~isempty(config.align_channels)
    c_idx = [1,ismember(2:length(markers),config.align_channels)];
    align_channels = config.align_channels;
else
    c_idx = ones(1,length(markers));
    align_channels = 2:length(markers);
end    

total_images = height(path_table)/length(markers);
if ~isempty(config.align_slices) && isequal(using_loaded_parameters,'true')
    % Subset by slice range   
    for i = 1:length(config.align_slices)
        chunk_start(i) = min(config.align_slices{i});
        chunk_end(i) =  max(config.align_slices{i});
        chunk_start_adj(i) =  max(1,chunk_start(i)-chunk_pad);
        chunk_end_adj(i) =  min(chunk_end(i)+chunk_pad,total_images);
        fprintf('Aligning slices %d to %d\n', chunk_start(i), chunk_end(i))
    end
    
    init_tform = cell(1,length(markers)-1);
    for i = 1:length(markers)-1
        if c_idx(i+1) == 0 
            continue 
        end
        init_tform{i} = out{1,i};
        init_tform{i}.TransformParameters = init_tform{i}.TransformParameters(1);
    end
    
    align_chunks = 1:length(chunk_start);
    z_min_adj = min(chunk_start_adj);
    z_max_adj = max(chunk_end_adj);
    path_table = path_table(ismember(path_table.z,z_min_adj:z_max_adj),:);
    z_range_save = unique([config.align_slices{:}]);
    
elseif  ~isempty(config.align_chunks) && isequal(using_loaded_parameters,'true')
    % Subset by chunk
    align_chunks = config.align_chunks;
    chunk_start =  cellfun(@(s) s.chunk_start,out(align_chunks,1));
    chunk_end =  cellfun(@(s) s.chunk_end,out(align_chunks,1));
    chunk_start_adj =  cellfun(@(s) s.chunk_start_adj,out(align_chunks,1));
    chunk_end_adj =  cellfun(@(s) s.chunk_end_adj,out(align_chunks,1));
    fprintf('Aligning chunks %d\n', align_chunks)

    z_min_adj = min(chunk_start_adj);
    z_max_adj = max(chunk_end_adj);
    path_table = path_table(ismember(path_table.z,z_min_adj:z_max_adj),:);
    
    slices = arrayfun(@(s,t) s:t, chunk_start, chunk_end, 'UniformOutput', false);
    z_range_save = unique([slices{:}]);
    
else
    % Use all images
    z_min_adj = 1;
    z_max_adj = total_images;
    align_chunks = NaN;
    z_range_save = 1:total_images;
    init_tform = cell(1,length(markers)-1);
end

%% Preconfigure matrixes
tic
% Get basic image and file information 
n_images = height(path_table)/length(markers);
z_range_adj = z_min_adj:z_max_adj;
tempI = imread(path_table.file{1});
[nrows, ncols] = size(tempI);
cropping_flag =  [0 0 0]; % Flag for cropping images to reference channel
dim_adj = round([nrows ncols total_images]./s); % Resampled image dimensions

% Check if channels have the same sized images
for i = [1,align_channels]
    path_sub = path_table(path_table.markers == markers(i),:);
    tempI2 = imread(path_sub.file{1});    
    if size(tempI) ~= size(tempI2)
        cropping_flag(i) = 1;
    end
end

% Initialize matrices for storing raw images
I_raw = cell(1,length(markers));
I_raw(:) = {zeros(nrows,ncols,total_images,'uint16')};

%Start parallel pool
try
    p = gcp('nocreate');
catch
    p = 1;
end
if isempty(p)
    parpool(length(config.markers))
end

%% Read the images
for i = 1:length(markers)
    % Read image slices for only the channels that will be aligned
    if i > 1 && c_idx(i) == 0
        continue
    end
    fprintf('Reading images and pre-processing marker %s \n', markers(i))
    
    % If using pre-aligned, account for a shift in z by reading images from
    % alignment table
    if ~isempty(alignment_table) && i>1
        path_sub = cell(total_images,1);
        path_sub(alignment_table.Reference_Z) = alignment_table{:,"file_" + num2str(i)};
    else
        path_sub = path_table(path_table.markers == markers(i),:).file;
    end
    
    for j = 1:length(z_range_adj)
        if ~isempty(path_sub{z_range_adj(j)})
            I_raw{i}(:,:,z_range_adj(j)) = read_slice(path_sub{z_range_adj(j)},...
                res_adj{i},cropping_flag(i));
        end
    end
end

% Apply intensity adjustments for each tile
if any(arrayfun(@(s) isequal(s,"true"),config.adjust_intensity))
    I_raw = apply_intensity_adjustments_tile(I_raw, config, markers, c_idx);
end

% Apply pre-alignments if necessary
if ~isempty(alignment_table)
    for i = 2:length(markers)
        fprintf('Pre-aligning marker %s \n', markers(i))
        path_sub = path_table(path_table.markers == markers(i),:);
        
        % Get respective translations for marker
        if ~isempty(alignment_table)
            files = alignment_table{:,"file_" + num2str(i)};
            x_shift = alignment_table{:,"X_Shift_" + markers(i)};
            y_shift = alignment_table{:,"Y_Shift_" + markers(i)};
        end
        for j = 1:length(z_range_adj)
            % Find the correct image
            idx = find(arrayfun(@(s) isequal(s,path_sub.file(z_range_adj(j))),files),1);
            if ~isempty(idx)
                % Apply translation
                I_raw{i}(:,:,z_range_adj(j)) = imtranslate(I_raw{i}(:,:,z_range_adj(j)),...
                    [x_shift(idx) y_shift(idx)]);
            end
        end
    end
end


if isequal(using_loaded_parameters,'false') || isequal(config.load_alignment_params,"update")
    I = cell(1,length(I_raw));
    for i = 1:length(markers)
        if i > 1 && c_idx(i) == 0
            I{i} = [];
            continue
        else
            I{i} = I_raw{i};
        end
        
        % Smooth blobs using Guassian filter to reduce effects of
        % noise. Not a big impact if resampling
        for j = 1:length(z_range_adj)
            if isequal(smooth_blobs,"true")
                I{i}(:,:,j) = imgaussfilt(I{i}(:,:,j),1);
            end
        end
        
        % Adjust intensity and convert to int16 to run in elastix
        I{i} = im2int16(imadjustn(I{i},[0 upperThresh(i)],[],img_gamma_adj(i)));
        % Resample image
        I{i} = imresize3(I{i},dim_adj);
    end
    
    % Histogram matching
    for i = align_channels
        if ~isempty(histogram_bins) && histogram_bins(i) > 0 && c_idx(i-1) > 0
            I{i} = imhistmatch(I{i},I{1},histogram_bins(i));
        end 
    end
    
    % Generate mask using reference channel
    if isempty(config.align_slices)
        mask = generate_sampling_mask(I{1},mask_int_threshold,lowerThresh,upperThresh);
    else
        % If specific slices selected, use all the features in the slice
        % objects
        mask = zeros(size(I{1}));
        mask(:,:,chunk_start_adj(1):chunk_end_adj(1)) = 1;
        fprintf('Using masked ROI for slices \n')
    end
    
    % Determine chunk positions
    if isempty(config.align_chunks) && isempty(config.align_slices)
        % If total number of images is less than the max chunk size, adjust
        % padding
        if max_chunk_size>total_images
           chunk_pad = 0; 
        end
        [chunk_start, chunk_end, chunk_start_adj, chunk_end_adj] = get_chunk_positions(mask, n_images, max_chunk_size, chunk_pad);
        align_chunks = 1:length(chunk_start);
    end
end
fprintf('Images loaded and pre-processed in %d seconds\n', round(toc))

%% Perform the registration
if isequal(using_loaded_parameters,'false') || isequal(config.load_alignment_params,"update")
    %Initial registration by translation on whole downsampled stack
    fprintf('Performing intial registration \n')
    parfor i = 1:length(markers)-1
       if c_idx(i+1) == 0 || isempty(transform_path{i,1}) || ~isempty(config.align_slices) 
            fprintf("Skipping intial whole stack registration for channel %d\n", c_idx(i+1))
            init_tform{i} = [];
            continue
        elseif i ~= 1
            pause(1)
       end
       [init_tform{i},~] = elastix(I{i+1}(:,:,z_range_adj),I{1}(:,:,z_range_adj),...
            outputDir{i},transform_path(i,1),'s',s,'threads',[]);
        fprintf("Intial transform parameters: %s\n",sprintf("%.3f\t",init_tform{i}.TransformParameters{1}.TransformParameters))
        % Preview registration
        % imshowpair((img(:,:,500)),(I{1}(:,:,500)))
    end
        
    % Save transform parameters into cell array
    if isequal(using_loaded_parameters, 'false')
        out = cell(length(align_chunks),length(markers)-1);
    end
    
    fprintf('Performing chunk-wise registration \n')
    for i = align_chunks
        fprintf('Aligning chunk %d out of %d\n', i, length(align_chunks))

        % Take the mask for the respective chunk to align
        mask_chunk = mask(:,:,chunk_start_adj(i):chunk_end_adj(i));
        
        % Clip the bottom and top sections in padded regions
        mask_chunk(:,:,1:chunk_pad-1) = 0;
        mask_chunk(:,:,size(mask_chunk,3)-chunk_pad-1:size(mask_chunk,3)) = 0;
        
        % Take chunk images from reference channel
        I_ref = I{1}(:,:,chunk_start_adj(i):chunk_end_adj(i));
        
        % Save chunk start and end location with and without accounting for
        % padding
        chunk_start1 = chunk_start_adj(i); chunk_end1 = chunk_end_adj(i);
        chunk_start2 = chunk_start(i); chunk_end2 = chunk_end(i);
        
        % Chain transform parameters
        transforms = [transform_path(:,2), transform_path(:,3)];
        transforms = transforms(cellfun(@(s) ~isempty(s),transforms));
        
        tic
        parfor j = 1:length(markers)-1
           if c_idx(j+1) == 0
                continue
            elseif j ~= 1
                pause(1)
           end
           
            if ~isempty(init_tform{j})
                % Take initial transform of whole stack and save it as an
                % elastix parameters file
                init_tform{j}.TransformParameters{1}.Size(3) = size(mask_chunk,3);
                init_tform_path = sprintf('%s/init_transform_C%d_%d_%d_%d_%s.txt',outputDir{j},j,y,x,i,config.sample_name);
                elastix_paramStruct2txt(init_tform_path,init_tform{j}.TransformParameters{1})

                % Perform registration on the chunk for this channel
                [out{i,j},~] = elastix(I{j+1}(:,:,chunk_start1:chunk_end1),I_ref,...
                    outputDir{j},transforms(j,:),'fMask', mask_chunk,...
                    't0',init_tform_path,'s',s,'threads',[]);

                % Save initial transform filename into structure
                out{i,j}.TransformParameters = horzcat(init_tform{j}.TransformParameters,...
                    out{i,j}.TransformParameters);
            else
                % Perform registration on the chunk for this channel
                [out{i,j},~] = elastix(I{j+1}(:,:,chunk_start1:chunk_end1),I_ref,...
                    outputDir{j},transforms(j,:),'fMask', mask_chunk,'s',s,'threads',[]);
            end
            
            % Save chunk positions
            out{i,j}.chunk_start_adj = chunk_start1;
            out{i,j}.chunk_end_adj = chunk_end1;
            out{i,j}.chunk_start = chunk_start2;
            out{i,j}.chunk_end = chunk_end2;
            out{i,j}.x = x;
            out{i,j}.y = y;

            % Reset temporary directory
            rmdir(outputDir{j},'s')
            mkdir(outputDir{j})
            
            %Preview registration?            
            %output = transformix(I{j+1}(:,:,chunk_start1:chunk_end1),out{i,j},s,[]);            
            %output = int16(output);
            %slices = 1:5:size(output,3);
            %I_sample1 = im2uint16(I_ref(:,:,slices));
            %I_sample2 = im2uint16(output(:,:,slices));
            %options.overwrite = true;
            %saveastiff(I_sample1,sprintf('topro_%d.tif',i),options)
            %saveastiff(I_sample2,sprintf('transform_%d.tif',i),options)
        end
    fprintf('Finished registration in %d seconds\n', round(toc))
    end
    % Cleanup temporary directories
    cellfun(@(s) rmdir(s,'s'), outputDir)
end

% Return if not aplying transformations
if ~warp_images || isequal(config.save_images,'false')
    return
end

%% Apply transformations to adjusted images
% Subset chunks that we're interested in aligning
out = out(align_chunks,:);

I_raw = apply_transformations(I_raw,out,total_images,ncols,nrows,markers,c_idx);

%% Save the images
% Create directory to store images
if exist(fullfile(output_directory,'aligned'),'dir') ~= 7
    mkdir(fullfile(output_directory,'aligned'));
end

if ~isempty(config.align_slices)
    z_range_save = unique([config.align_slices{:}]);
    fprintf("%s\t Saving only slices %d to %d \n",datetime('now'),min(z_range_save),max(z_range_save))
end

if ~isempty(config.align_slices) && isequal(config.load_alignment_params,'true')
    z_range_save = config.align_slices{1};
end

% Save aligned images
for j = 1:length(markers)
   if j > 1 && c_idx(j) == 0
        continue
   end
    fprintf("%s\t Writing aligned images for marker %s\n",datetime('now'),markers(j))
    for i = 1:length(z_range_save)
        img_name = sprintf('%s_%s_C%d_%s_0%d_0%d_aligned.tif',sample_name,num2str(z_range_save(i),'%04.f'),j,markers(j),y,x);
        img_path = fullfile(output_directory,'aligned',img_name);        
        %fprintf("%s\t Writing aligned image %s\n",datetime('now'),img_name)
        imwrite(I_raw{j}(:,:,z_range_save(i)),img_path)
    end
end

end


function img = read_slice(file, res_adj, cropping_flag)
% Read image and apply adjustments
img = imread(file);

% Resample image to image reference image if different
% resolutions
if isequal(res_adj,ones(1,2))
    img = imresize(img,round(size(img)./res_adj),'Method','bicubic');
end

% Crop or pad image if it is not the same size as reference image
if cropping_flag == 1
    img= crop_to_ref(tempI,img);
end

end


function I_raw = apply_transformations(I_raw,out,n_images,ncols,nrows,markers,c_idx)
% Apply elastix transform parameters

% Get chunk sizes and locations in the stack with/without padding
n_chunks = size(out,1);
chunk_start_adj = cellfun(@(s) s.chunk_start_adj,out);
chunk_end_adj = cellfun(@(s) s.chunk_end_adj,out);

chunk_start = cellfun(@(s) s.chunk_start,out);
chunk_end = cellfun(@(s) s.chunk_end,out);

% Apply transformatin for each chunk
for i = 1:n_chunks
    if i == n_chunks && chunk_end(i) == n_images
        a = 0;
    else
        a = 1;
    end
    tic    
    for j = 1:length(markers)-1
       if c_idx(j+1) == 0
            continue
        elseif j ~= 1
            pause(1)
        end
        
       % Get locations for current chunk
       chunk_start1 = chunk_start_adj(i,j);
       chunk_end1 = chunk_end_adj(i,j);

        top_trim = abs(chunk_start1-chunk_start(i,j));
        bottom_trim = abs(chunk_end1-chunk_end(i,j))+a;

        chunk_start2 = chunk_start(i,j);
        chunk_end2 = chunk_end(i,j)-a; 

        % out contains transform parameters
        out2 = out{i,j};
        
        % Adjust size back from rescaled resolution 
        size1 = out2.TransformParameters{end}.Size;
        size1 = [ncols nrows size1(3)];

        % Backward compatibility. Move init_transform field as 1st in 3
        % transform parameters
        if length(out2.TransformParameters) == 2
            out2.TransformParameters = horzcat({out2.InitTransformParameters},out2.TransformParameters);
        end
        
        % Note using int32. Weird pixel errors if using uint/double data
        % types. Adjust sizes and adjust default settings.
        for k = 1:length(out2.TransformParameters)
            out2.TransformParameters{k}.Size = size1;
            out2.TransformParameters{k}.Spacing = [1 1 1];
            out2.TransformParameters{k}.DefaultPixelValue = -32768;
            out2.TransformParameters{k}.InitialTransformParametersFileName = 'NoInitialTransform';
        end
        out2.TransformParameters{end}.Transform = 'RecursiveBSplineTransform';        

        I_chunk = im2int16(I_raw{j+1}(:,:,chunk_start1:chunk_end1));
        
        % Apply transform
        I_chunk = transformix(I_chunk,out2,[1 1 1],[]);
        
        %Trim ends
        z_range = 1+top_trim:size(I_chunk,3)-bottom_trim;
        I_chunk = im2uint16(int16(I_chunk(:,:,z_range)));
        
        % Replace in image stack
        I_raw{j+1}(:,:,(chunk_start2:chunk_end2)) = I_chunk;   
    end
    fprintf('Finished transformation in %d seconds\n', round(toc))
end

end


function mask = generate_sampling_mask(I,mask_int_threshold,lowerThresh,upperThresh)

% Generate intensity threshold if not provided
%if isempty(mask_int_threshold)
%    mask_int_threshold(1) = (lowerThresh(1)+upperThresh(1))/2;
%    mask_int_threshold(2) = otsuthresh(im2uint16(mask(:)))*0.9;
%    mask_int_threshold(3) = mask_int_threshold1+mask_int_threshold2/2;
%end

%mask_int_threshold = 0.05;
    
% Generate mask
% Downsample to 20% resolution
I2 = imresize3(im2uint16(I),0.20,'linear'); 

% Calculate mask intensity threshold
% These are empirically determined based on adjusted otsu thresholding and
% lowerThresh determined from examining all images
if isempty(mask_int_threshold)
    mask_int_threshold = min(graythresh(I2)*0.75,(lowerThresh(1)/upperThresh(1))*1.25);
end

% Binarize mask and fill holes
mask = imbinarize(I2,mask_int_threshold);
mask = imfill(mask,26,'holes');
    
% Keep only largest compnent. Disconnected components will give errors
labels = bwconncomp(mask);
sizes = cellfun(@(s) length(s), labels.PixelIdxList);
[~, idx] = max(sizes);
for i = 1:length(sizes)
    if i ~= idx
        mask(labels.PixelIdxList{i}) = 0;
    end
end
    
% Resize back to original size
mask = imresize3(double(mask), size(I),'method', 'nearest');
signal = sum(mask(:))/numel(mask);

fprintf('\n Using a mask intensity threshold of %.4f \n', mask_int_threshold)
fprintf('\n Using %.4f percent of pixels \n', signal*100)

end


function [chunk_start, chunk_end, chunk_start_adj, chunk_end_adj] = get_chunk_positions(mask, num_images, max_chunk_size, chunk_pad)
% Get chunk position from binary mask

ind = find(sum(sum(mask,1),2) > 0);
z_start = min(ind);
z_end = max(ind);
z_pos = z_start:z_end;

n_chunks = 1;
chunk_size = length(z_pos);

if chunk_size < max_chunk_size
    chunk_start = max(z_start - chunk_pad,1);
    chunk_end = min(z_end + chunk_pad,num_images);

    chunk_start_adj = chunk_start;
    chunk_end_adj = chunk_end;
else
    while chunk_size > max_chunk_size
        n_chunks = n_chunks + 1;
        chunk_size = ceil(length(z_pos)/n_chunks);
    end
    chunk_start = z_start:chunk_size:z_end;
    chunk_end = chunk_start + chunk_size;

    chunk_start(1) = z_start;
    chunk_end(end) = z_end;

    chunk_start_adj = chunk_start;
    chunk_start_adj(2:end) = chunk_start_adj(2:end) - chunk_pad;
    chunk_end_adj = chunk_end;
    chunk_end_adj(1:(n_chunks-1)) = chunk_end_adj(1:(n_chunks-1)) + chunk_pad;
end

end


function I_raw = apply_intensity_adjustments_tile(I_raw, config, markers, c_idx)
% Specific to elastix channel alignment
% Apply intensity adjustments prior to aligning images

adj_params = config.adj_params;

% Crop y_adj and flatfield to match reference 
for i = 1:length(markers)
   if c_idx(i) == 0
        continue
   else
      fprintf('Applying intensity adjustments for marker %s\n',markers(i));
   end
    adj_params.y_adj{i} = crop_to_ref(adj_params.y_adj{1},adj_params.y_adj{i});
    adj_params.flatfield{i} = crop_to_ref(adj_params.flatfield{1},adj_params.flatfield{i});
    adj_params.darkfield{i} = crop_to_ref(adj_params.darkfield{1},adj_params.darkfield{i});

    % Adjust intensities
    for j = 1:size(I_raw{i},3)
       if sum(I_raw{i}(:,:,j) == 0)
           continue
       end
       if isequal(config.adjust_tile_shading(c_idx(i)),'basic')
           I_raw{i}(:,:,j) = apply_intensity_adjustment(I_raw{i}(:,:,j),...
               'flatfield', adj_params.flatfield{i},...
               'darkfield', adj_params.darkfield{i});
       elseif isequal(config.adjust_tile_shading(c_idx(i)),'manual')
          I_raw{i}(:,:,j) = apply_intensity_adjustment(I_raw{i}(:,:,j),...
              'y_adj',adj_params.y_adj{i});
       end
    end
end

end