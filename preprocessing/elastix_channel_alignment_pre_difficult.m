function elastix_channel_alignment_pre_difficult(x,y,path_table,config)
% Aligning multiple channels to a reference nuclei channel using a coarse
% non-linear B-Spline registration. This uses a modified version of
% melastix, a MATLAB wrapper for the elastix registration toolbox. Melastix
% will need to be in your MATLAB path and Elastix libraries will need to be
% properly setup.In its current implementation, this function contains many
% parameters that were empirically adjusted to account for subtle
% deformations. Further fine-tuning is likley necessary if the expected
% deformations are large. Also note, this has relatively high memory
% requirements if the image stacks are large. Outputs are saved,
% multi-channel 2D tiff images

% Parameters not found in TCp_template 
max_chunk_size = 200; %Important parameters. Decreasing might improve precision but takes longer. Higher values usually gives a more reproducible results
chunk_pad = 30; %Important parameter. If misalignment is very large, you will have to increase to prevent zero areas from occuring in the output chunk edges
img_gamma_adj = [0.8 0.8 0.8]; %Important parameter (0-1). Decrease for channels with low background. Different than gamma applied to image during intensity adjustment
enhance_blobs = "false"; %Enhance blob structures for more robust alignment

% Unpack image information variables
sample_name = config.sample_name;
markers = config.markers;
output_directory = config.output_directory;
number_of_cores = config.number_of_cores;
home_path = config.home_path;

% Registration-specific parameters
h_bins = config.h_bins;
mask_int_threshold = config.mask_int_threshold;
s = config.resample_s;

% Unpack intensity adjustment variables
lowerThresh = config.lowerThresh;
upperThresh = config.upperThresh;

% Paths to elastix parameter files
% Chose transform parameters. Use 32 bin parameters by default
parameterDir = fullfile(home_path,'elastix_parameter_files');
transform_path_init = fullfile(parameterDir,'parameters_Translation_3D_init.txt');
transform_path1 = cell(1,length(markers)-1); transform_path2 = transform_path1;
outputDir = transform_path1;
for i = 1:length(markers)-1
    % Create empty tmp directory for saving images/transforms
    outputDir{i} = fullfile(parameterDir,sprintf('tmp%d',i));
    if ~exist(outputDir{i},'dir')
        mkdir(outputDir{i})
    else
        rmdir(outputDir{i},'s')
        mkdir(outputDir{i})
    end

   % Locations of parameter files
   if h_bins(i) == 16
       transform_path1{i} = fullfile(parameterDir,'parameters_Rigid_3D_16.txt');
       transform_path2{i} = fullfile(parameterDir,'parameters_BSpline_3D_16.txt');
   else
       transform_path1{i} = fullfile(parameterDir,'parameters_Rigid_3D_32.txt');
       transform_path2{i} = fullfile(parameterDir,'parameters_BSpline_3D_32.txt');
   end
end

%s Downsampling factor for downsizing images [y,x,z];
histogram_bins = [0 0 128]; %Match histogram bins to reference channel? If so, specify number of bins. Otherwise leave at 0. This can be useful for low contrast images

%% Check if previous registration parameters exist
if isequal(config.load_alignment_params,"true")
    fprintf('Loading previous channel registraton parameters \n')
    if ~isequal(exist(fullfile(output_directory,'variables',sprintf('out_0%d_0%d.mat',y,x)),'file'),2)
        error('%s\t Could not locate elastix parameters in variables folder \n',datetime('now'))
    else
        load(fullfile(output_directory,'variables',sprintf('out_0%d_0%d.mat',y,x)), 'out')
    end
    using_loaded_parameters = 'true';
else
    using_loaded_parameters = 'false';
end

%% Read images
tic
num_images = height(path_table)/length(markers);
I_raw = cell(1,length(markers));

tempI = imread(path_table.file{1});
[nrows, ncols] = size(tempI);
cropping_flag =  [0 0 0]; %Flag for cropping images to reference channel
se = strel('disk', config.nuc_radius);

% Check if channels have the same sized images
for i = 1:length(markers)
    path_sub = path_table(path_table.markers == markers(i),:);
    tempI2 = imread(path_sub.file{1});    
    if size(tempI) ~= size(tempI2)
        cropping_flag(i) = 1;
    end
end

% Initialize matrices for storing raw images
I_raw(:) = {uint16(zeros(nrows,ncols,num_images))};

%Start parallel pool
p = gcp('nocreate');
if isempty(p) && number_of_cores>1
    parpool(min(length(markers),number_of_cores));
end

dim_adj = round([nrows ncols num_images]./s);

fprintf('Reading images and pre-processing \n')
if isequal(using_loaded_parameters,'true')
    parfor i = 1:length(markers)
        path_sub = path_table(path_table.markers == markers(i),:);
        
        for j = 1:num_images
            img = imread(path_sub.file{j});
            if cropping_flag(i) == 1
                I_raw{i}(:,:,j) = crop_to_ref(tempI,img);
            else
                I_raw{i}(:,:,j) = img;
            end
        end
    end
else
    I = I_raw;
    parfor i = 1:length(markers)
        path_sub = path_table(path_table.markers == markers(i),:);
        for j = 1:num_images
            img = imread(path_sub.file{j});
            if cropping_flag(i) == 1
                I{i}(:,:,j) = imresize(img,[nrows ncols]);
            else
                I{i}(:,:,j) = img;
            end
            I_raw{i}(:,:,j) = I{i}(:,:,j);
            if isequal(enhance_blobs,"true") && i > 1
                I_filt = imtophat(I{i}(:,:,j),se);
                I{i}(:,:,j) = I{i}(:,:,j) + I_filt*2;
            end
        end
        I{i} = im2int16(imadjustn(I{i},[0 upperThresh(i)],[],img_gamma_adj(i)));
        I{i} = imresize3(I{i},dim_adj);
    end

    % Histogram matching
    for i = 2:length(markers)
        if histogram_bins(i) > 0
            I{i} = imhistmatch(I{i},I{1},histogram_bins(i));
        end 
    end
    
    % Generate mask using reference channel
    mask = generate_sampling_mask(I{1},mask_int_threshold,lowerThresh,upperThresh);
    
    % Determine chunk positions
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

fprintf('Images loaded and pre-processed in %d seconds\n', round(toc))

%% Now perform the registration
if isequal(using_loaded_parameters,'false')
    %Initial registration on whole downsampled stack
    init_tform = cell(1,length(markers)-1);
    fprintf('Performing intial registration \n')
    fprintf('Using %d chunks\n', n_chunks)
    for i = 2:length(markers)-1 
        [init_tform{i},reg_img] = elastix(I{i+1}(:,:,300:600),I{1}(:,:,270:630),outputDir{i},{transform_path_init},'s',s,'threads',[]);
    end
    imshowpair(I{1}(:,:,350),reg_img(:,:,80))
    figure
    imshowpair(I{1}(:,:,100),I{3}(:,:,100))
        
    % Save transform parameters into cell array
    out = cell(n_chunks,length(markers)-1);
    for i = 1:n_chunks
        % Take the mask for the respective chunk to align
        mask_chunk = mask(:,:,chunk_start_adj(i):chunk_end_adj(i));
        
        % Clip the bottom and top sections in padded regions
        mask_chunk(:,:,1:round(chunk_pad/2)) = 0;
        mask_chunk(:,:,size(mask_chunk,3)-round(chunk_pad/2):size(mask_chunk,3)) = 0;
        
        % Take chunk images from reference channel
        I_ref = I{1}(:,:,chunk_start_adj(i):chunk_end_adj(i));
        
        % Save chunk start and end location with and without accounting for
        % padding
        chunk_start1 = chunk_start_adj(i);
        chunk_end1 = chunk_end_adj(i);
        
        chunk_start2 = chunk_start(i);
        chunk_end2 = chunk_end(i);

        tic
        for j = 2:length(markers)-1
            if j ~= 1
                pause(1)
            end
            % Take initial transform of whole stack and save it as an
            % elastix parameters file
            init_tform{j}.TransformParameters{1}.Size(3) = size(mask_chunk,3);
            init_tform_path = sprintf('%s/init_transform_C%d_%d_%d_%d_%s.txt',outputDir{j},j,y,x,i,sample_name);
            elastix_paramStruct2txt(init_tform_path,init_tform{j}.TransformParameters{1})

            % Perform registration on the chunk for this channel
            out{i,j} = elastix(I{j+1}(:,:,chunk_start1:chunk_end1),I_ref,...
                outputDir{j},{transform_path1{j} transform_path2{j}},'fMask', mask_chunk,...
                't0',init_tform_path,'s',s,'threads',[]);
        
            % Save chunk information and initial transform filename into out
            % structure
            out{i,j}.TransformParameters = horzcat(init_tform{j}.TransformParameters,...
                out{i,j}.TransformParameters);
            
            out{i,j}.chunk_start_adj = chunk_start1;
            out{i,j}.chunk_end_adj = chunk_end1;
            
            out{i,j}.chunk_start = chunk_start2;
            out{i,j}.chunk_end = chunk_end2;
            
            % Reset temporary directory
            rmdir(outputDir{j},'s')
            mkdir(outputDir{j})
            
            %Preview registration?
            %output = transformix(I{j+1}(:,:,chunk_start1:chunk_end1),out{i,j},s,[]);            
            %output = int16(output);
            slice = 150;
            I_topro_slice = im2uint16(I_ref(:,:,slice));
            I_transform_slice = im2uint16(output(:,:,slice));
            figure
            imshowpair(I_topro_slice*2, I_transform_slice*1.5)
        end
    fprintf('Finished registration in %d seconds\n', round(toc))
    end
    
    % Save transform parameters to variables folder
    save(fullfile(char(output_directory),'variables',sprintf('out_0%d_0%d.mat',y,x)), 'out')
    
    % Cleanup temporary directories
    cellfun(@(s) rmdir(s,'s'), outputDir)
end

%% This is where intensity adjustment should happen.
if isequal(config.adjust_intensity, "true") || isequal(config.adjust_intensity, "load")
    adj_params = config.adj_params;
   
    % Crop y_adj to match reference 
    for i = 1:length(markers)
        adj_params.y_adj{i} = crop_to_ref(adj_params.y_adj{1},adj_params.y_adj{i});
        adj_params.flatfield{i} = crop_to_ref(adj_params.flatfield{1},adj_params.flatfield{i});
        adj_params.darkfield{i} = crop_to_ref(adj_params.darkfield{1},adj_params.darkfield{i});
    end
    
    % Adjust intensities
    parfor i = 1:length(markers)
       fprintf(strcat(char(datetime('now')),'\t Applying intensity adjustments for marker: %s\n'),markers(i));
       for j = 1:num_images
           I_raw{i}(:,:,j) = apply_intensity_adjustment(I_raw{i}(:,:,j), adj_params,...
				'r',y,'c',x,'idx',i);
       end
    end
end

%% Now apply transformations to original
I_raw = apply_transformations(I_raw,out,num_images,ncols,nrows,markers);

%% Save the images
fprintf(strcat(char(datetime('now')),'\t Writing aligned images\n'));

% Create directory to store images
if exist(fullfile(char(output_directory),'aligned'),'dir') ~= 7
    mkdir(fullfile(char(output_directory),'aligned'));
end

% Create directory to store sample images
if exist(fullfile(char(output_directory),'samples'),'dir') ~= 7
    mkdir(fullfile(char(output_directory),'samples'));
end

% Save some sample images
for i = 100:100:num_images
    I_merge = cat(3,I_raw{1}(:,:,i),I_raw{2}(:,:,i),I_raw{3}(:,:,i));

    img_name = sprintf('merge_%s_0%d_0%d_aligned.tif',num2str(i,'%04.f'),y,x);
    img_path = fullfile(char(output_directory),'samples',img_name);
    imwrite(I_merge,img_path)
end

% Save aligned images
parfor j = 1:length(markers)
    for i = 1:num_images
        img_name = sprintf('%s_%s_C%d_%s_0%d_0%d_aligned.tif',sample_name,num2str(i,'%04.f'),j,markers(j),y,x);
        img_path = fullfile(char(output_directory),'aligned',img_name);
        imwrite(I_raw{j}(:,:,i),img_path)
    end
end

end


function I_raw = apply_transformations(I_raw,out,num_images,ncols,nrows,markers)
% Apply elastix transform parameters

% Get chunk sizes and locations in the stack with/without padding
n_chunks = size(out,1);
chunk_start_adj = cellfun(@(s) s.chunk_start_adj,out(:,1));
chunk_end_adj = cellfun(@(s) s.chunk_end_adj,out(:,1));

chunk_start = cellfun(@(s) s.chunk_start,out(:,1));
chunk_end = cellfun(@(s) s.chunk_end,out(:,1));

% Apply transformatin for each chunk
for i = 1:n_chunks
    if i == n_chunks && chunk_end(i) == num_images
        a = 0;
    else
        a = 1;
    end

    % Get locations for current chunk
    chunk_start1 = chunk_start_adj(i);
    chunk_end1 = chunk_end_adj(i);

    top_trim = abs(chunk_start1-chunk_start(i));
    bottom_trim = abs(chunk_end1-chunk_end(i))+a;
    
    chunk_start2 = chunk_start(i);
    chunk_end2 = chunk_end(i)-a;

    tic
    parfor j = 1:length(markers)-1
        if j ~= 1
            pause(1)
        end
        % out contains transform parameters
        out2 = out{i,j};
        
        % Adjust size back from rescaled resolution 
        size1 = out2.TransformParameters{end}.Size;
        size1 = [ncols nrows size1(3)];

        % Backward compatibility. Move init_transform field as 1st in 3
        % transform parameters
        %if length(out2.TransformParameters) == 2
        %    out2.TransformParameters = horzcat({out2.InitTransformParameters},out2.TransformParameters);
        %end
        
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
