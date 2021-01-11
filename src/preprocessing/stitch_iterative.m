function stitch_iterative(config, path_table)
%--------------------------------------------------------------------------
% Perform iterative 2D image stitching using phase correlation with
% optional refinement using SIFT. Image slices positions are presumed to be
% adjusted along the z dimension to find the optimal correspondence.
% Therefore this program will start from roughly the middle z positions and
% work its way out towards the first and last slice.
%--------------------------------------------------------------------------

usfac = 10;
peaks = 5;

% Create directory for stitched images
if ~exist(fullfile(config.output_directory,'stitched'),'dir')
    mkdir(fullfile(config.output_directory,'stitched'))
end

% If using SIFT, check if vl_feat toolbox exists
if isequal(config.sift_refinement,"true")
    if exist('vl_sift','file') == 0
        error("Could not call vl_feat toolbox. Download here: https://www.vlfeat.org/install-matlab.html")
    else
        fprintf("%s\t Using SIFT refinement \n",datetime('now'))
    end
end

% Check if only certain channels to be stitched
if isempty(config.stitch_sub_channel)
    config.stitch_sub_channel = 1:length(config.markers);
end

% Generate image name grid
img_name_grid = cell(max(path_table.y),max(path_table.x),length(config.stitch_sub_channel),max(path_table.z_adj));
path_table = sortrows(path_table,["z_adj","channel_num","x","y"],'ascend');

% Arrange images into correct positions
try
    img_name_grid = reshape(path_table.file,size(img_name_grid));
catch ME
    if isequal(ME.identifier,'MATLAB:getReshapeDims:notSameNumel')
        error("Inconsistent image file information. Recalculate adjusted z and/or check configuration\n")
    end
end

% Count number of sections
nb_sections = size(img_name_grid,4);
[nrows,ncols] = size(img_name_grid(:,:,1));
fprintf("%s\t Begin stitching %d slices \n",datetime('now'),nb_sections)

%Check if only certain sub-section is to be stitched
if ~isempty(config.stitch_sub_stack)
    img_name_grid = img_name_grid(:,:,:,config.stitch_sub_stack);
    z_range = config.stitch_sub_stack;
    nb_sections = length(z_range);
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

% Create .mat file for storing stitching information
stitch_file = fullfile(config.output_directory,'variables','stitch_tforms.mat');
h_pos = (ncols-1)*nrows*2;
v_pos = (nrows-1)*2;
if ~exist(stitch_file,'file')
    h_stitch_tforms = zeros(h_pos,nb_sections);
    v_stitch_tforms = zeros(v_pos,nb_sections);
    save(stitch_file,'h_stitch_tforms','v_stitch_tforms','-v7.3');    
end

% Determine order of images being stitched based on which z position has
% the most bright features in the tile with lowest signal. Otherwise start
% from middle or from other specified z 
if isempty(config.stitch_start_slice)
    fprintf("%s\t Calculating stitching start position \n",datetime('now'))
    % Check 10% of tiles
    pos = max(1,floor(linspace(1,nb_sections,round(nb_sections*0.1)))); 
    signal = zeros(length(pos),numel(img_name_grid(:,:,1)));
    for i = 1:length(pos)
        img_grid = img_name_grid(:,:,1,pos(i));
        for j = 1:numel(img_grid)
            tempI = imread(img_grid{j});
            if size(tempI,3)>1
                tempI = tempI(:,:,1);
            end
            signal(i,j) = sum(tempI(:)>config.signalThresh(1)*65535)/numel(tempI);
        end
    end
    % Find which tile has the lowest mean signal 
    signal_ave = mean(signal,1);
    [~,z_idx] = max(signal(:,signal_ave == min(signal_ave)));
    start_z = pos(z_idx);    
    
    % Find which tiles have some kind of signal
    %sig_tiles = sum(signal>0.05,2);
    %sig_pos = pos(sig_tiles == max(sig_tiles));
    % Choose closest to the middle from these positions
    %if sig_tiles < nb_img_tiles
    %    small_signal = signal(:,~any(signal>0.10));
    %    [~,z_idx] = max(small_signal);
    %    z_idx = round(mean(z_idx));
    %    start_z = pos(z_idx);
    %else
    %    [~,start_z] = min(abs(nb_sections/2-sig_pos));
    %    start_z = sig_pos(start_z);
    %end
elseif isnumeric(config.stitch_start_slice) && isempty(config.stitch_sub_stack)
    % Start from specified z position
    start_z = config.stitch_start_slice;
else
    % Otherwise start from the middle
    start_z = round(size(img_name_grid,4)/2);
end 

% Start parallel pool
if nb_sections > 20
try
    p = gcp('nocreate');
catch
    p = 1;
end
else
    p = 1;
end

if isempty(p)
    parpool(2)
end

%Begin stitching from the starting position. Stitching proceeds iteratively
%from here to top and bottom. Previous translations are used as thresholds
%for maximum translations for the current section. If threshold is exceeded
%(i.e. images are moving more than expected) the previous translation is
%used as the current section is lacking enough features (likely because
%it's at the edge of the sample
for idx2 = 1:2
    m = matfile(stitch_file,'Writable',true);
    if ~isempty(config.stitch_sub_stack)
        %%%%%%%Wrong
       h_tform = m.h_stitch_tforms(:,start_z);
       h_tform = {reshape(h_tform, 2, length(h_tform)/2)};
       v_tform = m.v_stitch_tforms(:,start_z)';
       v_tform = {reshape(v_tform, 2, length(v_tform)/2)};
    else
        h_tform = []; v_tform = [];
    end
    h_tform = []; v_tform = [];
    
    if idx2 == 1
        %From middle to top
        for i = fliplr(1:start_z)
            z_pos = z_range(i);
            pre_h_tform1 = h_tform;
            pre_v_tform1 = v_tform;
            
            %Print image being stitched
            fprintf(strcat(char(datetime('now')),'\t Stitching image %d \n'),z_range(i));
            [h_tform,v_tform]=stitch_worker(img_name_grid(:,:,:,i),pre_h_tform1,pre_v_tform1,...
                config,z_pos,adj_params,usfac,peaks);
            
            %Store translations
            h_tform1 = h_tform'; v_tform1 = v_tform';
            m.h_stitch_tforms(:,z_pos) = [h_tform1{:}]'; 
            m.v_stitch_tforms(:,z_pos) = [v_tform1{:}]';
        end
    else
        m = matfile(stitch_file,'Writable',true);
        %From middle to bottom
        for j = start_z+1:nb_sections
            z_pos = z_range(j);
            pre_h_tform2 = h_tform;
            pre_v_tform2 = v_tform;

            %Print image being stitched
            fprintf(strcat(char(datetime('now')),'\t Stitching image %d \n'),z_range(j));
            [h_tform,v_tform]=stitch_worker(img_name_grid(:,:,:,j),pre_h_tform2,pre_v_tform2,...
                config,z_pos,adj_params,usfac,peaks);
           
            %Store translations
            h_tform1 = h_tform'; v_tform1 = v_tform';
            m.h_stitch_tforms(:,z_pos) = [h_tform1{:}]'; 
            m.v_stitch_tforms(:,z_pos) = [v_tform1{:}]';
        end
    end
end

end

function [pre_h_tform,pre_v_tform] = stitch_worker(img_grid,pre_h_tform,pre_v_tform,config,z_idx,adj_params,usfac,peaks)
% Worker for stitching_iterative function

% Defaults
border_pad = config.border_pad; % Border cropping along edges
min_overlap = 50;    % Minimum overlapping region in pixels

% Image grid info
[nrows,ncols,nchannels] = size(img_grid);

% Read images, adjust intensities, apply translations for multichannel
A = read_stitching_grid(img_grid,config.stitch_sub_channel,config.markers,...
    adj_params,config.alignment_table);
[img_height,img_width] = size(A{1});

% Calculate overlaps in pixels
h_overlap = round(img_width*config.overlap);
v_overlap = round(img_height*config.overlap);
x0 = round(h_overlap/2);
y0 = round(v_overlap/2);

% Calculate extended overlap if presumed overlap is small
if h_overlap < min_overlap
    ext_adj_h = min([min_overlap, img_width]) - h_overlap;
else
    ext_adj_h = 0;
end
if v_overlap < min_overlap
    ext_adj_v = min([min_overlap, img_height]) - v_overlap;
else
    ext_adj_v = 0;
end
overlap_h_min = 1:h_overlap + ext_adj_h;
overlap_v_min = 1:v_overlap + ext_adj_v;

% Sizes of optimal, fully stitched image
full_width = img_width*ncols-h_overlap*(ncols-1);
full_height = img_height*nrows-v_overlap*(nrows-1);

%Check for previous translations and set limits
if ~iscell(pre_h_tform)
    limit_x = NaN;
    limit_y = NaN;
    pre_h_tform = repmat({[NaN,NaN]},[nrows,ncols-1]);
else
    limit_x = 5;
    limit_y = 5; 
end

%Save first column before horizontal stitching
B = A(:,1,:);

% Calculate horizontal translations
for i = 1:nrows 
    if ncols == 1
        continue
    end
    for j = 1:ncols-1
        % Update overlap region of the left image 
        overlap_h_max = size(B{i,1},2)-length(overlap_h_min)+1:size(B{i,1},2);
        %overlap_h_max = img_width - length(overlap_h_min)+1:img_width;

        % Load overlapped regions
        ref_img = B{i,:,1}(:,overlap_h_max);
        %ref_img = A{i,j,1}(:,overlap_h_max);
        mov_img = A{i,j+1,1}(:,overlap_h_min);

        % Check number of bright pixels
        signal = sum(ref_img(:)>config.lowerThresh(1))/numel(mov_img);

        if signal <0.005
            try 
                final_tform = affine2d([1 0 0; 0 1 0; pre_h_tform{i,j}(1) pre_h_tform{i,j}(2) 1]);
            catch
                if isnan(pre_h_tform{i,j}(1)) || isnan(pre_h_tform{i,j}(2))
                    error("%s\t Calling NaN as expected transform. Check thresholds or recalculate start z position. \n",datetime('now'))
                end
            end
        else
            % Perform phase correlation
            [pc_img,ref_img,tformPC] = calculate_phase_correlation(mov_img,ref_img,peaks,usfac);

            % Warn user if large translation
            if isempty(tformPC) || abs(tformPC.T(3)-pre_h_tform{i,j}(1))>limit_x+ext_adj_h || abs(tformPC.T(6)-pre_h_tform{i,j}(2))>limit_y
                fprintf(strcat(char(datetime('now')),'\t Warning: large horizontal displacement at %d x %d \n'),i,j);
                tformPC = affine2d([1 0 0; 0 1 0; pre_h_tform{i,j}(1) pre_h_tform{i,j}(2) 1]);
                pc_img = imtranslate(mov_img, [pre_h_tform{i,j}(1) pre_h_tform{i,j}(2)]);
            end

            % Refine using SIFT
            if isequal(config.sift_refinement,'true')
                tformSIFT = sift_refinement_worker(pc_img,ref_img);
                final_tform = affine2d(tformPC.T*tformSIFT.T);
            else
                final_tform = tformPC;
            end

            % If not able to calculate transform, use previous transform
            if abs(final_tform.T(3)-pre_h_tform{i,j}(1))>limit_x+ext_adj_h || abs(final_tform.T(6)-pre_h_tform{i,j}(2))>limit_y
                final_tform = affine2d([1 0 0; 0 1 0; pre_h_tform{i,j}(1) pre_h_tform{i,j}(2) 1]);
            end
        end
        
        % If overlap is extended, resize to to expected overlapping region
        if ext_adj_h > 0
            final_tform.T(3) = final_tform.T(3) - ext_adj_h;
            overlap_h_min1 = overlap_h_min(1:h_overlap);
        else
            overlap_h_min1 = overlap_h_min;
        end
        
        % Create adjusted blending weights
        x1 = x0 - final_tform.T(3);
        if any(arrayfun(@(s) isequal(s,"sigmoid"),config.blending_method))
            w_h_adj = 1-1./(1 + exp(config.sd*(overlap_h_min1-x1)));
        else
            w_h_adj = overlap_h_min1/h_overlap;
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
        
        % Save translation
        pre_h_tform{i,j} = [final_tform.T(3), final_tform.T(6)];
        
        % Transform and merge images (faster to use for loop on each channel)
        ref_fixed2 = imref2d([img_height img_width+floor(final_tform.T(3))]);
        for k = 1:nchannels
            reg_img = imwarp(A{i,j+1,k},final_tform,'OutputView',ref_fixed2,'FillValues',0,'SmoothEdges',true);
            
            %Adjust intensity again?
            if isequal(config.adjust_tile_position,"true")
                adj_factor = prctile(B{i,k}(:,overlap_h_max),75,'all')/prctile(reg_img(:,overlap_h_min),75,'all');
                reg_img = reg_img * adj_factor;
            end            
            B{i,k} = blend_images(reg_img,B{i,k},false,config.blending_method(k),w_h_adj);  
        end
        
    end
end  

%disp(final_tform.T)
%figure; imshow(imadjust(uint16(B{1,1}))*2)

% Crop horizontally stitched images to minimum width
min_width = min(cellfun(@(s) size(s,2),B(:,1)));
B = cellfun(@(s) s(:,1:min_width), B,'UniformOutput',false);
I = B(1,:);

% Check for previous translations and set limits
if ~iscell(pre_v_tform)
    limit_x = NaN;
    limit_y = NaN;
    pre_v_tform = repmat({[NaN,NaN]},[1 length(B)-1]);
else
    limit_x = 10;
    limit_y = 10; 
end

% Do blended row tiles
for i = 1:length(B)-1
    % Update overlap region of the top image 
    overlap_v_max = size(I{1},1)-length(overlap_v_min)+1:size(I{1},1);
    %overlap_v_max = img_height - length(overlap_v_min)+1:img_height;
    
    % Load overlapped regions
    ref_img = I{1}(overlap_v_max,1:min_width);
    %ref_img = B{i,1}(overlap_v_max,1:min_width);
    mov_img = B{i+1,1}(overlap_v_min,1:min_width);

    signal = sum(ref_img(:)>config.lowerThresh(1))/numel(mov_img);
    
    % When there is no intensity, use previous translation
    if signal <0.02
        final_tform = affine2d([1 0 0; 0 1 0; pre_v_tform{i}(1) pre_v_tform{i}(2) 1]);
    else
        % Perform phase correlation and refine with SIFT
        [pc_img,ref_img,tformPC] = calculate_phase_correlation(mov_img,ref_img,peaks,usfac);
        
        if isempty(tformPC) || abs(tformPC.T(3)-pre_v_tform{i}(1))>limit_x || abs(tformPC.T(6)-pre_v_tform{i}(2))>limit_y+ext_adj_v
            fprintf(strcat(char(datetime('now')),'\t Warning: large vertical displacement at %d\n'),i);
            tformPC = affine2d([1 0 0; 0 1 0; pre_v_tform{i}(1) pre_v_tform{i}(2) 1]);
            pc_img = imtranslate(mov_img,[pre_v_tform{i}(1) pre_v_tform{i}(2)]);
        end
        
        %Refine using SIFT
        if isequal('sift_refinement','true')
            [tformSIFT] = sift_refinement_worker(pc_img,ref_img);
            final_tform = affine2d(tformPC.T*tformSIFT.T);
        else
            final_tform = tformPC;
        end

        %If not able to calculate transform, use previous transform
        if abs(final_tform.T(3)-pre_v_tform{i}(1))>limit_x || abs(final_tform.T(6)-pre_v_tform{i}(2))>limit_y+ext_adj_v
            final_tform = affine2d([1 0 0; 0 1 0; pre_v_tform{i}(1) pre_v_tform{i}(2) 1]);
        end
    end
    
    % Adjust horizontally overlapped pixels based on translations
    %ref_fixed2 = imref2d([img_height+ceil(final_tform.T(6)) size(I{1},2)]);
    %ref_fixed2 = imref2d([img_height size(I{1},2)]);
    
    % If overlap is extended, resize to to expected overlapping region
    if ext_adj_v > 0
        final_tform.T(6) = final_tform.T(6) - ext_adj_v;
        overlap_v_min1 = overlap_v_min(1:v_overlap);
    else
        overlap_v_min1 = overlap_v_min;
    end
    
    % Create adjusted blending weights
    y1 = y0 - final_tform.T(6);
    if isequal(config.blending_method(k),"sigmoid")
        w_v_adj = 1-1./(1 + exp(config.sd*(overlap_v_min1-y1)))';
    else
        w_v_adj = (overlap_v_min1/v_overlap)';
    end
    min_w_v_adj = min(w_v_adj(w_v_adj>0));
    w_v_adj(w_v_adj>=min_w_v_adj) = (w_v_adj(w_v_adj>=(min_w_v_adj)) - min_w_v_adj)./(1-min_w_v_adj);

    % Clip ends based on vertical translation where image intensity
    % is 0
    if final_tform.T(6) > 0
        disp(final_tform.T(6))
        adj_top = max(border_pad,ceil(final_tform.T(6)));
        w_v_adj(1:adj_top) = 0;
    else
        disp(final_tform.T(6))
        adj_bottom = max(border_pad,ceil(abs(final_tform.T(6))));
        w_v_adj(1:adj_bottom) = 0;
    end
        
    %Transform images
    ref_fixed2 = imref2d([img_height+floor(final_tform.T(6)) size(I{1},2)]);
    
    for k = 1:nchannels
        reg_img = imwarp(B{i+1,k},final_tform,'OutputView',ref_fixed2,'FillValues',0,'SmoothEdges',true);
        %Adjust intensity again?
        if isequal(config.adjust_tile_position,"true")
            adj_factor = prctile(I{k}(overlap_v_max,:),75,'all')/prctile(reg_img(overlap_v_min,:),75,'all');
            reg_img = reg_img * adj_factor;
        end
        I{k} = blend_images(reg_img,I{k},false,config.blending_method(k),w_v_adj); 
    end
    
    %Save translation
    pre_v_tform{i} = [final_tform.T(3), final_tform.T(6)];
end

% Return if only calculating parameters
if isequal(config.save_images,'false')
    return
end

% Postprocess the image with various filters, background subtraction, etc.
c_idx = config.stitch_sub_channel;
for i = 1:length(c_idx)
    I{i} = postprocess_image(config, I{i}, c_idx(i));
end

%Crop or pad images based on ideal size
I = cellfun(@(s) crop_to_ref(zeros(full_height,full_width),s),I,'UniformOutput',false);

figure; imshow(imadjust(uint16(I{1}*50)))

%Save images as individual channels (will be large)
for i = 1:length(c_idx)
    img_name = sprintf('%s_%s_C%d_%s_stitched.tif',...
        config.sample_name,num2str(z_idx,'%04.f'),c_idx(i),config.markers(c_idx(i)));
    img_path = fullfile(char(config.output_directory),'stitched',img_name);
    imwrite(uint16(I{i}),img_path)
end

end

function tform = sift_refinement_worker(mov_img,ref_img)

% Defaults
min_distance = 10;

% Crop to overlapping region
if size(ref_img,1)>size(ref_img,2)
    % Horizontal
    idx = find(any(mov_img,1));
    mov_img = mov_img(:,min(idx):max(idx));
    ref_img = ref_img(:,min(idx):max(idx));
    ncols = size(mov_img,2);
    x = -floor(ncols/2):floor(ncols/2);
    weight = -x.^2+1;
    weight = (weight - min(weight))./(max(weight) - min(weight));
else
    % Vertical
    idx = find(any(mov_img,2));
    mov_img = mov_img(min(idx):max(idx),:);
    ref_img = ref_img(min(idx):max(idx),:);
    nrows = size(mov_img,1);
    x = -floor(nrows/2):floor(nrows/2);
    weight = (-x.^2+1)';
    weight = (weight - min(weight))./(max(weight) - min(weight));
end

% Detect SIFT features
[f1, d1] = vl_sift(ref_img,'PeakThresh',10,'EdgeThresh',2);
[f2, d2] = vl_sift(mov_img,'PeakThresh',10,'EdgeThresh',2);

% Match SIFT features
matches = vl_ubcmatch(d1, d2);

% Take x,y positions of matched points
x1 = f1(1:2,matches(1,:));
x2 = f2(1:2,matches(2,:));

% Calculate distance between matches
x3  = sqrt(sum((x1 - x2).^2));

% Remove point far away from each other
matches(:,abs(x3)>min_distance)=[];

% Display number of matches
numMatches = size(matches,2);

X1 = f1(1:2,matches(1,:)); X1(3,:) = 1;
X2 = f2(1:2,matches(2,:)); X2(3,:) = 1;

% Add non-linear blend weight to weigh matches in the middle of the images,
% (where the blending seams occur) more strongly
w = 1+weight.*((1-weight)/0.25);

if size(ref_img,1)>size(ref_img,2)
    w = w(round(X1(1,:)));
else
    w = w(round(X1(2,:)))';
end

%For plotting points
%imshow(imadjust(uint16(ref_img)))
%h1 = vl_plotframe(X1)';
%h2 = vl_plotframe(X2)';

%set(h1,'color','k');
%set(h2,'color','y');

if numMatches > 3
    % Instead of RANSAC, use just the average feature locations since we're
    % already removing outliers based on distance
    X3(1,:) = (X1(1,:)-X2(1,:)).*w;
    X3(2,:) = (X1(2,:)-X2(2,:)).*w;

    x = sum(X3(1,:))/sum(w);
    y = sum(X3(2,:))/sum(w);

    tform = affine2d([1 0 0; 0 1 0; x y 1]);
else
    fprintf(strcat(char(datetime('now')),'\t Not enough matches to use SIFT, attempting image registration\n'))
    metric = registration.metric.MattesMutualInformation;
    optimizer = registration.optimizer.RegularStepGradientDescent;
       
    metric.NumberOfSpatialSamples = 500;
    metric.NumberOfHistogramBins = 50;
       
    optimizer.GradientMagnitudeTolerance = 1.00000e-04;
    optimizer.MinimumStepLength = 1.00000e-05;
    optimizer.MaximumStepLength = 1.00000e-01;
    optimizer.MaximumIterations = 100;
    optimizer.RelaxationFactor = 0.5;

    %tform = imregtform(mov_img,ref_img,'translation',optimizer,metric,'PyramidLevels',2);
    tform = affine2d([1 0 0; 0 1 0; 0 0 1]);
end
end
