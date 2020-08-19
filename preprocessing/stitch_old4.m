function stitch_old4(path_table,output_directory,overlap,sd,config_vars,adj_params,markers,stitch_sub_stack,alignment_table,use_middle,subtract_img_background,sift_refinement,blending_method,number_of_cores)

fprintf(strcat(char(datetime('now')),'\t Begin stitching\n'))
usfac = 10;
peaks = 5;

%Create directory for stitched images
if ~exist(char(strcat(output_directory,'/stitched')),'dir')
  mkdir(char(strcat(output_directory,'/stitched')))
end

%Generate image name grid
img_name_grid = cell(max(path_table.y),max(path_table.x),length(markers),max(path_table.z_adj));
path_table = sortrows(path_table,["z_adj","channel_num","x","y"],'ascend');

try
    img_name_grid = reshape(path_table.file,size(img_name_grid));
catch ME
    if isequal(ME.identifier,'MATLAB:getReshapeDims:notSameNumel')
        msg = strcat(char(datetime('now')),sprintf('\t Inconsistent image file information. Recalculate adjusted z and/or check configuration\n'));
        error(msg)
    end
end

nb_sections = size(img_name_grid,4);
[nrows,ncols] = size(img_name_grid(:,:,1));
start_z = 0;

%Check if only certain sub-section is to be stitched
if ~isempty(stitch_sub_stack)
    img_name_grid = img_name_grid(:,:,:,stitch_sub_stack);
    z_range = stitch_sub_stack;
    nb_sections = length(z_range);
else
    z_range = 1:size(img_name_grid,4);
end

% Create .mat file for storing stitching information
stitch_file = fullfile(output_directory,'variables','stitch_tforms.mat');
h_pos = (ncols-1)*nrows*2;
v_pos = (nrows-1)*2;
if ~exist(stitch_file,'file')
    h_stitch_tforms = zeros(h_pos,nb_sections);
    v_stitch_tforms = zeros(v_pos,nb_sections);
    save(stitch_file,'h_stitch_tforms','v_stitch_tforms','-v7.3');
elseif ~isempty(stitch_sub_stack) && exist(stitch_file,'file')
    
end

nb_img_tiles = size(img_name_grid,1)*size(img_name_grid,2);

%Determine order of images being stitched based on which Z position has the
%most signal from all tiles. If use_middle == 1, start from the middle
if isequal(use_middle,'false')
    fprintf(strcat(char(datetime('now')),'\t Calculating stitching start position\n'));
    %Check 20% of tiles
    pos = max(1,floor(linspace(1,nb_sections,ceil(nb_sections/100)))); 
    signal = zeros(length(pos),numel(img_name_grid(:,:,1)));
    for i = 1:length(pos)
        img_grid = img_name_grid(:,:,:,pos(i));
        for j = 1:numel(img_grid)
            tempI = imread(img_grid{j});
            if size(tempI,3)>1
                tempI = tempI(:,:,1);
            end
            signal(i,j) = sum(tempI(:)>config_vars.lowerThresh(1))/numel(tempI);
        end
    end
    %Find which tiles have some kind of signal
    sig_tiles = sum(signal>0.05,2);
    sig_pos = pos(sig_tiles == max(sig_tiles));
    %Choose closest to the middle from these positions
    if sig_tiles< nb_img_tiles
        small_signal = signal(:,~any(signal>0.10));
        [~,z_idx] = max(small_signal);
        z_idx = round(mean(z_idx));
        start_z = pos(z_idx);
    else
        [~,start_z] = min(abs(nb_sections/2-sig_pos));
        start_z = sig_pos(start_z);
    end
else
    %Otherwise start from the middle
    start_z = round(size(img_name_grid,4)/2);
end 

%Start parallel pool
p = gcp('nocreate');
if isempty(p) && number_of_cores>1
    %parpool(2)
end

%Begin stitching from the starting position. Stitching proceeds iteratively
%from here to top and bottom. Previous translations are used as thresholds
%for maximum translations for the current section. If threshold is exceeded
%(i.e. images are moving more than expected) the previous translation is
%used as the current section is lacking enough features (likely because
%it's at the edge of the sample
for idx2 = 1:2
    m = matfile(stitch_file,'Writable',true);
    if ~isempty(stitch_sub_stack) && exist(stitch_file,'file')
       h_tform = m.h_stitch_tforms(:,start_z+1);
       h_tform = {reshape(h_tform, 2, length(h_tform)/2)};
       v_tform = m.v_stitch_tforms(:,start_z+1)';
       v_tform = {reshape(v_tform, 2, length(v_tform)/2)};
    else
       
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
            [h_tform,v_tform]=stitch_worker(img_name_grid(:,:,:,i),output_directory,...
            overlap,pre_h_tform1,pre_v_tform1,sd,config_vars,adj_params,markers,...
            z_pos,alignment_table,subtract_img_background,sift_refinement,usfac,blending_method);
            
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
            [h_tform,v_tform]=stitch_worker(img_name_grid(:,:,:,j),output_directory,...
            overlap,pre_h_tform2,pre_v_tform2,sd,config_vars,adj_params,markers,...
            z_pos,alignment_table,subtract_img_background,sift_refinement,usfac,blending_method);
           
            %Store translations
            h_tform1 = h_tform'; v_tform1 = v_tform';
            m.h_stitch_tforms(:,z_pos) = [h_tform1{:}]'; 
            m.v_stitch_tforms(:,z_pos) = [v_tform1{:}]';
        end
    end
end

end

function [pre_h_tform,pre_v_tform] = stitch_worker(img_grid,output_directory,...
    overlap,pre_h_tform,pre_v_tform,sd,config_vars,adj_params,markers,...
    z_idx,alignment_table,subtract_img_background,sift_refinement,usfac,blending_method)

peaks1 = 3;
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
            %Read image
            A{i,j,k} = imread(img_grid{i,j,k});
            
            %Crop or pad images
            if ~isequal(size(A{i,j,k}),[img_height,img_width])
                A{i,j,k} = crop_to_ref(A{1},A{i,j,k});
            end
            
            %Background subtraction
            if isequal(subtract_img_background,'true')
                t_adj = cellfun(@(s) s(i,j), adj_params.t_adj(k));
                y_adj = adj_params.y_adj{k};
                se = strel('disk',adj_params.nuc_radius);
                upperThresh = adj_params.upperThresh(k);
                lowerThresh = adj_params.lowerThresh(k);
                gamma = adj_params.gamma(k);
                filter = adj_params.filter(k);
                
               % Crop laser width adjustment if necessary
                if length(y_adj) ~= length(adj_params.y_adj{1})
                    y_adj = crop_to_ref(adj_params.y_adj{1},adj_params.y_adj{i});
                    A{i,j,k} = crop_to_ref(A{1,1,1},A{i,j,k});
                end
                A{i,j,k} = subtract_background(A{i,j,k},se,t_adj,y_adj,lowerThresh,upperThresh,gamma,filter);
            end
            
            if ~isempty(adj_params) && ~isequal(config_vars.adjust_intensity,'false')
                A{i,j,k} = apply_intensity_adjustment(A{i,j,k},'p',adj_params,'r',i,'c',j,'idx',k);
            end
            
            %Apply alignment transforms
            if ~isempty(alignment_table) && k>1                
                alignment_table_sub = alignment_table{i,j}(alignment_table{i,j}.Reference_Z == z_idx,:);
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
v_overlap = round(img_height*overlap);
overlap_v_min = 1:v_overlap;
h_overlap = round(img_width*overlap);
overlap_h_min = 1:h_overlap;

%Sizes of optimal, fully stitched image
full_width = img_width*ncols-h_overlap*(ncols-1);
full_height = img_height*nrows-v_overlap*(nrows-1);

%Generate pixel merge weights using sigmoid function
%Horizontal
w_h = linspace(-sd,sd,h_overlap);
w_h = 1./(1 + exp(-(w_h)));

%Vertical
w_v = linspace(-sd,sd,v_overlap);
w_v = 1./(1 + exp(-(w_v)))';

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

%Calculate Horizontal Translations
for i = 1:nrows 
for j = 1:ncols-1
    %Update overlap region of the left image 
    overlap_h_max = size(B{i,1},2)-h_overlap+1:size(B{i,1},2);
    
    if i == 3 && j == 2
        aaaa =1;
    end

    %Load overlapped regions
    ref_img = B{i,:,1}(:,overlap_h_max);
    mov_img = A{i,j+1,1}(:,overlap_h_min);
    
    %Check number of bright pixels
    signal = sum(ref_img(:)>config_vars.lowerThresh(1))/numel(mov_img);

    if signal <0.005
        try 
            final_tform = affine2d([1 0 0; 0 1 0; pre_h_tform{i,j}(1) pre_h_tform{i,j}(2) 1]);
        catch
            if isnan(pre_h_tform{i,j}(1)) || isnan(pre_h_tform{i,j}(2))
                msg = sprintf('%s\t Calling NaN as expected transform. Check thresholds or recalculate start z position. \n',datetime('now'));
                error(msg)
            end
        end
    else
        %Perform phase correlation and refine with SIFT
        [pc_img,ref_img,tformPC] = calculate_phase_correlation(mov_img,ref_img,peaks1,usfac);

        if isempty(tformPC) || abs(tformPC.T(3)-pre_h_tform{i,j}(1))>limit_x || abs(tformPC.T(6)-pre_h_tform{i,j}(2))>limit_y
            fprintf(strcat(char(datetime('now')),'\t Warning: large horizontal displacement at %d x %d \n'),i,j);
            tformPC = affine2d([1 0 0; 0 1 0; pre_h_tform{i,j}(1) pre_h_tform{i,j}(2) 1]);
            pc_img = imtranslate(mov_img, [pre_h_tform{i,j}(1) pre_h_tform{i,j}(2)]);
        end
        
        %Refine using SIFT
        if isequal(sift_refinement,'true')
            [tformSIFT] = sift_refinement_worker(pc_img,ref_img,single(w_h));
            final_tform = affine2d(tformPC.T*tformSIFT.T);
        else
            final_tform = tformPC;
        end
        
        %If not able to calculate transform, use previous transform
        if abs(final_tform.T(3)-pre_h_tform{i,j}(1))>limit_x || abs(final_tform.T(6)-pre_h_tform{i,j}(2))>limit_y
            final_tform = affine2d([1 0 0; 0 1 0; pre_h_tform{i,j}(1) pre_h_tform{i,j}(2) 1]);
        end
    end

    %Adjust horizontally overlapped pixels based on translations
    ref_fixed2 = imref2d([img_height img_width+ceil(final_tform.T(3))]);

    %Transform and merge images (faster to for loop on each channel)
    for k = 1:nchannels
        reg_img = imwarp(A{i,j+1,k},final_tform,'OutputView',ref_fixed2,'FillValues',0);
        B{i,k} = blend_images(reg_img,B{i,k},w_h,overlap_h_min,overlap_h_max,...
            blending_method(k),'horizontal');  
    end
    
    %Save translation
    pre_h_tform{i,j} = [final_tform.T(3), final_tform.T(6)];
end
end

%Crop horizontally stitched images to minimum width
[min_width] = min(cellfun(@(s) size(s,2),B(:,1)));
B = cellfun(@(s) s(:,1:min_width), B,'UniformOutput',false);
I = B(1,:);

%Check for previous translations and set limits
if ~iscell(pre_v_tform)
    limit_x = NaN;
    limit_y = NaN;
    pre_v_tform = repmat({[NaN,NaN]},[1 length(B)-1]);
else
    limit_x = 10;
    limit_y = 10; 
end

for i = 1:length(B)-1
    %Update overlap region of the top image 
    overlap_v_max = size(I{1},1)-v_overlap+1:size(I{1},1);
    
    %Load overlapped regions
    ref_img = I{1}(overlap_v_max,1:min_width);
    mov_img = B{i+1,1}(overlap_v_min,1:min_width);

    signal = sum(ref_img(:)>config_vars.lowerThresh(1))/numel(mov_img);
    
    %When there is no intensity, use previous translation
    if signal <0.02
        final_tform = affine2d([1 0 0; 0 1 0; pre_v_tform{i}(1) pre_v_tform{i}(2) 1]);
    else
        %Perform phase correlation and refine with SIFT
        [pc_img,ref_img,tformPC] = calculate_phase_correlation(mov_img,ref_img,peaks1,usfac);
        
        if isempty(tformPC) || abs(tformPC.T(3)-pre_v_tform{i}(1))>limit_x || abs(tformPC.T(6)-pre_v_tform{i}(2))>limit_y
            fprintf(strcat(char(datetime('now')),'\t Warning: large vertical displacement at %d\n'),i);
            tformPC = affine2d([1 0 0; 0 1 0; pre_v_tform{i}(1) pre_v_tform{i}(2) 1]);
            pc_img = imtranslate(mov_img,[pre_v_tform{i}(1) pre_v_tform{i}(2)]);
        end
        
        %Refine using SIFT
        if isequal('sift_refinement','true')
            [tformSIFT] = sift_refinement_worker(pc_img,ref_img,w_v);
            final_tform = affine2d(tformPC.T*tformSIFT.T);
        else
            final_tform = tformPC;
        end

        %If not able to calculate transform, use previous transform
        if abs(final_tform.T(3)-pre_v_tform{i}(1))>limit_x || abs(final_tform.T(6)-pre_v_tform{i}(2))>limit_y
            final_tform = affine2d([1 0 0; 0 1 0; pre_v_tform{i}(1) pre_v_tform{i}(2) 1]);
        end
    end
    
    %Adjust horizontally overlapped pixels based on translations
    ref_fixed2 = imref2d([img_height+ceil(final_tform.T(6)) size(I{1},2)]);
        
    %Transform image
    for k = 1:nchannels
        reg_img = imwarp(B{i+1,k},final_tform,'OutputView',ref_fixed2,'FillValues',0,'SmoothEdges',true);
        
        %Adjust intensity again?
        %adj_factor = median(I{k}(overlap_v_max,:),'all')/median(reg_img(overlap_v_min,:),'all');
        %adj_factor = prctile(I{k}(overlap_v_max,:),75)/prctile(reg_img(overlap_v_min,:),75);
        
        %reg_img = reg_img * adj_factor;
        
        I{k} = blend_images(reg_img,I{k},w_v,overlap_v_min,overlap_v_max,...
            blending_method(k),'vertical'); 
    end
    
    %Save translation
    pre_v_tform{i} = [final_tform.T(3), final_tform.T(6)];
end

%Subtract background
if isequal(subtract_img_background,'true')
    I = cellfun(@(s) subtract_background(s,strel('disk',adj_params.nuc_radius*2),1,1,0,1,1,'true'),I,'UniformOutput',false);
end

%Crop or pad images based on ideal size
I = cellfun(@(s) crop_to_ref(zeros(full_height,full_width),s),I,'UniformOutput',false);

%Save images as individual channels (will be large)
for i = 1:nchannels
    img_name = sprintf('%s_%s_C%d_%s_stitched.tif',config_vars.sample_name,num2str(z_idx,'%04.f'),i,markers(i));
    img_path = fullfile(char(output_directory),'stitched',img_name);
    imwrite(uint16(I{i}),img_path)
end

end

function [tform] = sift_refinement_worker(mov_img,ref_img,weight)
%Detect SIFT features
[f1, d1] = vl_sift(ref_img,'PeakThresh',10,'EdgeThresh',2);
[f2, d2] = vl_sift(mov_img,'PeakThresh',10,'EdgeThresh',2);

%Match SIFT features
[matches] = vl_ubcmatch(d1, d2);

%Take x,y positions of matched points
x1 = f1(1:2,matches(1,:));
x2 = f2(1:2,matches(2,:));

%Calculate distance between matches
x3  = sqrt(sum((x1 - x2).^2));

%Remove point far away from each other
matches(:,abs(x3)>10)=[];

%Display number of matches
numMatches = size(matches,2);

X1 = f1(1:2,matches(1,:)); X1(3,:) = 1;
X2 = f2(1:2,matches(2,:)); X2(3,:) = 1;

%Add non-linear blend weight
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
    
%Instead of RANSAC, use just the average feature locations since we're
%already removing outliers based on distance
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

function ref_img = blend_images(mov_img,ref_img,w,overlap_min,overlap_max,...
    blending_method,direction)

switch blending_method
    case 'sigmoid'
        if isequal(direction,'horizontal')
            %Perform non-linear weight
            w(sum(mov_img(:,overlap_min))==0)=0;
            inv_w = 1-w;
            ref_img(:,overlap_max) = ref_img(:,overlap_max).*inv_w +...
                mov_img(:,overlap_min).*w;
            mov_img(:,overlap_min) = [];

            %Concatanate images
            ref_img = horzcat(ref_img,mov_img);
        else

            %Perform non-linear weight
            w(sum(mov_img(overlap_min),2)==0)=0;
            inv_w = 1-w;            
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
    case 'mean'
        if isequal(direction,'horizontal')
            %Take max in the overlapping region
            ref_img(:,overlap_max) = (ref_img(:,overlap_max)+mov_img(:,overlap_min))/2;
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

