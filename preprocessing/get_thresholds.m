function [lowerThresh, upperThresh] = get_thresholds(stack,config)
%--------------------------------------------------------------------------
% Get thresholds without creating a structure for intensity adjustments
%--------------------------------------------------------------------------

% Needs updating. This procedure doesn't match measure_images

low_prct = 5;   % Low percentile for sampling background pixels
high_prct = 95; % High percentile for sampling bright pixels
image_sampling = 0.02;   % Fraction of all images to sample

% Get upper and lower intensity thresholds from a small subset of images in
% a stack
overlap = config.overlap;

x_tiles = length(unique(stack.x));
y_tiles = length(unique(stack.y));

%Count number of images and measure image dimensions
nb_images = height(stack)/(x_tiles*y_tiles);

% Get image positions
s = round(image_sampling*nb_images);
img_range = round(linspace(1,nb_images,s));
overlaps = length(img_range)*((x_tiles-1)*(y_tiles)+(y_tiles-1)*(x_tiles));

%Read image size
tempI = imread(stack.file{1});
[nrows, ncols] = size(tempI); 

%Determine overlap region
%Horizontal image overlap regions
overlap_min_h = {[1,nrows],[1,round(ncols*overlap)]};
overlap_max_h = {[1,nrows],[ncols-overlap_min_h{2}(2)+1,ncols]};

%Vertical image overlap regions
overlap_min_v = {[1,round(nrows*overlap)],[1,ncols]};
overlap_max_v = {[nrows-overlap_min_v{1}(2)+1,nrows],[1,ncols]};

%Initialize measurement vectors
p_low = zeros(1,overlaps);
p_high = zeros(1,overlaps);
stdev = zeros(1,overlaps);
img_idx = 1;

%Measure pairwise horizontal intensity differences
fprintf(strcat(char(datetime('now')),'\t Measuring Between Tile Differences Horizontally\n'));
for i = 1:y_tiles
for j = 1:x_tiles-1
    for k = img_range
        %Read image regions where tiles should overlap
        file_left = stack(stack.y == i & stack.x == j & stack.z == k,:);
        file_right = stack(stack.y == i & stack.x == j+1 & stack.z == k,:);
        
        I_left = double(imread(file_left.file{1},'PixelRegion',overlap_max_h));
        I_right = double(imread(file_right.file{1},'PixelRegion',overlap_min_h));     
        
        %Measure 1 percentile of all pixels. This gives rough background
        %and upper intensity
        p_low(img_idx) = mean([prctile(I_left(:),1) prctile(I_right(:),1)]);
        p_high(img_idx) = max(max([I_left(:) I_right(:)]));
        stdev(img_idx) = std2([I_left I_right]);
        
        %Update vector indices
        img_idx = img_idx+1;
    end
end
end

%Measure pairwise vertical intensity differences
fprintf(strcat(char(datetime('now')),'\t Measuring Between Tile Differences Vertically\n'));

for i = 1:y_tiles-1
    for k = img_range
        %Read image regions where tiles should overlap
        file_top = stack(stack.y == i & stack.z == k,:);
        file_bottom = stack(stack.y == i+1 & stack.z == k,:);
        
        for kk = 1:height(file_top)
        
        I_top = double(imread(file_top.file{kk},'PixelRegion',overlap_max_v));
        I_bottom = double(imread(file_bottom.file{kk},'PixelRegion',overlap_min_v));
                  
        %Measure 1 percentile of all pixels. This gives rough background
        p_low(img_idx) = mean([prctile(I_top(:),1) prctile(I_bottom(:),1)]);
        p_high(img_idx) = max(max([I_top(:), I_bottom(:)]));
        stdev(img_idx) = std2([I_top I_bottom]);
            
        %Update vector indices
        img_idx = img_idx+1;
        end
    end
end

%Calculate lower threshold
lowerThresh = round(median(p_low)+1*median(stdev));
upperThresh = min(prctile(p_high,90),65535);

end
