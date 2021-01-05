function preprocess_3dunet(res,file_delimiter,n)
% Preprocess images for running 3dunet
% Normalize images
% Identify nuclei borders, fills, and centroids
% Save images in preprocessed directory
%img_location = './Updated Training Samples/f121_images_224_64';
%cen_location = './Updated Training Samples/c121_centroids_361_102/finished';
%cen_location = [];

img_location = sprintf('./data/raw/%s/c%s_images_final_224_64',res,res);
cen_location = sprintf('./data/raw/%s/c%s_cen_final_224_64',res,res);

patch_size = [224 224 64];
split_size = [112 112 32];
final_size = [112 112 32];

patch_resolution = [1.21, 1.21, 4];
capture_resolution = [0.75, 0.75, 2.5];
save_directory = './data/preprocessed';
split_images = true;
rewrite_centroids = false;
write_edges = false;
erode_blobs = true;
normalize_images = true;
remove_background = false;
%file_delimiter = 'condensed';
size_threshold = 0;
pix_threshold = 0.4;

% Create preprocessed directory
if ~exist(save_directory,'dir')
    mkdir(save_directory)
elseif ~isempty(dir(save_directory))
    rmdir(save_directory,'s')
    mkdir(save_directory)
end

% Load images, normalize, and save images
files = dir(img_location);
file_idx = arrayfun(@(s) contains(s.name,'.nii'),files);
img_files = files(file_idx);

if ~isempty(cen_location)
    files = dir(cen_location);
    file_idx = arrayfun(@(s) contains(s.name,'.nii'),files);
    cen_files = files(file_idx);
else
    cen_files = [];
end

if ~isempty(file_delimiter)
    file_idx = arrayfun(@(s) contains(s.name,file_delimiter),files);
    cen_files = files(file_idx);
end
cen_files = cen_files(1:n);

n_cells = 0;
if ~isempty(cen_files)
    for n = 1:length(cen_files)
    cen_name = strsplit(cen_files(n).name,'.');
    fprintf('%s\n',cen_name{1})
    
    % Read centroid image
    cen_path = fullfile(cen_location,cen_files(n).name);
    L = read_nii(cen_path);
    
    % Read image
    img_idx = arrayfun(@(s) string(s.name(1:5)) == string(cen_name{1}(1:5)),...
        img_files);
    img_path = fullfile(img_location,img_files(img_idx).name);
    fprintf('%s\n',img_path)

    I = read_nii(img_path);
    I = uint16(I);
    
    if size(I,1) > patch_size(1)
        I = I(1:patch_size(1),1:patch_size(2),1:patch_size(3));
        L = L(1:patch_size(1),1:patch_size(2),1:patch_size(3));
    end
    
    % Normalize image
    if normalize_images
        I = normalize_image(I,remove_background);
    end
        
    % Expand centroids
    if rewrite_centroids
        L = rewrite_centroid(L,I,capture_resolution,patch_resolution);
    else
        L = imresize3(L,patch_size,'Method','nearest');
    end
    L = uint8(L>0);
    
    % Check orientatioin of centroids to image
    [L,~] = check_orientation(L,I);
    
    % Erode centroid blobs
    if erode_blobs
        L = erode_blob(L,I,size_threshold,pix_threshold);
    end
    
    % Save edge pixel values as 2
    if write_edges
       L = write_edge(L); 
    end
    
    % Count cells
    labels = bwconncomp(L,6);
    n_cells = n_cells + labels.NumObjects;
    
    % Count percent of positive pixels
    pixels(n) = sum(L(:)>0)/numel(L(:));
    
    % Write images
    if split_images
        n_chunks = patch_size./split_size;
        
        x_chunks = floor(linspace(0,patch_size(1),n_chunks(1)+1));
        y_chunks = floor(linspace(0,patch_size(2),n_chunks(2)+1));
        z_chunks = floor(linspace(0,patch_size(3),n_chunks(3)+1));

        a = 1;
        for x = 1:n_chunks(1)
            for y = 1:n_chunks(2)
                for z = 1:n_chunks(3)
                    xi = x_chunks(x)+1:x_chunks(x+1);
                    yi = y_chunks(y)+1:y_chunks(y+1);
                    zi = z_chunks(z)+1:z_chunks(z+1);
                    
                    L_chunk = L(xi,yi,zi);
                    I_chunk = I(xi,yi,zi);
                    
                    pad_size = (final_size - split_size)/4;
                    L_chunk = padarray(L_chunk,pad_size,0);
                    I_chunk = padarray(I_chunk,pad_size,0);

                    
                    % Create new folder in directory
                    dir_name = fullfile(save_directory,...
                        sprintf('s%d%s',a,cen_name{1}));
                    mkdir(dir_name)
    
                    niftiwrite(L_chunk,fullfile(dir_name,'truth'),'Compressed',true)
                    niftiwrite(I_chunk,fullfile(dir_name,'t1'),'Compressed',true)
                    a = a+1;
                end
            end
        end
    else
        % Create new folder in directory
        dir_name = fullfile(save_directory,cen_name{1});
        mkdir(dir_name)
    
        niftiwrite(L,fullfile(dir_name,'truth'),'Compressed',true)
        niftiwrite(I,fullfile(dir_name,'t1'),'Compressed',true)
    end
end
else
    for n = 1:length(img_files)
    % Read image
    img_path = fullfile(img_location,img_files(n).name);
    fprintf('%s\n',img_path)

    I = read_nii(img_path);
    I = uint16(I);
    
    if size(I,1) > patch_size(1)
        I = I(1:patch_size(1),1:patch_size(2),1:patch_size(3));
    end
    
    % Normalize image
    if normalize_images
        I = normalize_image(I,remove_background);
    end
    
    % Write images
    if split_images
        n_chunks = patch_size./split_size;
        
        x_chunks = floor(linspace(0,patch_size(1),n_chunks(1)+1));
        y_chunks = floor(linspace(0,patch_size(2),n_chunks(2)+1));
        z_chunks = floor(linspace(0,patch_size(3),n_chunks(3)+1));

        a = 1;
        for x = 1:n_chunks(1)
            for y = 1:n_chunks(2)
                for z = 1:n_chunks(3)
                    xi = x_chunks(x)+1:x_chunks(x+1);
                    yi = y_chunks(y)+1:y_chunks(y+1);
                    zi = z_chunks(z)+1:z_chunks(z+1);
                    
                    I_chunk = I(xi,yi,zi);
                    
                    pad_size = (final_size - split_size)/4;
                    I_chunk = padarray(I_chunk,pad_size,0);

                    % Create new folder in directory
                    dir_name = fullfile(save_directory,...
                        sprintf('s%d%s',a,img_files(n).name));
                    mkdir(dir_name)
    
                    niftiwrite(I_chunk,fullfile(dir_name,'t1'),'Compressed',true)
                    a = a+1;
                end
            end
        end
    else
        % Create new folder in directory
        dir_name = fullfile(save_directory,img_name{1});
        mkdir(dir_name)
    
        niftiwrite(I,fullfile(dir_name,'t1'),'Compressed',true)
    end
end

end

fprintf('Preprocessed %d cells in %d images \n',n_cells,n)
fprintf('Percent of pixels: %f \n',mean(pixels)*100)

end

function I = read_nii(path)
    try
        I = niftiread(path);
    catch
        I = niftiread(sprintf('%s.gz',path));
    end
end

function img_adj = normalize_image(img, remove_background)
% Remove background
if remove_background
    se = strel('disk',20);
    for z = 1:size(img,3)
        img(:,:,z) = imtophat(img(:,:,z),se);
    end
end

% Adjust intensity
int_low = double(min(img(:)))/65535;
int_high = double(max(img(:)))/65535;
img_adj = imadjustn(img,[int_low, int_high(1)]);
        
end

function L_new = rewrite_centroid(L,I,capture_resolution,patch_resolution)
% Save centroids based in rescaled dimensions
res = capture_resolution./patch_resolution;
L_new = zeros(size(L));
cc = bwconncomp(L,6);
label_matrix = labelmatrix(cc);

rp = regionprops(label_matrix);
        
for j = 1:length(rp)
    pos = round(rp(j).Centroid.*res);
    L_new(pos(2),pos(1),pos(3)) = j;
end

L_new = expand_centroids3(label_matrix,I,rp);
end

function L_new = write_edge(L)
nhood = [0 1 0;1 1 1;0 1 0];
L_bound = zeros(size(L));
for z = 1:size(L,3)
    L_slice = imdilate(L_bound(:,:,z),nhood);
    boundary = bwboundaries(L(:,:,z));
    pos = cell2mat(boundary);
    idx = sub2ind(size(L_slice),pos(:,1),pos(:,2));
    L_slice(idx) = 1;
    L_bound(:,:,z) = L_slice;
end

L_new = L + uint8(L_bound);
end

function [L,I] = check_orientation(L,I)

z = round(linspace(1,size(L,3),5));

a = 1;
score = zeros(1,6);
for i = 1:2
   for j = 1:3
       cc = zeros(1,5);
       for k = 1:5
           if i == 2
            lmod = flip(L(:,:,z(k)),1);
           else
            lmod = L(:,:,z(k));
           end
           
           if j == 2
               lmod = imrotate(lmod,90);
           elseif j == 3
               lmod = imrotate(lmod,-90);
           end

            cc(k) = corr2(lmod,I(:,:,z(k)));
       end
       score(a) = mean(cc);
    a = a+1;
   end    
end

[~,idx] = max(score);

if any(idx == 4:6)
    L = flip(L,1);
end

if any(idx == [2,4])
   L = imrotate(L,90);
elseif any(idx == [3,6])
   L = imrotate(L,-90);
end

end

function L_new = erode_blob(L,I,size_threshold,pix_threshold)

L_new = zeros(size(L));
dim = size(L(:,:,1));

for z = 1:size(L,3)
   L_slice = L(:,:,z);
   I_slice = I(:,:,z);
   cc = bwconncomp(L_slice,4);
   rp = regionprops(cc);
   
   % Erode blobs a certain size
   if size_threshold > 0 
       for k = 1:cc.NumObjects
         L_object = zeros(dim);
         L_object(cc.PixelIdxList{k}) = 1;
        if length(cc.PixelIdxList{k}) > size_threshold
             L_object = L_object - bwperim(L_object); 
        end
        L_new(:,:,z) = L_new(:,:,z) + L_object;
       end
   else
       % Erode some percent of pixels around the centroid
      for k = 1:cc.NumObjects
         L_object = zeros(dim);
         L_object(cc.PixelIdxList{k}) = 1;
         npix = length(cc.PixelIdxList{k});
         pix_to_erode = floor(npix*pix_threshold);
         
         [x1,y1] = ind2sub(dim,cc.PixelIdxList{k});
         cen = rp(k).Centroid;
         D = sqrt(sum(([cen(2),cen(1)] - [x1,y1]).^2, 2));
         [~, idx] = sort(D,'descend');
         
         erode_idx = cc.PixelIdxList{k}(idx(1:pix_to_erode));
         L_object(erode_idx) = 0;

        L_new(:,:,z) = L_new(:,:,z) + L_object;
      end
   end
end

L_new = uint8(L_new);
end

