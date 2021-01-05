function [patches, ftable, cen_sub, cen_idx] = get_centroid_patches(centroids, path_table_stitched, config, patch_size, k, low_thresh)
%--------------------------------------------------------------------------
% Get 2D image pacthes around centroid positions
%--------------------------------------------------------------------------
% Inputs:
% centroids - matrix with centroids positions + mask annotations.
%
% path_table_stitched - input image table.
%
% config - config structure from NM_analyze.
%
% patch_size - (optional) 1x2 element vector. 1st element is 0.5xlength of
% saved patch image for viewing. 2nd element is 0.5xlength of the actual 
% sampling window. Should be less than 1st element and ideally set to the
% approximate radius of a nucleus. 
%
% k - (optional) number of patches to generate (default: 1000 patches)
%
% low_thresh - (optional) threshold for minimum intensity for all channels
% set as fraction (i.e. 0-1). For example, if set to 0.5, nuclei with
% intensities below the 50th percentile in all channels will be presumed
% negative and ignored from patch generation.
%
% Outputs:
%
%
%--------------------------------------------------------------------------

% Defaults
if nargin < 4
    p = 50;
    s = 6;
    fprintf('%s\t Using default sampling window of %d pixels\n',datetime('now'),s)
elseif length(patch_size)<2
    s = 6;
    fprintf('%s\t Using default sampling window of %d pixels\n',datetime('now'),s)
else
    p = patch_size(1);
    s = patch_size(2);
end

if nargin < 5
    k = 1000;
end

if nargin < 6
    low_thresh = 0.5;
end

% Make an output directory
save_directory = fullfile(config.output_directory,'classifier');
if exist(save_directory,'dir') ~= 7
    mkdir(save_directory);
end

markers = config.markers;
if length(markers)>3
    markers = markers(1:3);
end

if ~isequal(centroids,'load')
    % Subset cells above minimum intensity threshold
    k_idx = zeros(size(centroids,1),1);
    for i = 1:length(markers)-1
       idx = i + 5;
       thresh = prctile(centroids(:,idx),low_thresh*100);
       k_idx = k_idx | centroids(:,idx)>thresh;
    end
    fprintf('%s\t Retaining %d from %d that are above the threshold %f \n',...
        datetime('now'),sum(k_idx),length(k_idx),low_thresh)
    cen_sub = centroids(k_idx,:);
    cen_idx = find(k_idx)';

    % Take random sample from indexes
    s_idx = randsample(size(cen_sub,1),k);
    fprintf('%s\t Selecting %d random patches \n',datetime('now'),k)

    cen_sub = cen_sub(s_idx,:);
    cen_idx = cen_idx(s_idx);
    
    [~,i] = sort(cen_sub(:,3));
    cen_sub = cen_sub(i,:);
    cen_idx = cen_idx(i)';
else
   % Load previous patch info list
   fprintf('%s\t Loading previous patch info \n',datetime('now'))
   patch_name = fullfile(save_directory,sprintf('%s_patch_info.csv',config.sample_name));
   patch_info = readmatrix(patch_name); 
   cen_idx = patch_info(:,1);
   cen_sub = patch_info(:,2:end);
end

% Get intensity thresholds for rescaling
adj1 = zeros(3,2);
adj1(:,1) = 90/65535;
for i = 1:length(markers)
   idx = i + 4;
   adj1(i,2) = prctile(cen_sub(:,idx),95)/65535;
end

% Get layer and structure info
structures = string(bin_annotation_structures(cen_sub(:,4),'cortex_large'));
layers = string(bin_annotation_structures(cen_sub(:,4),'layers'));

patches = zeros([2*p+1,2*p+1,3,k],'uint8');
ftable = cell(1,length(markers));
for i = 1:size(cen_sub,1)
    % Get 
    z = cen_sub(i,3)+1;
    pos = cen_sub(i,1:2)+1;
    
    file = path_table_stitched(path_table_stitched.z == z,1).file;
    ranges = {[pos(1)-p,pos(1)+p], [pos(2)-p,pos(2)+p]};
    for j = 1:length(markers)
        img = imread(file{j},'PixelRegion',ranges);

        % Get features from patch
        f = measure_patch_features(img, s, true, config.markers(j));
        if i == 1
            ftable{j} = f;
        else
            ftable{j} = vertcat(ftable{j},f);
        end
        
        % Adjust patch intensity and add to stack
        idx = img>adj1(j,2);
        img = im2uint8(imadjust(img,[adj1(j,1) adj1(j,2)]));
        patches(:,:,j,i) = img;
    end
end

% Save
options.overwrite = true;
options.color = true;
img_name = fullfile(save_directory,sprintf('%s_patches.tif',config.sample_name));
saveastiff(patches,char(img_name),options)

% Save
if ~isequal(centroids,'load')
    patch_info = horzcat(cen_idx,cen_sub);
    patch_name = fullfile(save_directory,sprintf('%s_patch_info.csv',config.sample_name));
    writematrix(patch_info,patch_name)
end

% Save
feature_table = array2table([layers,structures],'VariableNames',{'Layer', 'Structure'});
for i = 1:length(ftable)
   feature_table = horzcat(feature_table,ftable{i});
end

table_name = fullfile(save_directory,sprintf('%s_patch_features.csv',config.sample_name));
writetable(feature_table,table_name)

end