function I_final = plot_centroid_stack(config, spacing, plot_classes, save_image)
%--------------------------------------------------------------------------
% Visualize image centroids on a corresponding z slice.
%--------------------------------------------------------------------------
% Usage:
% I_final = plot_centroid_stack(config, spacing, plot_classes, save_image)
%
%--------------------------------------------------------------------------
% Inputs:
% config: Analyis configuration structure. 
%
% spacing: (1x3 integer) Amount of downsampling for each dimension. 
% (default: [10,10,100])
% 
% plot_classes: (logical) Color by classificaiton if they exist. 
% (default: true)
% 
% save_image: (logical) Save image to 'samples' directory. (default: true,
% if no output, false otherwise)
% 
%--------------------------------------------------------------------------
% Outputs:
% I_final: (mxnx2 uint) Contains raw image in slice 1 and associated
% centroid calls in slice 2. To create overlayed image run:
% img = labeloverlay(I_final(:,:,1),I_final(:,:,2),'ColorMap','prism');
%
%--------------------------------------------------------------------------

if nargin<3; plot_classes = true; end
if nargin<4
    if nargout == 1
        save_image = false;
    else
        save_image = true;
    end
end

% Check for centroids
res_path = fullfile(config.output_directory,strcat(config.sample_id,'_results.mat'));
if ~isfile(res_path)
    error("Could not locate results structure in %s output directory",config.sample_id)
else
    var_names = who('-file',res_path);
    if ismember('centroids',var_names)
        load(res_path,'centroids')
    else
        error("No centroids detected in results structure")
    end
    if plot_classes
        if ismember('classes',res_path)
            load(res_path,'classes')
            assert(length(classes) == size(centroids,1),...
                "Number of centroids and classes do not match")
        else
            warning("No classes detected in results structure")
            plot_classes = false;
        end
    end
end
    
% Get image paths
path_table = path_to_table(config);
path_table = path_table(path_table.channel_num == 1,:);

% Get z position at current slice and slices above/below
z_position = z_position-1:z_position+1;
idx = ismember(centroids(:,3),z_position);
cen = centroids(idx,:);
if plot_classes
    class = classes(idx,:);
else
    class = ones(1,size(cen,1));
end
I = read_img(path_table,[1,z_position(2)]);

% Create label image
L = zeros(size(I));
se = strel([0 1 0;1 1 1;0 1 0]);
for i = 1:3
    img = zeros(size(I));
    idx = cen(:,3) == z_position(i);
    val = class(idx);
    idx = sub2ind(size(I),cen(idx,1),cen(idx,2));
    img(idx) = val;
    if i == 2
        img = imdilate(img,se);
    end
    L = L + img;
end

% Save image
if save_image
    I = imadjust(I);
    save_img = labeloverlay(I,L,'ColorMap','prism');
    save_name = sprintf('%s_slice_%d.png',config.sample_id,z_position(2));
    sample_dir = fullfile(config.output_directory,'samples');
    if ~isfolder(sample_dir)
        mkdir(sample_dir)
    end
    save_name = fullfile(sample_dir,save_name);
    fprintf('%s\t Writing image %s \n',datetime('now'),save_name)
    imwrite(save_img,save_name)
end

% Create nxmx2 image
if nargout == 1
    I_final = cat(3,I,L);
end

end