function [S, B] = estimate_flatfield(path_table, config)

sampling_freq = 0.2;

if ~isempty(config.shading_correction_tiles)
    max_x = max(unique(path_table.x));
    max_y = max(unique(path_table.y));
    
    [x,y] = ind2sub([max_x, max_y],config.shading_correction_tiles);
    
    path_table = path_table(ismember(path_table.x,x),:);
    path_table = path_table(ismember(path_table.y,y),:);
end

x_tiles = unique(path_table.x);
y_tiles = unique(path_table.y);

disp(x_tiles)
disp(y_tiles)

n_tiles = length(x_tiles)*length(y_tiles);
n_samples = round(sampling_freq*height(path_table)/n_tiles); 

z_idx = floor(linspace(1,max(path_table.z),n_samples));
path_sub = path_table(ismember(path_table.z,z_idx),:);
n_images = height(path_sub);

fprintf('%s\t Loading %d images for estimating flatfield \n',datetime('now'),n_images)

tempI = loadtiff(path_sub.file{1});
I = zeros(size(tempI),'single');

tic
parfor i = 1:n_images
    I(:,:,i) = loadtiff(path_sub.file{i});
    if i == round(n_images/2)
        fprintf('%s\t Halfway there... \n',datetime('now'))
    end
end
toc

fprintf('%s\t Running BaSIC \n',datetime('now'),height(path_sub))

[S,B] = BaSiC(I,'darkfield','true'); 

end