function I_final = create_image_slices2(config, path_table, filename, save_flag)
% Create overlay of centroids with cell-type

classes = 1:3;      % Classes to display. Note can only display 3
spacing = 0.1;      % 
resampling = 0.1;

% Check config stage
if ~isequal(config.stage,"analyze")
    error("Analysis configuration structure required")
end

% Get centroids filename
if nargin<3 || isempty(filename)
    % Use classes if present
    files = dir(fullfile(config.output_directory,"*.csv"));
    names = {files.name};
    if any(contains(names,'classes'))
        files = files(contains(names,'classes'));
        s1 = arrayfun(@(s) strsplit(s.name,{'classes','.csv'}),files,...
            'UniformOutput',false);
        use_classes = true;
        
    elseif any(contains(names,'centroids'))
        files = files(contains(names,'centroids'));
        s1 = arrayfun(@(s) strsplit(s.name,{'centroids','.csv'}),files,...
            'UniformOutput',false);
        use_classes = false;
    else
        error("Could not locate centroids or classes file in the output directory")
    end
    
    s1 = cellfun(@(s) str2double(s(2)),s1);
    files = files(s1 == max(s1));
    filename = fullfile(files(1).folder,files(1).name);
else
    [~,file] = fileparts(filename);
    if contains(file,'classes')
        use_classes = true;
    else
        use_classes = false;
    end
end

if nargin<4
    save_flag = false;
end

% Read centroids
centroids = readmatrix(filename);

% Get z positions
min_z = min(path_table.z);
max_z = max(path_table.z);

z_positions = length(unique(path_table.z));
z_positions = ceil(z_positions*spacing);
z_positions = linspace(min_z,max_z,z_positions);

% Subset centroids
centroids(:,1:3) = centroids(:,1:3);
idx = any(centroids(:,3) == z_positions,2);
centroids = centroids(idx,:);

% Get position and value
if use_classes
    centroids = centroids(ismember(centroids(:,end),classes),:);
    val = centroids(:,end);
else
    val = ones(size(centroids,1),1)+1;
end
pos = centroids(:,1:3);

% Get dimensions
tempI = imread(path_table.file{1});
dims = round(size(tempI)*resampling);

I_final = zeros([dims(1), dims(2), 3, length(z_positions)],'uint8');
se = strel('disk',10);

for z = 1:length(z_positions)
    R = zeros(size(tempI),'uint8');
    G = zeros(size(tempI),'uint8');
    B = zeros(size(tempI),'uint8');
    
    
    pos1 = pos(pos(:,3) == z_positions(z),:);
    val1 = val(pos(:,3) == z_positions(z),:);
    
    for i = 1:size(pos1,1)
        if val1(i) == 2
            R(pos1(i,1),pos1(i,2)) = 254;
        end
        if val1(i) == 3
            G(pos1(i,1),pos1(i,2)) = 254;
        end
        if val1(i) == 4
            B(pos1(i,1),pos1(i,2)) = 254;
        end
    end
    R = imdilate(R,se);
    G = imdilate(G,se);
    B = imdilate(B,se);
    
    R = imresize(R,[dims(1) dims(2)]);
    G = imresize(G,[dims(1) dims(2)]);
    B = imresize(B,[dims(1) dims(2)]);
        
    I_final(:,:,:,z) = cat(3,R,G,B);
end

if save_flag
    fprintf('%s\t Writing .tif image \n',datetime('now'))
    save_name = sprintf('%s_%d_%d.tif',config.sample_id,spacing);
    options.color = true;
    options.overwrite = true;
    options.verbose = false;
    saveastiff(I_final,save_name,options);
end

end