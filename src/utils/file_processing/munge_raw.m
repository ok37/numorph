function [table_series_final, path_table_nii] = munge_raw(config)
%--------------------------------------------------------------------------
% Read filename information from raw image directories. Can contain .tif
% files from light-sheet or multiple .nii file formats for MRI images. 
%
% These are all possible image extensions we're able to read:
% ext = [".tif",".tiff",".nii",".nrrd",".nhdr",".mhd"];
%--------------------------------------------------------------------------
% Usage:
% [path_table_series, path_table_nii] = munge_raw(config);
%
% Inputs:
% config: Configuration structure generated by NM_config using any stage.
%
% Outputs:
% path_table_series: Table with munged .tif series.
%
% path_table_nii: Table with munged .nii files.
%--------------------------------------------------------------------------

% These are all possible image extensions we're able to read
ext = [".tif",".tiff",".nii",".nrrd",".nhdr",".mhd"];

% Get list of all files in all sub_directories
if length(config.img_directory) == 1
    all_files = dir(fullfile(config.img_directory,'**/*.*'));
else
    assert(length(config.img_directory)<=length(config.markers),"More image "+...
        "directories specified than there are markers")
    all_files = cell(1,length(config.img_directory));
    for i = 1:length(config.img_directory)
        all_files{i} = dir(fullfile(config.img_directory(i),'*.*'));
    end
    all_files = cat(1,all_files{:});
    all_files = all_files(arrayfun(@(s) contains(s.name,ext),all_files),:);
end

% Remove output directory folder
if config.img_directory == config.output_directory
    error("Input image directory cannot match output directory")
else
    all_files = all_files(arrayfun(@(s) ~contains(s.folder,config.output_directory,'IgnoreCase',true),all_files));
end

% Get which image extensions are present and subset these images
idx = false(length(ext),length(all_files));
for i = 1:length(ext)
    idx(i,:) = arrayfun(@(s) contains(s.name,ext(i)),all_files);
end
ext_here = ext(any(idx,2));

% Read nifiti file formats if present
nifti_ext = ext_here(~ismember(ext_here,[".tif",".tiff"]));
if ~isempty(nifti_ext) && isfield(config,'mri_markers')
    nifti_files = all_files(any(idx(3:end,:)));
    path_table_nii = munge_nifti_raw(nifti_files,config);
else
    path_table_nii = [];
end

% Subset tif files
% Divide into different folders
tiff_files = all_files(any(idx(1:2,:)));

% Return if no .tif files present
if isempty(tiff_files)
    warning("No raw .tif files detected in image directories") 
    path_table_series = [];
    return
end

% If specified, remove channels without channel number
if isfield(config,'channel_num') && ~isempty(config.channel_num)
    idx = false(length(tiff_files),1);
    c_idx = cell(1,length(config.channel_num));
    for i = 1:length(config.channel_num)
        a = arrayfun(@(s) contains(s.name,config.channel_num(i)),tiff_files);
        b = a & arrayfun(@(s) contains(s.name,config.markers(i)),tiff_files);
        if ~isempty(b) && any(b)
            c_idx{i} = b;
            idx = idx | b;
        else
            c_idx{i} = a;
            idx = idx | a;
        end
    end
    tiff_files = tiff_files(idx);
    c_idx = cellfun(@(s) s(idx),c_idx,'UniformOutput',false);
end

tiff_folders = unique({tiff_files.folder});
assert(length(tiff_folders) <= length(config.markers), "More .tif folders than image markers detected.")

% Try matching folders by channel number, then by folder name, then by
% marker name in file
path_table_series = cell(1,length(tiff_folders));
files_read = false;
if length(tiff_folders) == length(config.markers)
    % Unique marker per folder
    for i = 1:length(tiff_folders)
        path_idx = [];
        % Check channel num
        if isfield(config,'channel_num') && ~isempty(config.channel_num)
            path_idx = tiff_files(c_idx{i});
        end
        % Check folder name
        if isempty(path_idx)
            path_idx = tiff_files(arrayfun(@(s) contains(s.folder,config.markers(i)),tiff_files));
        end
        
        % Check for unique marker names in filename
        if isempty(path_idx)
           path_idx = tiff_files(arrayfun(@(s) contains(s.name,config.markers(i)),tiff_files));
        end
        
        % Finally, just assign this folder to marker based on index
        if isempty(path_idx)
            path_idx = tiff_files(arrayfun(@(s) isequal(s.folder,tiff_folders{i}),tiff_files),:);
        end

        % Use regular expressions to extract information and save 
        path_table_series{i} = get_table_struct(path_idx, config.sample_id, config.markers(i), i, config.position_exp);
    end
    files_read = true;
end

% More complicated scenario
% Different number of folders and markers present
if ~files_read && isfield(config,'channel_num') && ~isempty(config.channel_num)
    % Channel numbers provided to read multiple channels in the same folder
    for i = 1:length(c_idx)
        path_idx = tiff_files(c_idx{i});
        
        % Use regular expressions to extract information and save 
        path_table_series{i} = get_table_struct(path_idx, config.sample_id, config.markers(i), i, config.position_exp);
    end    
elseif ~files_read
    % Finally check if each folder and/or filename has unique marker, 
    % otherwise give error
    
    % Check files for marker name in the filename
    c_idx = zeros(length(config.markers),length(tiff_files));
    f_idx = zeros(1,length(config.markers));
    for i = 1:length(config.markers)
        c_idx(i,:) = arrayfun(@(s) contains(s.name,config.markers(i)),tiff_files);
        f_idx(i) = f_idx(i) + cellfun(@(s) contains(s,config.markers(i)),tiff_folders);
    end
    
    if ~any(sum(c_idx) > 1)
        % If unqiue marker present in each filename
        for i = 1:length(config.markers)
            path_idx = tiff_files(logical(c_idx(i,:)));
            
            % Use regular expressions to extract information and save 
            path_table_series{i} = get_table_struct(path_idx, config.sample_id, config.markers(i), i, config.position_exp);
        end
    elseif all(f_idx == ones(1,length(config.markers)))
        for i = 1:length(config.markers)
            f = cellfun(@(s) contains(s,config.markers(i)),tiff_folders);
            path_idx = tiff_files(arrayfun(@(s) isequal(s.folder,tiff_folders{f}),tiff_files));
            
            % Use regular expressions to extract information and save 
            path_table_series{i} = get_table_struct(path_idx, config.sample_id, config.markers(i), i, config.position_exp);
        end
    else
        error("Channel number indexes or unique folders for each marker are needed to import these images")
    end
end

% Create table for each marker
for i = 1:length(path_table_series)
   table_series = struct2table(path_table_series{i});

    % Set x,y,z positions to start at 1
    if ~iscell(table_series.y)
        table_series.y = table_series.y - min(table_series.y)+1;
    else
        table_series.y = ones(height(table_series),1);
    end
    if ~iscell(table_series.x)
        table_series.x = table_series.x - min(table_series.x)+1;
    else
        table_series.x = ones(height(table_series),1);
    end
        
    table_series.z = table_series.z - min(table_series.z)+1; 
    
    table_series_final.(config.markers(i)) = table_series;
end

end


function paths_nifti = munge_nifti_raw(nifti_files, config)
% Create structure for nifiti files. Assumed to single tile, single
% channel, multi-slice. Must contain a marker in the filename

file = cell(1,length(nifti_files));
marker = repmat("",1,length(nifti_files));
channel_num = zeros(1,length(nifti_files));
y_res = zeros(1,length(nifti_files));
x_res = zeros(1,length(nifti_files));
z_res = zeros(1,length(nifti_files));

% Subset only images with marker present
idx=1;
for i = 1:length(config.mri_markers)
    sub = nifti_files(arrayfun(@(s) contains(s.name,config.mri_markers(i)),nifti_files));
    for j = 1:length(sub)
       file{idx} = fullfile(sub(j).folder,sub(j).name);
       marker(idx) = config.mri_markers(i);
       channel_num(idx) = i;
       y_res(idx) = config.mri_resolution(1);
       x_res(idx) = config.mri_resolution(2);
       z_res(idx) = config.mri_resolution(3);
       idx = idx+1;
    end
end

% Create table
paths_nifti = table('Size',[length(1:idx-1), 7],...
    'VariableTypes',{'cell','string','string','string','double','double','double'},...
    'VariableNames',{'file','sample_id','marker','channel_num','x_res','y_res','z_res'});

paths_nifti.file = file(1:idx-1)';
paths_nifti.sample_id = repmat(config.sample_id,idx-1,1);
paths_nifti.marker = marker(1:idx-1)';
paths_nifti.channel_num = channel_num(1:idx-1)';
paths_nifti.y_res = y_res(1:idx-1)';
paths_nifti.x_res = x_res(1:idx-1)';
paths_nifti.z_res = z_res(1:idx-1)';

end


function path_sub = get_table_struct(path_idx, sample_id, marker, channel, position_exp)

% Save full file path
path_sub = cell2struct(fullfile({path_idx.folder},{path_idx.name}),'file');

% Scan through each filename for relevant information
for j = 1:length(path_idx)
    try
        path_sub(j).sample_id = sample_id;
        if position_exp(1) ~= ""
            path_sub(j).y = str2double(regexprep(regexp(path_idx(j).name,position_exp(1),'match'),'[^\d+]',''));
        else
            path_sub(j).y =1;
        end
        if position_exp(2) ~= ""
            path_sub(j).x = str2double(regexprep(regexp(path_idx(j).name,position_exp(2),'match'),'[^\d+]',''));
        else
            path_sub(j).x = 1;
        end
        if position_exp(3) ~= ""
            path_sub(j).z = str2double(regexprep(regexp(path_idx(j).name,position_exp(3),'match'),'[^\d+]',''));
        else
            path_sub(j).z = 1;
        end

        if isempty([path_sub(j).z])
            if contains(path_sub(j).name,'stitched.tif')
                path_sub = munge_stitched(path_idx);
                break
            elseif contains(path_sub(j).name,'aligned.tif')
                path_sub = munge_aligned(path_idx);
                break
            else
                error("Error reading file: %s",path_sub(j).file)
            end
        end
        path_sub(j).channel_num = channel;
        path_sub(j).markers = marker;
    catch ME
        error("Error reading file: %s",path_sub(j).file)
    end
end

end