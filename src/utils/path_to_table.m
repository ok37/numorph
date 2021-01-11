function paths_table_main = path_to_table(config,location,quick_load,save_table)
%--------------------------------------------------------------------------
% Convert a structure of image directories, identifies tile positions and 
% marker names by regular expressions, and outputs a final table organized 
% with relevant position and channel information for futher processing. The
% file identifiers used here are specific to those used by the ImSpector 
% software used by LaVision Biotech. Images must be 2D, single channel
% series.
%--------------------------------------------------------------------------

% Check file location if not provided
if nargin<2
    if isstring(config) || ischar(config)
        % Directory location and not config structure
        fprintf("%s\t Reading image filename from single folder using default settings \n",datetime('now'))
        paths = string(config);
        assert(isfolder(paths) && length(paths) == 1,...
            "Input must single folder in current directory")
        paths_table_main = munge_single(paths);
        return
    elseif isequal(config.img_directory,fullfile(config.output_directory,'aligned'))
        % Start from after multi-channel alignment    
        fprintf("%s\t Reading image filename information from aligned directory \n",datetime('now'))
        location = "aligned";
    elseif isequal(config.img_directory, fullfile(config.output_directory,'stitched'))
        % Start from after stitching
        fprintf("%s\t Reading image filename information from stitched directory \n",datetime('now'))
        location = "stitched";
    elseif isequal(config.img_directory, fullfile(config.output_directory,'resampled'))
        % Start from after resampled
        fprintf("%s\t Reading image filename information from resampled directory \n",datetime('now'))
        location = "resampled";
    elseif contains(config.img_directory,'.csv')
        % Read from csv file
        fprintf("%s\t Reading image filename information from .csv file \n",datetime('now'))
        location = "csv";
    else
        %Start from raw image directory
        fprintf("%s\t Reading image filename information from raw image directory \n",datetime('now'))
        location = "raw";
    end
elseif isempty(location)
    %Start from raw image directory by default
    fprintf("%s\t Reading image filename information from raw image directory \n",datetime('now'))
    location = "raw";
end

% Quick load from saved path_table variable
if nargin<3
    quick_load = true;
end

% Overwrite saved path_table variable
if nargin<4
    save_table = false;
end

% Load table if variable exists
var_location = fullfile(config.output_directory,'variables','path_table.mat');
if exist(var_location,'file') == 2
    load(var_location,'path_table')
    if ~isfield(path_table,location)
        quick_load = false;
        save_table = true;
    end
else
    path_table = [];
end

% Unpack paths
if isequal(location,"resampled")
    paths = {dir(fullfile(config.output_directory,'resampled'))};
    quick_load = false; 
elseif ~isequal(location,"csv") && ~isequal(location,"single")
    for i = 1:length(config.img_directory)
        paths{i} = dir(config.img_directory(i));
        % Check for empty paths fields
        if isempty(paths{i})
            error("Directory %s does not exist",config.img_directory(i))
        end
    end
end

% Munge paths or read filename information from previously saved variable
switch location
    case 'aligned'
        if ~isempty(path_table) && quick_load
            paths_table_main = path_table.aligned;
            return            
        end
        paths_new = munge_aligned(paths);
    case 'stitched'
        if ~isempty(path_table) && quick_load
            paths_table_main = path_table.stitched;
            return
        end
        paths_new = munge_stitched(paths);
    case 'resampled'
        paths_new = munge_resampled(paths);
        paths_table_main = struct2table(reshape([paths_new{1}],[],1));
        return
    case 'raw'
        if ~isempty(path_table) && quick_load
            paths_table_main = path_table.raw;
            if all(ismember(paths_table_main.markers,config.markers))
                return
            end
        end        
        paths_new = munge_raw(config,paths);
    case 'csv'
        if ~isempty(path_table) && quick_load
            paths_table_main = path_table.raw;
            return
        end   
        paths_table_main = readtable(config.img_directory);
        paths_new = [];
    otherwise
        error("Unrecognized location selected")
end

% Convert to table
if ~isempty(paths_new)
    for i = 1:length(paths_new)
        paths_table = struct2table(reshape([paths_new{i}],[],1));

        % Subset positions and channel numbers to base 1 index
        if iscell(paths_table.x)
            paths_table.x = [paths_table.x{:}] + 1-min([paths_table.x{:}]);
            paths_table.y = [paths_table.y{:}] + 1-min([paths_table.y{:}]);
            paths_table.z = [paths_table.z{:}] + 1-min([paths_table.z{:}]);
        else
            paths_table.x = paths_table.x + 1-min(paths_table.x);
            paths_table.y = paths_table.y + 1-min(paths_table.y);
            paths_table.z = paths_table.z + 1-min(paths_table.z);
        end
        
        if i > 1
            paths_table_main = cat(1,paths_table_main,paths_table);
        else
            paths_table_main = paths_table;
        end
    end
end

% Subset positions and channel numbers to base 1 index
paths_table_main.x = paths_table_main.x + 1-min(paths_table_main.x);
paths_table_main.y = paths_table_main.y + 1-min(paths_table_main.y);
paths_table_main.z = paths_table_main.z + 1-min(paths_table_main.z);
paths_table_main.channel_num = paths_table_main.channel_num + 1-min(paths_table_main.channel_num);


% Check table for correct size
if height(unique(paths_table_main(:,2:end),'rows')) ~= height(paths_table_main)
    error("Duplicate entries found in assmebled image file table")
end

%total_images = max(paths_table_main.x)*max(paths_table_main.y)*max(paths_table_main.z)*max(paths_table_main.channel_num);
%if total_images > height(paths_table_main)
%    error("Extra entries detected for %d tiles in assembled image file table",...
%        max(unique(paths_table_main.x))*max(unique(paths_table_main.y)))
%elseif total_images < height(paths_table_main)
%    error("Missing entries detected for %d tiles in assembled image file table",...
%        max(unique(paths_table_main.x))*max(unique(paths_table_main.y)))
%end

% Check for multiple channels in single image
if length(imfinfo(paths_table_main.file{1})) > 1
    error("Multi-page .tif detected. NuMorph currently does not support multi-channel .tif images")
end

% Save path_table for quicker loading next time
if save_table
    fprintf("%s\t Saving path table \n",datetime('now'))
    path_table.(location) = paths_table_main;
    save(var_location,'path_table')
end

end

function paths_new = munge_single(paths)
% Try to create path table for single folder with whatever information is
% available.

position_exp = ["[\d*", "\d*]","Z\d*"]; % Default position expression vector

paths = dir(paths);
paths = paths(arrayfun(@(s) contains(s.name, '.tif'),paths));
files = {paths.name};

% Fill in table
paths_new = table('Size',[length(files), 7],...
    'VariableTypes',{'cell','string','string','string','double','double','double'},...
    'VariableNames',{'file','sample_name','markers','channel_num','x','y','z'});
paths_new.file = arrayfun(@(s,t) fullfile(s.folder,t),paths,files');

% Can only do channel number here
paths_new.channel_num = cellfun(@(s) str2double(regexprep(regexp(s,'_C\d*_','match'),{'_C','_'},'')),files)';

% Do slice positions
paths_new.y = cellfun(@(s) str2double(regexprep(regexp(s,position_exp(1),'match'),'[^\d+]','')),files)';
paths_new.x = cellfun(@(s) str2double(regexprep(regexp(s,position_exp(2),'match'),'[^\d+]','')),files)';
paths_new.z = cellfun(@(s) str2double(regexprep(regexp(s,position_exp(3),'match'),'[^\d+]','')),files)';

% Set channels and positions to start at 1
paths_new.channel_num = paths_new.channel_num - min(paths_new.channel_num) + 1;
paths_new.x = paths_new.x- min(paths_new.x) + 1;
paths_new.y = paths_new.y - min(paths_new.y) + 1;
paths_new.z = paths_new.z - min(paths_new.z) + 1;

% Fill in remainig columns
paths_new.sample_name = repmat("SAMPLE",height(paths_new),1);
paths_new.markers = "marker"+num2str(paths_new.channel_num);

end


function paths_new = munge_stitched(paths)
%From stitching step  

% Remove unnecessary fields
fields_to_remove = {'folder','date','isdir','bytes','datenum'};

% Get directory contents
if ~isstruct(paths)
    paths_sub = dir(fullfile(paths{1}(1).folder));
else
    paths_sub = paths;
end

% Check .tif in current folder
paths_sub = paths_sub(arrayfun(@(x) contains(x.name,'.tif'),paths_sub));

% Create new field for file location
C = arrayfun(@(x) fullfile(paths_sub(1).folder,x.name),paths_sub,'UniformOutput',false);
[paths_sub.file] = C{:};

% Remove unused fields
paths_sub = rmfield(paths_sub,fields_to_remove);

% Generate table components
components = arrayfun(@(s) strsplit(s.name,{'_','.'}), paths_sub, 'UniformOutput', false);
components = vertcat(components{:});

% Take image information
[paths_sub.sample_name] = components{:,1};
channel_num = cellfun(@(s) str2double(s(2)),components(:,3),'UniformOutput',false);
[paths_sub.channel_num] = channel_num{:,1};
[paths_sub.markers] = components{:,4};
positions = cellfun(@(s) str2double(s),components(:,[6,5,2]),'UniformOutput',false);
x = num2cell(ones(1,length(paths_sub)));

% Set x,y tiles as 1 and save z positions
[paths_sub.x] = x{:};
[paths_sub.y] = x{:};
[paths_sub.z] = positions{:,3};

paths_new = {rmfield(paths_sub,{'name'})};

end


function paths_new = munge_aligned(paths)
%From alignment step  

% Remove unnecessary fields
fields_to_remove = {'folder','date','isdir','bytes','datenum'};

% Get directory contents
if ~isstruct(paths)
    paths_sub = dir(fullfile(paths{1}(1).folder));
else
    paths_sub = paths;
end

%Check .tif in current folder
paths_sub = paths_sub(arrayfun(@(x) contains(x.name,'.tif'),paths_sub));

%Create new field for file location
C = arrayfun(@(x) fullfile(paths_sub(1).folder,x.name),paths_sub,'UniformOutput',false);
[paths_sub.file] = C{:};

paths_sub = rmfield(paths_sub,fields_to_remove);

%Generate table components
components = arrayfun(@(s) strsplit(s.name,{'_','.'}), paths_sub, 'UniformOutput', false);
components = vertcat(components{:});

assert(size(components,1) > 0, "Could not recognize any aligned files by filename structure")
assert(length(unique(components(:,1))) == 1, "Multiple sample ids found in image path")
%assert(length(unique(components(:,4))) == length(config.markers), "Number of markers detected" +...
%    " does not match number of markers specified for this sample")

%Take image information
for i = 1:length(paths_sub)
    paths_sub(i).sample_name = string(components{i,1});
    paths_sub(i).markers = string(components{i,4});
end

positions = cellfun(@(s) str2double(s),components(:,[6,5,2]),'UniformOutput',false);
[paths_sub.x] = positions{:,1};
[paths_sub.y] = positions{:,2};
[paths_sub.z] = positions{:,3};
channel_num = cellfun(@(s) str2double(s(2)),components(:,3),'UniformOutput',false);
[paths_sub.channel_num] = channel_num{:,1};

paths_new = {rmfield(paths_sub,{'name'})};

end


function paths_new = munge_resampled(paths)
% From resampled step

% Remove unnecessary fields
fields_to_remove = {'folder','date','isdir','bytes','datenum'};

paths_sub = dir(fullfile(paths{1}(1).folder));
    
% Check .nii in current folder
paths_sub = paths_sub(arrayfun(@(x) contains(x.name,'.nii'),paths_sub));

% Create new field for file location
C = arrayfun(@(x) fullfile(paths_sub(1).folder,x.name),paths_sub,'UniformOutput',false);
[paths_sub.file] = C{:};

paths_sub = rmfield(paths_sub,fields_to_remove);

% Generate table components
components = arrayfun(@(s) strsplit(s.name,{'_','.'}), paths_sub, 'UniformOutput', false);
components = vertcat(components{:});

if ~isequal(components{1,1}, 'C1') && ~isequal(components{1,2}, 'C1')
    error("Could not recognize .nii filename structure. First file should start as SAMPLE_C1_..."+...
        ".nii")
elseif isequal(components{1,1}, 'C1')
    components = horzcat(repmat({'SAMPLE'},size(components,1),1),components);
end

% Take image information
[paths_sub.sample_name] = components{:,1};
[paths_sub.channel_num] = components{:,2};
[paths_sub.markers] = components{:,3};

positions = cellfun(@(s) str2double(s),components(:,[4,5,6]),'UniformOutput',false);
[paths_sub.x] = positions{:,1};
[paths_sub.y] = positions{:,2};
[paths_sub.z] = positions{:,3};

paths_new = {rmfield(paths_sub,{'name'})};

end


function paths_new = munge_raw(config, paths)
% Read filename information from ImSpector Software in Lavision


% Search subdirectories if no .tifs in current folder
if length(paths) == 1 && all(arrayfun(@(s) ~contains(s.name,'.tif'), paths{1}))           
    % Subset folders with images 
    paths_full = paths{1}(arrayfun(@(s) contains(s.name,cellstr(config.markers),'IgnoreCase',true),paths{1}),:);            

    % Check for multiple folders with same marker name
    for i = 1:length(config.markers)
        idx = contains({paths_full.name},config.markers(i));
        if sum(idx) > 1
            error("Please specify image directories for each channel as string array")
        end
    end
    paths_full = fullfile(paths_full(1).folder,{paths_full.name});
    paths_full = paths_full(cellfun(@(s) ~endsWith(s,{'.','..'}),paths_full));
    processed_keys = {'stitched','aligned','resampled'};

    % Check for subdirectories
    sub_dir = cellfun(@(s) isfolder(s),paths_full);
    paths = {};
    for i = 1:length(paths_full)
        if ~sub_dir(i)
            continue
        else
            paths_sub = dir(paths_full{i});
            if ~any(arrayfun(@(s) contains(s.name,processed_keys),paths_sub))
                paths = cat(2,paths,paths_sub);
            end
        end
    end
end

% Use regular expressions to extract information
if length(config.img_directory) == length(config.markers)
    % Unique directory for each channel
    if ~isfield(config,'markers')
        error("Must provide marker names")
    else
        markers = config.markers;
        new_channels = 1:length(markers);
    end

    for i = 1:length(paths)
        % Take only tif files
        path_idx = paths{i}(arrayfun(@(x) contains(x.name,'.tif'),paths{i}));

        % Subset only this marker
        path_idx = path_idx(arrayfun(@(x) contains(x.name,markers(i)),path_idx));

        % Subset only this channel
        if isfield(config,'channel_num')
            path_idx = path_idx(arrayfun(@(x) contains(x.name,config.channel_num(i)),path_idx));
        end

        % Save full file path
        paths_sub = cell2struct(fullfile({path_idx.folder},{path_idx.name}),'file');

        % Scan through each filename for relevant information
        for j = 1:length(path_idx)
            try
                paths_sub(j).sample_name = config.sample_name;
                paths_sub(j).y = str2double(regexprep(regexp(path_idx(j).name,config.position_exp(1),'match'),'[^\d+]',''));
                paths_sub(j).x = str2double(regexprep(regexp(path_idx(j).name,config.position_exp(2),'match'),'[^\d+]',''));
                paths_sub(j).z = str2double(regexprep(regexp(path_idx(j).name,config.position_exp(3),'match'),'[^\d+]',''));

                if isempty([paths_sub(j).x]) || isempty([paths_sub(j).y]) || isempty([paths_sub(j).z])
                    if contains(path_idx(j).name,'stitched.tif')
                        paths_sub = munge_stitched(path_idx);
                        break
                    elseif contains(path_idx(j).name,'aligned.tif')
                        paths_sub = munge_aligned(path_idx);
                        break
                    else
                        error("Error reading file: %s",paths_sub(j).file)
                    end
                end
                paths_sub(j).channel_num = new_channels(i);
                paths_sub(j).markers = markers(i);
            catch ME
                error("Error reading file: %s",paths_sub(j).file)
            end
        end

        % Save into cell array
        if iscell(paths_sub)
            paths_new(i) = paths_sub;
        else
            paths_new{i} = paths_sub;
        end

    end
    % Check to see if equal number of z positions
    %%%%%%% Change in future update where z resolution is not
    %%%%%%% equal
    if ~all(cellfun(@(s) length(s),paths_new))
        error("Unequal z positions")
    end
else
    % Directories contain multiple channels
    if ~isfield(config,'channel_num')
        error("Channel number indicators must be specified to read these images correctly. "+...
            "Otherwise specify a unique directory for each marker.")
    end
    if length(config.channel_num) ~= length(config.markers)
        error("Number of markers must equal the number of channel numbers specified.")
    else
        markers = config.markers;
        new_channels = 1:length(markers);
    end

    % Combine paths
    for i = 2:length(paths)
        paths{1} = cat(1,paths{1},paths{i});
    end

    % Take only tif files
    paths = paths{1}(arrayfun(@(x) contains(x.name,'.tif'),paths{1}));
    if isempty(paths)
        error("No .tif files detected in input image directory")
    end

    % Save full file path
    paths_sub = cell2struct(fullfile({paths.folder},{paths.name}),'file');

    % Scan through each filename for relevant information
    for i = 1:length(paths_sub)
        try
            paths_sub(i).sample_name = config.sample_name;
            paths_sub(i).y = str2double(regexprep(regexp(paths(i).name,config.position_exp(1),'match'),'[^\d+]',''));
            paths_sub(i).x = str2double(regexprep(regexp(paths(i).name,config.position_exp(2),'match'),'[^\d+]',''));
            paths_sub(i).z = str2double(regexprep(regexp(paths(i).name,config.position_exp(3),'match'),'[^\d+]',''));
            channel_idx = cellfun(@(s) ~isempty(s),regexp(paths(i).name,config.channel_num,'match'));
            if sum(channel_idx) > 1
                % If non-unique channel number identifiers, parse by
                % marker
                channel_idx = channel_idx & cellfun(@(s) ~isempty(s),regexp(paths(i).name,markers,'match'));
            end
            paths_sub(i).channel_num = new_channels(channel_idx);
            paths_sub(i).markers = markers(channel_idx);
        catch ME
            if iscell(paths_sub(i))
                error("Error reading file: %s\n\n%s",paths_sub{i}.file,ME.message)
            else
                error("Error reading file: %s\n\n%s",paths_sub(i).file,ME.message)
            end
        end
    end
    paths_new{1} = paths_sub;
end

end