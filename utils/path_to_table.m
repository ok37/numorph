function paths_table_main = path_to_table(config,location)
%--------------------------------------------------------------------------
% Convert a structure of image directories, identifies tile positions and 
% marker names by regular expressions, and outputs a final table organized 
% with relevant position and channel information for futher processing. The
% file identifiers used here are specific to those used by the ImSpector 
% software used by LaVision Biotech. Images must be 2D, single channel
% series.
%--------------------------------------------------------------------------

if contains(config.img_directory,'.csv')
    % Read from csv file
    location = 'csv';
else
    % Unpack paths
    for i = 1:length(config.img_directory)
       paths{i} = dir(config.img_directory(i));
    end

    % Check for empty paths fields
    if any(cellfun(@(x) isempty(x),paths))
        error("Empty paths field detected. Check indicated channels and paths names for corrrect locations")
    end
end

% Remove unnecessary fields
fields_to_remove = {'folder','date','isdir','bytes','datenum'};

% Read filename information from previously done processing steps
switch location
    case 'aligned'
        % From alignment step
        paths_sub = dir(fullfile(paths{1}(1).folder));

        %Check .tif in current folder
        paths_sub = paths_sub(arrayfun(@(x) contains(x.name,'.tif'),paths_sub));

        %Create new field for file location
        C = arrayfun(@(x) fullfile(paths_sub(1).folder,x.name),paths_sub,'UniformOutput',false);
        [paths_sub.file] = C{:};

        paths_sub = rmfield(paths_sub,fields_to_remove);

        %Generate table components
        components = arrayfun(@(s) strsplit(s.name,{'_','.'}), paths_sub, 'UniformOutput', false);
        components = vertcat(components{:});
        
        assert(length(unique(components(:,1))) == 1, "Multiple sample ids found in image path")
        assert(length(unique(components(:,4))) == length(config.markers), "Number of markers detected" +...
            " does not match number of markers specified for this sample")

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

    case 'stitched'
        %From stitching step  
        paths_sub = dir(fullfile(paths{1}(1).folder));

        % Check .tif in current folder
        paths_sub = paths_sub(arrayfun(@(x) contains(x.name,'.tif'),paths_sub));

        % Create new field for file location
        C = arrayfun(@(x) fullfile(paths_sub(1).folder,x.name),paths_sub,'UniformOutput',false);
        [paths_sub.file] = C{:};

        paths_sub = rmfield(paths_sub,fields_to_remove);

        % Generate table components
        components = arrayfun(@(s) strsplit(s.name,{'_','.'}), paths_sub, 'UniformOutput', false);
        components = vertcat(components{:});

        % Take image information
        [paths_sub.sample_name] = components{:,1};
        [paths_sub.channel_num] = components{:,3};
        [paths_sub.markers] = components{:,4};
        positions = cellfun(@(s) str2double(s),components(:,[6,5,2]),'UniformOutput',false);
        x = num2cell(ones(1,length(paths_sub)));
        [paths_sub.x] = x{:};
        [paths_sub.y] = x{:};
        [paths_sub.z] = positions{:,3};

        paths_new = {rmfield(paths_sub,{'name'})};
        
    case 'resampled'
        paths_sub = dir(fullfile(paths{1}(1).folder));
    
        % Check .nii in current folder
        paths_sub = paths_sub(arrayfun(@(x) contains(x.name,'.nii'),paths_sub));
        
        if ~isequal(length(paths_sub), length(markers))
            error("Number of .nii files does not match number of specified markers");
        end
    
        % Create new field for file location
        C = arrayfun(@(x) fullfile(paths_sub(1).folder,x.name),paths_sub,'UniformOutput',false);
        [paths_sub.file] = C{:};
    
        paths_sub = rmfield(paths_sub,fields_to_remove);
    
        % Generate table components
        components = arrayfun(@(s) strsplit(s.name,{'_','.'}), paths_sub, 'UniformOutput', false);
        components = vertcat(components{:});
        
        % Take image information
        [paths_sub.sample_name] = components{:,1};
        [paths_sub.channel_num] = components{:,2};
        [paths_sub.markers] = components{:,3};
    
        paths_new = {rmfield(paths_sub,{'name'})};

    case 'raw'
        % Read filename information from ImSpector Software in Lavision
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
                
                % Save full file path
                paths_sub = cell2struct(fullfile({path_idx.folder},{path_idx.name}),'file');
                
                % Scan through each filename for relevant information
                for j = 1:length(path_idx)
                    try
                        paths_sub(j).sample_name = config.sample_name;
                        paths_sub(j).y = str2double(regexprep(regexp(path_idx(j).name,config.position_exp(1),'match'),'[^\d+]',''));
                        paths_sub(j).x = str2double(regexprep(regexp(path_idx(j).name,config.position_exp(2),'match'),'[^\d+]',''));
                        paths_sub(j).z = str2double(regexprep(regexp(path_idx(j).name,config.position_exp(3),'match'),'[^\d+]',''));
                        paths_sub(j).channel_num = new_channels(i);
                        paths_sub(j).markers = markers(i);
                    catch ME
                        error("Error reading file: %s",paths_sub(j).file)
                    end
                end
                % Save into cell array
                paths_new{i} = paths_sub;
            end
        else
            % Directories contain multiple channels
            if ~isfield(config,'markers')
                error("Must provide marker names")
            elseif isfield(config,'channel_num') && length(config.channel_num) ~= length(config.markers)
                error("Number of markers must equal the number of channel numbers provided")
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
                catch
                    error("Error reading file: %s",paths_sub{i}.file)
                end
            end
            paths_new{1} = paths_sub;
        end
    case 'csv'
        % Read directly from csv
        fprintf("%s\t Reading image filename information from csv file \n",datetime('now'))
        paths_table_main = readtable(config.img_directory);
        paths_new = [];
end

% Convert to table
if ~isempty(paths_new)
    paths_table_main = struct2table(reshape([paths_new{:}],[],1));
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
total_images = max(paths_table_main.x)*max(paths_table_main.y)*max(paths_table_main.z)*max(paths_table_main.channel_num);
if total_images > height(paths_table_main)
    error("Extra entries detected for %d tiles in assembled image file table",...
        max(unique(paths_table_main.x))*max(unique(paths_table_main.y)))
elseif total_images < height(paths_table_main)
    error("Missing entries detected for %d tiles in assembled image file table",...
        max(unique(paths_table_main.x))*max(unique(paths_table_main.y)))
end

% Check if images are single channel
if length(imfinfo(paths_table_main.file{1})) > 1
    error("Multi-page .tif detected. NuMorph currently does not support multi-channel .tif images")
end

end