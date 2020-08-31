function path_table_main = path_to_table(path,location,markers,channel_num)
%--------------------------------------------------------------------------
% Convert a structure of image directories, identifies tile positions and 
% marker names by regular expressions, and outputs a final table organized 
% with relevant position and channel information for futher processing. The
% file identifiers used here are specific to those used by the ImSpector 
% software used by LaVision Biotech. Images must be 2D, single channel
% series.
%--------------------------------------------------------------------------

% Check for empty path fields
if any(cellfun(@(x) isempty(x),path))
    error("Empty path field detected. Check indicated channels and path names for corrrect locations")
end

% Remove unnecessary fields
fields_to_remove = {'folder','date','isdir','bytes','datenum'};
path_new = cell(1,length(markers));
a=1;

% Read filename information from previously done processing steps
% From alignment step
switch location
    case 'aligned'
        path_sub = dir(fullfile(path{1}(1).folder));

        %Check .tif in current folder
        path_sub = path_sub(arrayfun(@(x) contains(x.name,'.tif'),path_sub));

        %Create new field for file location
        C = arrayfun(@(x) fullfile(path_sub(1).folder,x.name),path_sub,'UniformOutput',false);
        [path_sub.file] = C{:};

        path_sub = rmfield(path_sub,fields_to_remove);

        %Generate table components
        components = arrayfun(@(s) strsplit(s.name,{'_','.'}), path_sub, 'UniformOutput', false);
        components = vertcat(components{:});

        %Take image information
        [path_sub.sample_name] = components{:,1};
        [path_sub.channel_num] = components{:,3};
        [path_sub.markers] = components{:,4};
        positions = cellfun(@(s) str2double(s),components(:,[6,5,2]),'UniformOutput',false);
        [path_sub.x] = positions{:,1};
        [path_sub.y] = positions{:,2};
        [path_sub.z] = positions{:,3};

        path_new{1} = rmfield(path_sub,{'name'});

    case 'stitched'
        %From stitching step  
        path_sub = dir(fullfile(path{1}(1).folder));

        %Check .tif in current folder
        path_sub = path_sub(arrayfun(@(x) contains(x.name,'.tif'),path_sub));

        %Create new field for file location
        C = arrayfun(@(x) fullfile(path_sub(1).folder,x.name),path_sub,'UniformOutput',false);
        [path_sub.file] = C{:};

        path_sub = rmfield(path_sub,fields_to_remove);

        %Generate table components
        components = arrayfun(@(s) strsplit(s.name,{'_','.'}), path_sub, 'UniformOutput', false);
        components = vertcat(components{:});

        %Take image information
        [path_sub.sample_name] = components{:,1};
        [path_sub.channel_num] = components{:,3};
        [path_sub.markers] = components{:,4};
        positions = cellfun(@(s) str2num(s),components(:,[6,5,2]),'UniformOutput',false);
        x = num2cell(ones(1,length(path_sub)));
        [path_sub.x] = x{:};
        [path_sub.y] = x{:};
        [path_sub.z] = positions{:,3};

        path_new{1} = rmfield(path_sub,{'name'});
        
    case 'resampled'
        path_sub = dir(fullfile(path{1}(1).folder));
    
        % Check .nii in current folder
        path_sub = path_sub(arrayfun(@(x) contains(x.name,'.nii'),path_sub));
        
        if ~isequal(length(path_sub), length(markers))
            error("Number of .nii files does not match number of specified markers");
        end
    
        % Create new field for file location
        C = arrayfun(@(x) fullfile(path_sub(1).folder,x.name),path_sub,'UniformOutput',false);
        [path_sub.file] = C{:};
    
        path_sub = rmfield(path_sub,fields_to_remove);
    
        %Generate table components
        components = arrayfun(@(s) strsplit(s.name,{'_','.'}), path_sub, 'UniformOutput', false);
        components = vertcat(components{:});
        
        %Take image information
        [path_sub.sample_name] = components{:,1};
        [path_sub.channel_num] = components{:,2};
        [path_sub.markers] = components{:,3};
    
        path_new{1} = rmfield(path_sub,{'name'});

    case 'raw'
        % Read filename information from ImSpector Software in Lavision
        % For each folder in path
        for i = 1:length(path)
            %Remove hidden directories
            path_ = path{i}(~ismember({path{i}.name},{'.','..'}));

            %Determine if current path contains multiple channels
            path_cell = cellstr(char(path_(i).name));

            %Check channel names
            fun = @(s)contains(path_cell,s);
            out = cellfun(fun, cellstr(markers));
            %Make sure directories have indicated markers
            if ~out
                error('Indicated markers not found in directory %d',i)
            end

            %Check .tif in current folder
            path_sub = dir(fullfile(path{i}(1).folder));
            path_sub = path_sub(arrayfun(@(x) contains(x.name,'.tif'),path_sub));
            path_sub = path_sub(arrayfun(@(x) contains(x.name,channel_num),path_sub));

            %Create new field for file location
            C = arrayfun(@(x) fullfile(path_sub(1).folder,x.name),path_sub,'UniformOutput',false);
            [path_sub.file] = C{:};

            %Check that files contain channel numbers
            out2 = arrayfun(@(x) contains(x.name,channel_num),path_sub);
            if any(out2 ~= 1)
                error('Indicated channel numbers not found in directory %d',i)
            end

            %Remove unnecessary fields and preallocate fields
            path_sub = rmfield(path_sub,fields_to_remove);

            %Create new entries for x,y,z positions, markers values
            x_ind = strfind(path_sub(1).name,']'); %x tile
            y_ind = strfind(path_sub(1).name,'['); %x tile
            z_ind = strfind(path_sub(1).name,' Z'); % z tile

            markers_sub = markers(out);
            channels_sub = channel_num(out);
            channel_idx = find(out);

            for j = 1:length(channels_sub)
                path_sub2 = path_sub(arrayfun(@(x) contains(x.name,channels_sub(j)),path_sub));
                for z = 1:length(path_sub2)                
                    path_sub2(z).channel_num = channel_idx(j);

                    path_sub2(z).x = str2double(path_sub2(z).name(x_ind-2:x_ind-1))+1;

                    path_sub2(z).y= str2double(path_sub2(z).name(y_ind+1:y_ind+2))+1;

                    path_sub2(z).z = str2double(path_sub2(z).name(z_ind+2:z_ind+5))+1;

                    path_sub2(z).markers = markers_sub(j);

                    path_sub2(z).channel_num = channel_idx(j);
                end
                path_sub2 = rmfield(path_sub2,{'name'});
                path_new{a} = path_sub2;
                a=a+1;
            end
        end

        %Make sure correct markers are present again
        for i = 1:length(markers)
            if ~any(strcmp([path_new{i}.markers],markers(i)))
                error("Marker %s not found in directory %s",markers(i),path{i}.folder)
            end
        end
end

%Convert to table
path_table_main = struct2table(reshape([path_new{:}],[],1));

end