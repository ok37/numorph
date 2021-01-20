function [path_table_series, path_table_nii] = path_to_table(config,location,quick_load,save_table)
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
            "String input must be path specifying a single directory.")
        names = dir(paths);
        if any(arrayfun(@(s) contains(s.name,'aligned'),names))
            location = "aligned";
        elseif any(arrayfun(@(s) contains(s.name,'stitched'),names))
            location = "stitched";
        elseif all(arrayfun(@(s) contains(s.name,'.nii'),names)) && all(arrayfun(@(s) any(ismember(config.markers, s.name)),names))
            location = "resampled";
        else
            [path_table_series, path_table_nii] = munge_string(paths);
            return
        end
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
    elseif endsWith(config.img_directory,'.csv')
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
path_table = [];
if isstruct(config) && exist(fullfile(config.output_directory,'variables','path_table.mat'),'file') == 2
    var_location = fullfile(config.output_directory,'variables','path_table.mat');  
    load(var_location,'path_table')
    if ~isfield(path_table,location)
        quick_load = false;
        save_table = true;
    end
end

% Munge paths or read filename information from previously saved variable
switch location
    case 'raw'
        if ~isempty(path_table) && quick_load
            path_table_series = path_table.raw;
            if all(ismember(path_table_series.markers,config.markers))
                return
            end
        end        
        [path_table_series, path_table_nii] = munge_raw(config);
    case 'aligned'
        if ~isempty(path_table) && quick_load
            path_table_series = path_table.aligned;
            return            
        end
        path_table_series = munge_aligned(config);
    case 'stitched'
        if ~isempty(path_table) && quick_load
            path_table_series = path_table.stitched;
            return
        end
        path_table_series = munge_stitched(config);
    case 'resampled'
        path_table_series = [];
        path_table_nii = munge_resampled(config);
        return
    case 'csv'
        if ~isempty(path_table) && quick_load
            path_table_series = path_table.raw;
            return
        end   
        path_table_series = readtable(config.img_directory);
        return
    otherwise
        error("Unrecognized location selected")
end


% Save path_table for quicker loading next time
if save_table
    var_location = fullfile(config.output_directory,'variables','path_table.mat');  
    fprintf("%s\t Saving path table \n",datetime('now'))
    if isequal(location,"raw")
        path_table.tif_series = path_table_series;
        if ~isempty(path_table_nii)
            path_table.nii = path_table_nii;
        end
    end
    path_table.(location) = path_table_series;
    save(var_location,'path_table')
end

end
