function new_names = rename_tile_position(file_directory, overwrite)
%--------------------------------------------------------------------------
% Rename files based on adjusted file positions.
%--------------------------------------------------------------------------

position_exp = ["[\d*", "\d*]","Z\d*"];

if nargin <2
    overwrite = false;
end

files = dir(file_directory);
files = files(arrayfun(@(s) contains(s.name,'.tif'),files));

% Check if all file names the same length
assert(all(arrayfun(@(s) length(s.name) == length(s(1).name),files)), "Filename "+...
    "inconsistencies present")

tmp_name = files(1).name;
name_length = length(tmp_name);
num_files = length(files);
x = 1:name_length:(num_files+1)*name_length;
x = x(1:end-1);

% Get potential character indexes to replace
if ~isempty(position_exp)
    % By position exp
    for i = 1:length(position_exp)
        chars = regexp(tmp_name,position_exp(i),'start'):regexp(tmp_name,position_exp(i),'end');
        id{i} = chars(arrayfun(@(s) ~isnan(str2double(tmp_name(s))),chars));
    end
else
    % By incosistent characters
    name_array = [files.name];
    ic = [];
    for i = 1:name_length
       c = name_array(x+i-1);
       if ~all(arrayfun(@(s) s == c(1),c))
           ic = cat(1,ic,i);
       end
    end
    % Check for consecutive to mark numbers greater than 1 character long
    idx = arrayfun(@(s) str2double(s),files(1).name,'UniformOutput',false);
    idx = cellfun(@(s) isreal(s) && ~isempty(s) && ~isnan(s),idx);
    
    id = {};
    pre_i = false;
    for i = 1:name_length
        if idx(i) && ~pre_i
            id1 = i;
            pre_i = true;
        elseif idx(i) && pre_i
            id1 = [id1, i];
            pre_i = true;
        elseif ~idx(i) && pre_i
            id = cat(1,id,id1);
            pre_i = false;
        end
    end
    id = id(cellfun(@(s) any(ismember(s,ic)),id));
end

% Ask for user input
und = '  ^';
adj = zeros(1,length(id));
for i = 1:length(id)
    ic_range = id{i}(end)-10:id{i}(end)+10;
    if any(ic_range>name_length)
        ic_range = ic_range(1):name_length;
    end
    adj(i) = input(sprintf('Enter adjustment for this index?:\n\t%s\n\t\t%s\n',...
        files(1).name(ic_range),und));
end

% Apply adjustments
new_names = {files.name};
for i = 1:length(new_names)
    name = new_names(i);
    for j = 1:length(id)
        x1 = num2str(str2double(name{1}(id{j})) + adj(j));
        if length(x1) < length(id{j})
            % Pad with zero
            x1 = [repmat('0',1,length(id{j}) - length(x1)) x1];
        end
        name{1}(id{j}) = x1;
    end
    new_names(i) = name;
end

% Overwrite files
if overwrite
    for i = 1:length(files)
        a = fullfile(files(i).folder,files(i).name);
        b = [fullfile(files(i).folder,new_names{i}) '_tmp'];
        [status, msg] = movefile(a,b);
    end
    for i = 1:length(files)
        b = [fullfile(files(i).folder,new_names{i}) '_tmp'];
        c = fullfile(files(i).folder,new_names{i});
        [status, msg] = movefile(b,c);
    end
end

end