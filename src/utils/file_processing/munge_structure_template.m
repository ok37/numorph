function munge_structure_template(csv_path)
%--------------------------------------------------------------------------
% Munge a custom structure template to define structure annotations and
% save in a standardized format. The following variables (columns) are
% required in the input csv file:
% 
% 'index', 'acronym', 'name', 'depth', 'parent_structure_id'
% 
% Use /annotations/structure_template.csv as a starting point for defining
% these values.
%
%--------------------------------------------------------------------------
% Usage:
% munge_structure_template(csv_path, filename)
%
%--------------------------------------------------------------------------
% Inputs:
% csv_path: (string) Full path to csv containing structure annotation
% information.
%
%--------------------------------------------------------------------------

% Create default filename if not provided
[~, filename] = fileparts(csv_path);
filename = strcat(filename,'_structure_template.csv');

% Template csv must contain these columns
required_vars = {'index', 'acronym', 'name', 'depth', 'parent_structure_id'};
optional_vars = {'id', 'atlas_id', 'structure_id_path',	'color_hex_triplet'};
all_vars = horzcat(required_vars, optional_vars);

% Read table
assert(endsWith(csv_path, '.csv'), "Structure template must be a .csv file")
csv = readtable(csv_path);
col_names = csv.Properties.VariableNames;

% Check for missing columns
for i = 1:length(required_vars)
    if ~any(ismember(col_names, required_vars(i)))
        error("Structure template is missing column ""%s"". Check /annotations/structure_template.csv" +...
            " to see what this column represents", required_vars{i})

    end
end

% Note extra columns that won't be used
for i = 1:length(col_names)
    if ~any(ismember(all_vars, col_names(i)))
        warning("Structure template contains column ""%s"", which is not recognized during template munging", col_names{i})
    end
end

% Being building new structure template
df = table();
df.index = csv.index;

% Get real annotation ids that should match parent structure
if ~any(ismember(col_names, 'id'))
    df.id = csv.index;
    df.atlas_id = csv.index;
else
    df.id = csv.id;
    if any(ismember(col_names, 'atlas_id'))
        df.atlas_id = csv.atlas_id;
    else
        df.atlas_id = csv.id;
    end
end
csv.id = df.id;

% Generate structure_id_path
depth = csv.depth;
if min(csv.depth) < 0 
    depth = depth + abs(min(depth));
end

% Just use depth and parent info to create new tree for each
structure_id_path = repmat({sprintf('/%d',df.id(depth == 0))}, 1, length(df.id));
for i = 1:length(structure_id_path)
    parent_structure = csv.id(i);
    p_depth = depth(i);
    x = repmat("",1,depth(i));
    for j = 1:depth(i)
        if isempty(parent_structure) || ~ismember(parent_structure, csv.id) || p_depth == 0
            continue
        end
        x(p_depth) = parent_structure;
        parent_structure = csv.parent_structure_id(csv.id == parent_structure);
        p_depth = depth(csv.id == parent_structure);
    end
    structure_id_path{i} = sprintf("%s%s/", structure_id_path{i}, sprintf("/%s",x));
end
df.depth = depth;
df.structure_id_path = string(structure_id_path)';

% Generate colors if not defined
if any(ismember(col_names, 'color_hex_triplet'))
    df.color_hex_triplet = csv.color_hex_triplet;
else
    % Default colormap is hsv
    cmap = hsv(length(df.index));
    cmap = string(rgb2hex(cmap(randperm(size(cmap,1)),:)));
    cmap(1) = "#FFFFFF";
    df.color_hex_triplet = cmap;
end

% Add structure names
df.name = csv.name;
df.acronym = csv.acronym;

% Save table to /annotations
home_path = fileparts(which('NM_config'));
output_path = fullfile(home_path, 'annotations', filename);
writetable(df, output_path)
fprintf("Structure template %s saved to /annotations directory!\n", filename)

end