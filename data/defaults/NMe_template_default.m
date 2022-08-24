%% Evaluation template
% Outputs are saved to a 'results' folder located in numorph's working
% directory, unless specified otherwise. Note: this template only a single
% 'update' parameter to overwrite all statistics calculated.
update = "true";                    % Whether to overwrite previous measurements
results_directory = [];             % Output directory for summary (default: numorph/results)
calculate_voxel_stats = "true";     % Calculate voxel-level statistics                        
compare_groups = ["WT","KO"];      % Select groups to compare. Defined by the 2nd value in sample 'group' variable

% Specify which structures to compare
annotation_file = "ccfv3.mat";        % name of annotation data file from NMa_template
template_file = "structure_template.csv";  % table in /annotations containing all possible structure indexes present (if comparing by table)
compare_structures_by = "index";        % index, table; Compare all unique annotations (index), structures according to table (table)
structure_table = [];                   % For comparing by "table": .xls table in /annotations/custom_annotations indicating structures to evaluate. All indexes must exist in structure_template

% For cell counting, specify cell-type class information
use_classes = "true";               % Use cell-type classes. Otherwise raw centroid counts will be used
keep_classes = 1:3;                 % If using classes, select which ones to evaluate by specifying the class index
class_names = [];                   % Assign names to each class. Otherwise, leave empty and class indexes are prefixed with "CT"
minimum_cell_number = 10;           % Minimum number of cells for each structure to be compared
custom_class = [];                  % Specify a function to combine multiple classes (i.e. sum 2 classes or calculate cell-type proportions)

% Custom class examples:
% This calculates the proportion of cell-type 1 to 3 other cell-types:
% custom_class = ["CT1/(CT1+CT2+CT3"]; 
%
% This calculates the sum of cell-types 2+3
% custom_class = ["CT2+CT3"];
%
% To do both:
% custom_class = ["CT1/(CT1+CT2+CT3", "CT2+CT3"];