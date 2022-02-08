%% Evaluation template
% Outputs are saved to a 'results' folder located in numorph's working 
% directory, unless specified otherwise. Note: this template only a single
% 'update' parameter to overwrite all statistics calculated.

% Specify which structures to compare
update = "true";                                         % Whether to overwrite previous measurements
results_directory = [];                                  % Output directory for summary (default: numorph/results)
compare_structures_by = "table";                         % index, table; Compare all unique annotations (index), structures according to table (table)
structure_table = "harris_cortical_groupings_w_iso.xls"; % table in /annotations/custom_annotations indicating structures to evaluate
measure_cortex = "true";                                 % Optional: Take cortex measurments. Will only work with default ccfv3 annotations

% For cell counting, specify cell-type class information
use_classes = "true";                           % Use cell-type classes. Otherwise raw centroid counts will be used
keep_classes = 1:3;                             % If using classes, select which ones to evaluate by specifying the class index
class_names = [];                               % Assign names to classes
sum_all_classes = "true";                       % "true","false"; Create class that sums all classes detected (i.e. all nuclei)
custom_class = ["CT1./(CT1+CT2+CT3)",...        % Create a function to combine multiple classes (i.e. calculate cell-type proportions). Specify as class names
    "CT2./(CT1+CT2+CT3)",...
    "CT3./(CT1+CT2+CT3)",...
    "CT2+CT3"];                              
minimum_cell_number = 10;                        % Minimum number of cells for each structure to be compared

% Specify groups for comparisons
compare_groups = ["WW","FF"];                   % Select groups to compare. Defined by the 2nd value in sample 'group' variable
paired = "false";                               % Performed paired t-test. Samples paired by the 3rd value in sample 'group' variable
