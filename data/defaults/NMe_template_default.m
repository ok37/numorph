%% Specify Directory Information
groups = ["WT","TOP1"];                 % Labels for WT/KO groups in .csv files. 2 values required

% Path or paths to folder(s) containing individual sample cell counting 
% results as .csv files. 
% Cell counts csv files should be named: WT1_centroids.csv
% Structure volume .csv files should be named: WT1_volumes.csv
sample_directory = "/evaluation/Results_svm";

% Results directory. Important: should be different than sample_directory
% Cell count results saved as or appended to: TCe_summary_counts.csv
% Volume measurement saved as or appended to: TCe_summary_volumes.csv
% Statistics saved as: TCe_summary_stats.csv
results_directory = "./evaluation/";

% Sample names. Corresponding files should be named like
% sample_centroids.csv
sample_names = ["WT11L","WT1L","WT7R","WT8R",...
    "TOP110R","TOP11L","TOP14R","TOP16R"];

% Marker name for each unique cell-type annotation
% By default, marker 1 is given as all cells
markers = ["ToPro","Ctip2","Cux1"];
ct_column = 9;          % Column containing cell-type annotations. By default, last column is chosen

% Labels for distinguishing WT/KO groups in .csv files

combine_counts = 'load';               % Combine counts from multiple samples
combine_volumes = 'load';               % Combine volumes from multiple samples
sum_child_structures = 'true';          % Whether to sum counts for parent structures according to structure tree. Recommended true

calculate_stats = 'true';              % Create summary table with statistics
use_data = 'both';                      % Use counts or volumes or both for quantification
compare_structures_by = 'csv';         % Select which structures for multiple comparisons. Options are: tree, csv. csv will take a custom list of structure indexes

% Optional. Location of csv file containing custom list of structures to compare
%structure_csv_path = "/annotations/cortex_ontology_groupings.csv";
structure_csv_path = '/annotations/cortex_17regions.csv';

overwrite = 'true';                     % Overwrite existing table. Otherwise results will append to exisiting file
visualize_results = 'true';             % Create volume visualization from stats
orientation = 'coronal';                %

minimum_cell_number = 100;              % Minimum number of cells for each structure to be compared. Recommended minimum: 1 
structure_depth = [];                   % Minimum structure graph order for comparisons (1-10). Recommended minimum: 5
