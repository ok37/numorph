%% Specify Directory Information
% Path or paths to folder(s) containing individual sample cell counting 
% results as .csv files (i.e. SAMPLE_volumes.csv, SAMPLE_classes.csv, SAMPLE_centroids.csv)
sample_directory = [];

% Results directory 
% Statistics saved as: NMe_summary_stats.csv
results_directory = [];

% Marker name for each unique cell-type annotation
use_classes = "false";       % Use cell-type classes. Otherwise only raw number of cell will be evaluated
keep_classes = 1:3;          % If using classes, select which ones to evaluate
class_names = ["topro","ctip2","cux1"]; % Assign names to classes

% Labels for distinguishing WT/KO groups in .csv files
combine_counts = 'true';               % Combine counts from multiple samples
combine_volumes = 'true';               % Combine volumes from multiple samples

calculate_stats = 'true';              % Create summary table with statistics
comp_groups = ["WW","FF"];             % Select groups to compare. Leave empty to compare all groups
custom_comp = "c1./(c1+c2+c3)";        % Do custom comparison
paired = "false";                       % Perform paired t-test (for example to control for litter). Set pairing order as 3rd entry in 'group' in NM_samples
use_data = "both";                     % Use counts or volumes or both for quantification

compare_structures_by = 'csv';         % Select which structures for multiple comparisons. Options are: tree, csv. csv will take a custom list of structure indexes
structure_csv = 'harris_cortical_groupings_w_iso.csv';

visualize_results = 'true';             % Create volume visualization from stats
orientation = 'coronal';                

minimum_cell_number = 100;              % Minimum number of cells for each structure to be compared. Recommended minimum: 1 
structure_depth = [];                   % Minimum structure graph order for comparisons (1-10). Recommended minimum: 5
