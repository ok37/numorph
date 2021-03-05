function df_temp = calc_stats(df_results, config)
% Run statistics on cell count and volume data for each structure
% Current stats are mean, standard deviation, log2fold change, t test
% p-value, and FDR correced p-value

results_path = config.results_directory;
min_cell_number = config.minimum_cell_number;
compare_structures_by = config.compare_structures_by;
if isfield(config,"structure_csv_path") 
    structure_csv_path = config.structure_csv_path; 
end

markers = config.class_names(config.keep_classes);

% Keep only groups of interest
groups = cat(1,config.groups{:});
par_idx = ismember(groups(:,2),config.comp_groups);
config.samples = config.samples(par_idx);
config.groups = config.groups(par_idx);
group_delimiters = cat(2,config.samples,cat(1,config.groups{:}));

% Create new table with annotations of interest
new_template_path = fullfile(fileparts(which('NM_config')),'annotations','structure_template.csv');
df_temp = readtable(new_template_path);
s_idx = df_temp.index;

% Load any custom structure tables
if isequal(compare_structures_by,'csv')
    df_subset = readtable(fullfile(fileparts(which('NM_config')),'annotations',config.structure_csv));
    s_idx = df_subset.index;
    if ~any(ismember(s_idx,0))
        s_idx = cat(1,0,s_idx);
    end
end

% Check if paried t-test
if isequal(config.paired,"true")
    s1 = unique(group_delimiters(:,4),'rows');
    
    par_idx = false(1,length(s1));
    for i = 1:length(s1)
        idx = group_delimiters(:,4) == s1(i);
        if sum(idx) ~= 2
            par_idx(i) = 1;
        elseif length(unique(group_delimiters(idx,3))) ~= 1
            par_idx(i) = 1;
        end
    end
    group_delimiters = group_delimiters(ismember(group_delimiters(:,4),s1(par_idx)),:);
    group_delimiters = sortrows(group_delimiters,4);
    paired = true;
else
    paired = false;
end

% Get stats
df_stats = pairwise_comparison(df_results,s_idx,markers,group_delimiters,...
    config.comp_groups,min_cell_number,paired,config.custom_comp);

% Combine table
df_stats = horzcat(df_temp,cat(2,df_stats{:}));

% Write new results file
%stats_path = fullfile(results_path,'NMe_summary_stats.csv');
%writetable(df_stats, stats_path)

% Remove structures with no counts
if isequal(config.compare_structures_by,'csv')
   bin_header = strsplit(config.structure_csv,{'/','.csv'});
   bin_header = bin_header{end-1};
   
   % Write new results file
   stats_path = fullfile(results_path,sprintf('NMe_%s_stats.csv',bin_header));
   df_stats_binned = df_stats(ismember(df_stats.index,df_subset.index),:);
   writetable(df_stats_binned, stats_path)
end

end


function df_stats = pairwise_comparison(df_results,s_idx,markers,group_delimiters,group_names,min_cell_number,paired,custom_comp)

% Print
if paired
    fprintf('%s\t Using paired 2 sample t-test for groups %s and %s\n',...
        datetime('now'),group_names(1), group_names(2))    
else    
    fprintf('%s\t Using unpaired 2 sample t-test for groups %s and %s\n',...
        datetime('now'),group_names(1), group_names(2))
end
header_prefix = [repelem(["Mean","StdDev"],2),"PChange","p","p_adj","sig"];

% Split samples into 2 sets
set{1} = group_delimiters(group_delimiters(:,3) == group_names(1),1);
set{2} = group_delimiters(group_delimiters(:,3) == group_names(2),1);

df_stats = cell(1,4);
% Let's look at count stats

if ~isempty(df_results{1})
    df_counts = df_results{1};
    df_counts(~ismember(df_counts.index,s_idx),10:end) = {0};

    % Get column names
    colnames = df_counts.Properties.VariableNames;
    % For each channel, calculate stats
    stats = cell(1,length(markers));
    for i = 1:length(markers)
        group = cell(1,2);
        for j = 1:2
            % For individual cell-types
            idx = contains(colnames,markers(i)) & contains(colnames,set{j});
            group{j} = table2array(df_counts(:,idx));
        end
        
        % Get stats
        sub_stats = get_stats(group{1},group{2},paired,min_cell_number);
        
        % Create headers
        df_header = header_prefix + "_" + repelem(markers(i),8) + "_Counts";
        df_header(1:4) = repmat(group_names,1,2) + "_" + df_header(1:4);
        
        % Convert to table
        stats{i} = array2table(sub_stats,'VariableNames',df_header);
    end
    df_stats{1} = cat(2,stats{:});
end

% Now let's look at volume stats
if ~isempty(df_results{2})
    df_volumes = df_results{2};
    % Get column names
    colnames = df_volumes.Properties.VariableNames;
    % For each group, calculate stats
    group = cell(1,2);
    for j = 1:2
        idx = contains(colnames,[set{j}]);
        group{j} = table2array(df_volumes(:,idx));
    end
    
    % Get stats
    sub_stats = get_stats(group{1},group{2},paired,0);
    
    % Create headers
    df_header = header_prefix + "_" + "Volume";
    df_header(1:4) = repmat(group_names,1,2) + "_" + df_header(1:4);
    % Convert to table
    df_stats{2}  = array2table(sub_stats,'VariableNames',df_header);
end

% Now let's add density if both counts and volumes are present
if ~isempty(df_results{1}) && ~isempty(df_results{2})
    % Get column names
    colnames1 = df_counts.Properties.VariableNames;
    colnames2 = df_volumes.Properties.VariableNames;
    % For each channel, calculate stats
    stats = cell(1,length(markers));
    for i = 1:length(markers)
        group = cell(1,2);
        for j = 1:2
            idx_c = contains(colnames1,[set{j}]) & cellfun(@(s) contains(s,markers(i)),colnames1);            
            cols_c = string(colnames1(idx_c));      
            % Match sample counts and volumes columns
            s_idx = zeros(1,length(cols_c));
            for k = 1:length(cols_c)
                s = strsplit(cols_c(k),'_');
                s_idx(k) = find(startsWith(string(colnames2),s(1)));
            end
            group{j} = table2array(df_counts(:,idx_c));
            group{j} = group{j}./table2array(df_volumes(:,s_idx));
        end
        
        % Get stats
        sub_stats = get_stats(group{1},group{2},paired,0);
        
        % Create headers
        df_header = header_prefix + "_" + repelem(markers(i),8) + "_Density";
        df_header(1:4) = repmat(group_names,1,2) + "_" + df_header(1:4);
        
        % Convert to table
        stats{i} = array2table(sub_stats,'VariableNames',df_header);
    end
    df_stats{3} = cat(2,stats{:});
end

if ~isempty(custom_comp)
    stats = cell(1,length(custom_comp));
    for n = 1:length(custom_comp)
        % Get column names
        colnames = df_counts.Properties.VariableNames;

        % For each channel, calculate stats
        group = cell(length(markers),2);
        custom_func = '@(';
        a = 1; c_idx = false(1,length(markers));
        for i = 1:length(markers)
            if ~contains(custom_comp,strcat('c',string(i)))
                continue
            else
                if a == 1
                    custom_func = strcat(custom_func,string(strcat('c',string(i))));
                else
                    custom_func = strcat(custom_func,string(strcat(',c',string(i))));
                end
                a = 2; c_idx(i)=1;
            end
            for j = 1:2
                % For individual cell-types
                idx = contains(colnames,markers(i)) & contains(colnames,set{j});
                group{i,j} = table2array(df_counts(:,idx));
            end
        end
        
        % Create function
        fh = str2func(strcat(custom_func,')',custom_comp(n)));
        group = group(c_idx,:);
        
        % Apply function
        group{1,1} = fh(group{:,1});
        group{1,2} = fh(group{:,2});
        group = group(1,:);
        
        % Get stats
        sub_stats = get_stats(group{1,1},group{1,2},paired);
        
        % Create headers
        df_header = header_prefix + "_" + repelem(custom_comp,8) + "_Custom";
        df_header(1:4) = repmat(group_names,1,2) + "_" + df_header(1:4);
        
        % Convert to table
        stats{i} = array2table(sub_stats,'VariableNames',df_header);
    end
    df_stats{4} = cat(2,stats{:});
end


% For nuclei, need sum all cell-types
%idx = contains(colnames,[set{j}]) & cellfun(@(s) contains(s,markers),colnames);
%all_counts = table2array(df_counts(:,idx));
%n_samples = sum(idx)/length(markers);
%all_matrix = zeros(size(all_counts,1),n_samples);
%a = 1:length(markers);
%for k = 1:n_samples
%    all_matrix(:,k) = sum(all_counts(:,a),2);
%    a = a+length(markers);
%end
%group{j} = all_matrix;                

end


function stats = get_stats(group1,group2,paired,min_stat)

if nargin<4
    min_stat=0;
end
stats = zeros(size(group1,1),8);

% Mean
stats(:,1) = mean(group1,2);
stats(:,2) = mean(group2,2);
% Standard deviation
stats(:,3) = std(group1,0,2);
stats(:,4) = std(group2,0,2);
% Percent Relative Change
stats(:,5) = 100*(stats(:,2)-stats(:,1))./stats(:,1);    
% Set any cells in background to 0
stats(1,:) = 0;
% p-value
if paired
    [~,stats(:,6)] = ttest(group1',group2'); 
else            
    [~,stats(:,6)] = ttest2(group1',group2','Vartype','unequal'); 
end
% q-value
% First remove structures with low counts/no volume
if nargin<4
    s_pos = ~isnan(stats(:,6));
else    
    s_pos = max(stats(:,1:2),[],2) > min_stat;
end
p_val_thresholded = stats(s_pos,6);
[stats(s_pos,8), ~, ~, stats(s_pos,7)]=fdr_bh(p_val_thresholded);

stats(s_pos,8) = stats(s_pos,8) + single(stats(s_pos,7)<0.05);
stats(s_pos,8) = stats(s_pos,8) + single(stats(s_pos,7)<0.01);
stats(s_pos,8) = stats(s_pos,8) + single(stats(s_pos,7)<0.001);
stats(s_pos,8) = stats(s_pos,8) + single(stats(s_pos,7)<0.0001);

stats(~s_pos,1:5) = 0;
stats(~s_pos,6:7) = 1;
stats(~s_pos,8) = 0;

end
