function df_stats = calc_stats(df_results, config)
% Run statistics on cell count and volume data for each structure
% Current stats are mean, standard deviation, log2fold change, t test
% p-value, and FDR correced p-value

results_path = config.results_directory;
group_delimiters = config.group_delimiters;
markers = config.markers;
minimum_cell_number = config.minimum_cell_number;
overwrite = config.overwrite;
compare_structures_by = config.compare_structures_by;
if isfield(config,"structure_csv_path") structure_csv_path = config.structure_csv_path; end

% Create new table with containing the groups of interest
new_template_path = './supplementary_data/structure_template.csv';
df_stats = readtable(new_template_path);
s_idx = df_stats.index;

% Load any custom structure tables
if isequal(compare_structures_by,'csv')
    df_subset = readtable(structure_csv_path);
    s_idx = df_subset.index;
    if ~any(ismember(s_idx,0))
        s_idx = cat(1,0,s_idx);
    end
     %s_idx = ismember(df_stats.index,df_subset.index);
     %df_results = cellfun(@(s) s(s_idx,:),df_results,'UniformOutput',false);
end

% Let's look at count stats
if ~isempty(df_results{1})
    df_counts = df_results{1};
    df_counts(~ismember(df_counts.index,s_idx),10:end) = {0};
    
    % Get column names
    colnames = df_counts.Properties.VariableNames;
    % For each channel, calculate stats
    stats = cell(1,length(markers));
    for i = 1:length(markers)
        stats{i} = zeros(size(df_counts,1), 8);
        group = cell(1,length(group_delimiters));
        for j = 1:length(group_delimiters)
            if i == 1
                idx = contains(colnames,group_delimiters(j)) & cellfun(@(s) contains(s,markers),colnames);
                all_counts = table2array(df_counts(:,idx));
                n_samples = sum(idx)/length(markers);
                all_matrix = zeros(size(all_counts,1),n_samples);
                a = 1:length(markers);
                for k = 1:n_samples
                    all_matrix(:,k) = sum(all_counts(:,a),2);
                    a = a+length(markers);
                end
                group{j} = all_matrix;
            else
                idx = contains(colnames,markers(i)) & contains(colnames,group_delimiters(j));
                group{j} = table2array(df_counts(:,idx));
            end
        end
        % Mean
        stats{i}(:,1) = mean(group{1},2);
        stats{i}(:,2) = mean(group{2},2);
        % Standard deviation
        stats{i}(:,3) = std(group{1},0,2);
        stats{i}(:,4) = std(group{2},0,2);
        % Percent Change
        stats{i}(:,5) = 100*(stats{i}(:,2)-stats{i}(:,1))./stats{i}(:,1);

        %stats{i}(:,5) = log2(abs(stats{i}(:,2)- stats{i}(:,1))+1).*sign(stats{i}(:,2)- stats{i}(:,1));
        
        % Set any cells in background to 0
        stats{i}(1,:) = 0;
        % p-value
        [~,stats{i}(:,6)] = ttest2(group{1}',group{2}','Vartype','unequal'); 
    end
    % FDR adjustment 
    % First remove structures with less than minimum number of cells
    s_pos = max(stats{1}(:,1:2) > minimum_cell_number,[],2);
    for i = 1:length(stats)
        p_val_thresholded = stats{i}(s_pos,6);
        [~, ~, ~, stats{i}(s_pos,7)]=fdr_bh(p_val_thresholded);
        %[stats{i}(s_pos,7)]=mafdr(p_val_thresholded,'BHFDR','true');
        
        stats{i}(s_pos,8) = stats{i}(s_pos,8) + single(stats{i}(s_pos,7)<0.05);
        stats{i}(s_pos,8) = stats{i}(s_pos,8) + single(stats{i}(s_pos,7)<0.01);
        stats{i}(s_pos,8) = stats{i}(s_pos,8) + single(stats{i}(s_pos,7)<0.001);
        stats{i}(s_pos,8) = stats{i}(s_pos,8) + single(stats{i}(s_pos,7)<0.0001);
        
        stats{i}(~s_pos,1:5) = 0;
        stats{i}(~s_pos,6:7) = 1;
        stats{i}(~s_pos,8) = 0;
    end
    % Now append results to stats table
    header_prefix = [repelem(["Mean","StdDev"],length(group_delimiters)),"PChange","p","p_adj","sig"];
    for i = 1:length(markers)
        % Create headers
        df_header = header_prefix + "_" + repelem(markers(i),8) + "_Counts";
        df_header(1:4) = repmat(group_delimiters,1,2) + "_" + df_header(1:4);
        % Convert to table
        stats_table = array2table(stats{i},'VariableNames',df_header);
        % Concatenate stats to table
        df_stats = horzcat(df_stats, stats_table);
    end
end

% Now let's look at volume stats
if ~isempty(df_results{2})
    df_volumes = df_results{2};
    %df_volumes(~ismember(df_volumes.index,s_idx),10:end) = {0};

    % Get column names
    colnames = df_volumes.Properties.VariableNames;
    % For each group, calculate stats
    stats = zeros(size(df_volumes,1), 8);
    group = cell(1,length(group_delimiters));
    for j = 1:length(group_delimiters)
        idx = contains(colnames,group_delimiters(j));
        group{j} = table2array(df_volumes(:,idx));
    end
    % Mean
    stats(:,1) = mean(group{1},2);
    stats(:,2) = mean(group{2},2);
    % Standard deviation
    stats(:,3) = std(group{1},0,2);
    stats(:,4) = std(group{2},0,2);
    % Percent Relative Change
    stats(:,5) = 100*(stats(:,2)-stats(:,1))./stats(:,1);
    %stats(:,5) = log2(abs(stats(:,2)- stats(:,1))+1).*sign(stats(:,2)- stats(:,1));
    
    % Set any cells in background to 0
    stats(1,:) = 0;
    % p-value
    [~,stats(:,6)] = ttest2(group{1}',group{2}','Vartype','unequal');         
    % q-value
    % First remove structures with no volume
    s_pos = max(stats(:,1:2) > 0,[],2);
    p_val_thresholded = stats(s_pos,6);
    [stats(s_pos,8), ~, ~, stats(s_pos,7)]=fdr_bh(p_val_thresholded);
    
    stats(s_pos,8) = stats(s_pos,8) + single(stats(s_pos,7)<0.05);
    stats(s_pos,8) = stats(s_pos,8) + single(stats(s_pos,7)<0.01);
    stats(s_pos,8) = stats(s_pos,8) + single(stats(s_pos,7)<0.001);
    stats(s_pos,8) = stats(s_pos,8) + single(stats(s_pos,7)<0.0001);
    
    stats(~s_pos,1:5) = 0;
    stats(~s_pos,6:7) = 1;
    stats(~s_pos,8) = 0;
    
    % Now append results to stats table and save
    header_prefix = [repelem(["Mean","StdDev"],length(group_delimiters)),"PChange","p","p_adj","sig"];
    % Create headers
    df_header = header_prefix + "_" + "Volume";
    df_header(1:4) = repmat(group_delimiters,1,2) + "_" + df_header(1:4);
    % Convert to table
    stats_table = array2table(stats,'VariableNames',df_header);
    % Concatenate stats to table
    df_stats = horzcat(df_stats, stats_table);
end

% Now let's add density if both counts and volumes are present
if ~isempty(df_results{1}) && ~isempty(df_results{2})
    % Get column names
    colnames1 = df_counts.Properties.VariableNames;
    colnames2 = df_volumes.Properties.VariableNames;
    % For each channel, calculate stats
    stats = cell(1,length(markers));
    for i = 1:length(markers)
        stats{i} = zeros(size(df_counts,1), 8);
        group = cell(1,length(group_delimiters));
        for j = 1:length(group_delimiters)
            idx_c = contains(colnames1,markers(i)) & contains(colnames1,group_delimiters(j));            
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
        % Mean
        stats{i}(:,1) = mean(group{1},2);
        stats{i}(:,2) = mean(group{2},2);
        % Standard deviation
        stats{i}(:,3) = std(group{1},0,2);
        stats{i}(:,4) = std(group{2},0,2);
        % Percent Relative Change
        stats{i}(:,5) = 100*(stats{i}(:,2)-stats{i}(:,1))./stats{i}(:,1);
        %stats{i}(:,5) = log2(abs(stats{i}(:,2)- stats{i}(:,1))+1).*sign(stats{i}(:,2)- stats{i}(:,1));
        % Set any cells in background to 0
        stats{i}(1,:) = 0;
        % p-value
        %for j = 1:size(group{1},1)
        %    [~,stats{i}(j,6)] = ttest2(group{1}(j,:),group{2}(j,:),'Vartype','unequal'); 
        %end
        [~,stats{i}(:,6)] = ttest2(group{1}',group{2}','Vartype','unequal'); 
    end
    % q-value
    % First remove structures with no volume
    s_pos = max(table2array(df_stats(:,10:11)) > 0,[],2);
    for i = 1:length(stats)
        p_val_thresholded = stats{i}(s_pos,6);
        [~, ~, ~, stats{i}(s_pos,7)]=fdr_bh(p_val_thresholded);
        
        stats{i}(s_pos,8) = stats{i}(s_pos,8) + single(stats{i}(s_pos,7)<0.05);
        stats{i}(s_pos,8) = stats{i}(s_pos,8) + single(stats{i}(s_pos,7)<0.01);
        stats{i}(s_pos,8) = stats{i}(s_pos,8) + single(stats{i}(s_pos,7)<0.001);
        stats{i}(s_pos,8) = stats{i}(s_pos,8) + single(stats{i}(s_pos,7)<0.0001);
        
        stats{i}(~s_pos,1:5) = 0;
        stats{i}(~s_pos,6:7) = 1;
        stats{i}(~s_pos,8) = 0;
    end
    % Now append results to stats table and save
    header_prefix = [repelem(["Mean","StdDev"],length(group_delimiters)),"PChange","p","p_adj","sig"];
    for i = 1:length(markers)
        % Create headers
        df_header = header_prefix + "_" + repelem(markers(i),8) + "_Density";
        df_header(1:4) = repmat(group_delimiters,1,2) + "_" + df_header(1:4);
        % Convert to table
        stats_table = array2table(stats{i},'VariableNames',df_header);
        % Concatenate stats to table
        df_stats = horzcat(df_stats, stats_table);
    end
end

% Write new results file
stats_path = fullfile(results_path,'TCe_summary_stats.csv');
writetable(df_stats, stats_path)


% Remove structures with no counts
if isequal(config.compare_structures_by,'csv')
   bin_header = strsplit(config.structure_csv_path,{'/','.csv'});
   bin_header = bin_header{end-1};
   
   % Write new results file
   stats_path = fullfile(results_path,sprintf('TCe_%s_stats.csv',bin_header));
   df_stats_binned = df_stats(ismember(df_stats.index,df_subset.index),:);
   writetable(df_stats_binned, stats_path)
end




end