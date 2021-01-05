function ct = classify_cells_threshold(centroids, config)
%--------------------------------------------------------------------------
% Classify cell-types using intensity thresholding
%--------------------------------------------------------------------------
remove_l1 = 'true';             % Remove cortical layer 1 cells
stratify_structures = 'true';    % Perform GMM on each structure seperately

% Get parameters from config structure
markers = config.markers;

n_markers = length(markers);

% Create empty vector for storing cell-types
ct = ones(size(centroids,1),1);

% Remove any centroids without annotations
rm_idx = centroids(:,4) == 0;
fprintf('%s\t Removed %d nuclei with no annotation \n',datetime('now'),sum(rm_idx))

% Remove cells below minimum nuclei intensity threshold
low_idx = centroids(:,5) < config.min_intensity;
rm_idx = rm_idx | low_idx;
fprintf('%s\t Removed %d low intensity nuclei \n',datetime('now'),sum(low_idx))

% Ask if to remove layer 1 prior to gmm analysis
if isequal(remove_l1,"true")
   l1_table = readtable(fullfile(config.home_path,'annotations','cortical_layers',...
       'layer1.csv'));
   l1_idx = ismember(centroids(:,4),l1_table.index);
  
   s1 = centroids(:,6) > prctile(centroids(l1_idx,6),75);
   s2 = centroids(:,7) > prctile(centroids(l1_idx,7),75);

   l1_idx = s1 & s2 & l1_idx;
   rm_idx = rm_idx | l1_idx;
   fprintf('%s\t Removed %d cells \n',datetime('now'),sum(l1_idx))
else
    l1_idx = 0;
end

% Remove selected indexes
centroids = centroids(~rm_idx,:);

% Intitialize matrices and models
count = zeros(size(centroids,1),n_markers+1);
for i = 2:n_markers
    % Clustering using Gaussian Mixtures
    idx = 4+i;
    intensities = centroids(:,idx);

    % Load background intensities from background images
    if nargin > 1
       fprintf('%s\t Subtracting background from background image \n',datetime('now'))
        % If given background images, get background intensities from these
        %[back_thresh, ~] = get_background_at_centroids(centroids,...
        %    config,B{i},S{i});
        back_vals = readmatrix(fullfile(config.output_directory,'background','back_values.csv'));
        back_vals = back_vals(~rm_idx,:);
        back_thresh = back_vals(:,i-1);
        
    else
        % Otherwise use constant background
        back_thresh = 0;
    end
    
    % Subtract background
    intensities = intensities - back_thresh;
    intensities(intensities<0) = 0;
    values = log2(intensities +1);
    
    % Apply z normalization
    if isequal(config.z_normalization,"true")
       fprintf('%s\t Applying z normalization \n',datetime('now'))
        values = (values - mean(values))/std(values);
        if ~isempty(config.log_outliers) || config.log_outliers ~= 0
            fprintf('%s\t Supressing outliers \n',datetime('now'))
            l = config.log_outliers;
            values(values>l) = l + log2(values(values>l));
            values(values<-l) = -l - log2(abs(values(values<-l)));
        end
    end
    
    % Calculate threshold
    if ~isempty(config.expression)
        % From expression
        if iscell(config.expression)
            expression = config.expression{i-1};
        else
            expression = config.expression;
        end
        thresh = calculate_threshold_expression(values, expression);
        fprintf('%s\t Using expression threshold %f on marker %s \n',...
            datetime('now'),thresh, markers(i))
    else
        % From manually specified value
        assert(length(config.thresholds) == length(config.markers)-1,...
            "Found %d thresholds for %d markers to be classified",...
            length(config.thresholds),length(config.markers)-1)
        assert(length(config.thresholds(i-1)) < 65535,...
            "Threshold set above 16-bit limit")
        if config.thresholds(i-1) > 0 && config.thresholds(i-1) < 1
            thresh = config.thresholds(i-1)*65535;
        elseif config.threshold(i-1) > 1 
            thresh = config.thresholds(i-1);
        else
            thresh = 0;
        end
        fprintf('%s\t Using manual threshold %f on marker %s \n',...
            datetime('now'),thresh, markers(i))
    end
    
    % Count cells above threshold
    count(:,i) = values>thresh;
end

if n_markers > 2
    % Count cells positive for all markers and add to last column
    count(:,end) = all(count(:,2:end-1),2);
    count(:,2:end-1) = count(:,2:end-1) - count(:,end);
end

% Calculate all negative cells
count(:,1) = ~any(count(:,2:end),2);

% Save cell type counts
[r,c] = find(count);
v1(r) = c;
ct(~rm_idx) = v1;
ct(rm_idx) = 0;

% Calculate sums
max_ct = max(unique(ct));
sums = histcounts(ct,-0.5:max_ct+0.5);
sums = sums(2:end);
sums(1) = sums(1) + sum(l1_idx);
pct = sums/sum(sums(:));
disp(pct)

end


function threshold = calculate_threshold_expression(values, expression)
%--------------------------------------------------------------------------
% Calculate threshold from defined value or regular expression
%--------------------------------------------------------------------------

if isempty(expression)
    return
end

% Deconstruct regular expression
splitstring = regexp(expression,'\s','split');

threshold = 0;
for i = 1:length(splitstring)
    s1 = regexp(splitstring(i),'*','split');
    
    % Constant to multiply by
    if length(s1) < 2
        c = 1;
    else
        c = str2double(s1(1));
        s1(1) = [];
    end
    
    % Some stat from data
    s2 = regexp(s1,'(\w+)','match');
    if isempty(s2)
        continue
    end
    switch s2(1)
        case 'mean'
            stat = mean(values);
        case 'median'
            stat = median(values);
        case 'mode'
            stat = mode(values);
        case 'std'
            stat = std(values);
        case 'mad'
            stat = mad(values,1);
        case 'prctile'
            stat = prctile(values,str2double(s2(2)));
        case '+'
            continue
        otherwise
            warning("Could not recognize regular expression term %d\n",i)
            continue
    end
    % Add to threshold value
    threshold = threshold + stat*c;
end

end