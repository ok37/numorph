function [ct, p, gm] = count_celltypes(centroids, config, B, S)
% Count cells using local or global threshold
config.markers = ["ToPro" "Ctip2" "Cux1"];
remove_l1 = "false";
combined = "true";
use_structures = "true";
surpress_outliers = "false";

n_clusters = config.n_clusters;
confidence = config.confidence;

n_markers = length(config.markers);

% Create empty vector for storing cell-types
ct = ones(size(centroids,1),1);

% Remove any centroids without annotations
rm_idx = centroids(:,4) == 0;

% Ask if to remove layer 1 prior to gmm analysis
if isequal(remove_l1,"true")
   l1_table = readtable(fullfile(config.home_path,'annotations','cortical_layers',...
       'layer1.csv'));
   l1_idx = ismember(centroids(:,4),l1_table.index);
   
   s1 = centroids(:,5)<200;
   s2 = centroids(:,6)>prctile(centroids(:,6),95);
   s3 = centroids(:,7)>prctile(centroids(:,7),95);

   l1_idx =  s1;%| s2 & s3;   
   rm_idx = rm_idx | l1_idx;
   
  fprintf('%s\t Removed %d cells \n',datetime('now'),sum(rm_idx))
else
    l1_idx = 0;
end

% Remove selected indexes
centroids = centroids(~rm_idx,:);

% Intitialize matrices and models
gm = cell(1,n_markers-1);
p = cell(1,n_markers-1);
count = zeros(size(centroids,1),n_markers+1);
for i = 2:n_markers
    % Clustering using Gaussian Mixtures
    idx = 4+i;
    intensities = centroids(:,idx);

    % If thresholds aren't present in centroid list, load from background
    % images
    if i ~=1%idx+n_markers < size(centroids,2)
        % Check centroids table for extra columns
        back_thresh = centroids(:,end-(3-i));
        %std_thresh = centroids(:,idx+n_markers);
        %disp(mean(back_thresh(:)))
        b = 0;
    elseif nargin > 2 && ~isempty(B{i})
        % If given background images, get background intensities from these
        [back_thresh, std_thresh] = get_background_at_centroids(centroids,...
            config,B{i},S{i});
        b = 1;
    else
        % Otherwise use constant background
        back_thresh = repmat(100,size(centroids,1),1);
        std_thresh = repmat(5,size(centroids,1),1);
        b = 0;
    end
    
    % z-normalization
    intensities = intensities - back_thresh*b;	
    values = intensities;
%    values = (values - mean(values))/std(values);
%    values(values>10) = 10 + log10(values(values>10));
%    values(values<-10) = -10 - log10(abs(values(values<-10)));
    
    if isequal(combined,"false")
        fprintf('%s\t Running GMM on marker %s with %d clusters \n',datetime('now'),config.markers(i),n_clusters(i-1))
        gm{i-1} = fitgmdist(values,n_clusters(i-1),'Options',statset('MaxIter',400));
        
        % Find posteriors
        p{i-1} = posterior(gm{i-1},values);
        
        % Find background cluster
        [~,back_idx] = sort(gm{i-1}.mu);
        
        % Threshold by confidence of NOT being in background cluster
        back_calls = p{i-1}(:,back_idx(1))>(1-confidence(i-1));
        count(:,i) = ~back_calls;
    else
        val_save(:,i-1) = values;
        if i ~= n_markers
            continue
        end
        
        if ~isempty(config.mix_proportions)
            init_struct = generate_intial_conditions(val_save,config.mix_proportions,config.mix_markers);
        else
            init_struct = 'plus';
        end

        if ~isequal(use_structures,"true")
            fprintf('%s\t Running GMM all markers with %d clusters \n',datetime('now'),n_clusters(1))
            gm = fitgmdist(val_save,n_clusters(1),'Start',init_struct,'Options',statset('MaxIter',200));
            % Posteriors
            p = posterior(gm,val_save);
            %clust = cluster(gm,val_save); 
        else
            structures = bin_annotation_structures(centroids(:,4));
            [p,gm] = cluster_by_structure(val_save,structures,n_clusters(1),config);
        end

        % Finalize counts
        count = finalize_counts(p,config.confidence);
    end
        
end

if n_markers > 2 && isequal(combined,"false")
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


function [back_thresh, std_thresh] = get_background_at_centroids(centroids,config,B,S)

% Convert centroid positions to background image space
cen_location = centroids(:,1:3).*config.resolution;
cen_mask = round(cen_location./config.back_res);
res_dims = size(B);

% Limit any indexes that fall out of range of the background image
% after round
cen_mask(any(cen_mask==0)) = 1;
cen_mask(cen_mask(:,1) > res_dims(1)) = res_dims(1);
cen_mask(cen_mask(:,2) > res_dims(2)) = res_dims(2);
cen_mask(cen_mask(:,3) > res_dims(3)) = res_dims(3);

% Get linear indexes for centroid positions in the background image
cen_idx = sub2ind(res_dims, cen_mask(:,1),cen_mask(:,2), cen_mask(:,3));
        
% Append values for background and standard deviation
back_thresh = single(B(cen_idx));
std_thresh = single(S(cen_idx));

end

function init_struct = generate_intial_conditions(intensities,mix_proportions,mix_markers)
% Estimate GMM intial conditions using estimated mixing proportions

n_cells = size(intensities,1);

% Split cells into subsets defined by mixing proportions
subset = cell(1,length(mix_markers)); 
for i = 1:length(mix_markers)
    c_idx = mix_markers{i};
    if isempty(c_idx)
        % Background
        sorted = sort(intensities,'ascend');
        subset{i} = sorted(1:round(n_cells*mix_proportions(i)),:);
    elseif length(c_idx) > 1
        % Double channel
        sub1 = prod(intensities(:,c_idx-1),2);
        subset{i} = intensities(sub1 > prctile(sub1,(1-mix_proportions(i))*100),:);
    else
        % Single channel
        subset_idx = intensities(:,c_idx-1) > prctile(intensities(:,c_idx-1),(1-mix_proportions(i))*100);
        subset{i} = intensities(subset_idx,:);
    end    
end

% Determine parameters
for i = 1:length(subset)
    init_struct.ComponentProportion(i) = size(subset{i},1)/n_cells;
    init_struct.mu(i,:) = median(subset{i});
    init_struct.Sigma(:,:,i) = cov(subset{i});    
end

end

function [p,gm] = cluster_by_structure(val_save, structures, n_clusters,config)

p = zeros(size(val_save,1),n_clusters);
u_struct = unique(structures);

for i = 1:length(u_struct)
    if u_struct(i) == 0
        continue
    end
    fprintf('%s\t Running GMM on all markers with %d clusters on structure %d\n',...
        datetime('now'),n_clusters(1),i)
    
    idx = structures == u_struct(i);
    vals_sub = val_save(idx,:);

    %vals_sub = (vals_sub - mean(vals_sub,1))./std(vals_sub,[],1);
    %vals_sub(vals_sub>10) = 10 + log10(vals_sub(vals_sub>10));
    %vals_sub(vals_sub<-10) = -10 - log10(abs(vals_sub(vals_sub<-10)));


    init_struct = generate_intial_conditions(vals_sub,config.mix_proportions,config.mix_markers);
    
    try
        gm = fitgmdist(vals_sub,n_clusters(1),'Start',init_struct,'Options',statset('MaxIter',400));
    catch
	disp("ehllo")
        gm = fitgmdist(vals_sub,n_clusters(1),'Options',statset('MaxIter',400));
    end

    p_sub = posterior(gm,vals_sub);
    
    p(idx,:) = p_sub;
end

end
