function calc_voxel_stats(config)
%--------------------------------------------------------------------------
% Create 4D voxel volume representing fold change, p, q values
%--------------------------------------------------------------------------

% Load voxel volumes
n_samples = length(config.results_path);
imgs = cell(1,n_samples);
for i = 1:length(imgs)
    var_names = who('-file',config.results_path(i));
    if ~ismember('voxel_volume',var_names)
        results = load(config.results_path(i));
        results = create_voxel_volume(results);
        
        
        %error("No voxel volume found for sample %s",config.samples(i))
    end
    load(config.results_path(i),'voxel_volume')
    voxel_volume = voxel_volume(config.keep_classes);
    if isequal(config.contains_nuclear_channel,"true")
        imgs{i} = [voxel_volume,{sum(cat(4,voxel_volume{:}),4)}];
    else
        imgs{i} = voxel_volume;
    end
end

% Add custom classes
if ~isempty(config.custom_class)
    n_single = length(config.class_names(config.keep_classes));
    n_custom = length(config.custom_class);
    for i = 1:length(imgs)
        counts = cellfun(@(s) s(:),imgs{i},'UniformOutput',false);
        counts = cat(2,counts{:});
        counts = get_custom_class(counts,1:n_single,...
            config.custom_class);
        for j = 1:n_custom
            a = {reshape(counts(:,n_single+j),size(imgs{i}{1}))};
            imgs{i} = [imgs{i},a];
        end
    end
end

% Check all images are the same size
n_vox = cellfun(@(s) numel(s),[imgs{:}]);
assert(length(unique(n_vox)) == 1,"Voxel volume sizes are not the same for all volumes")

% Assign samples to groups
group_delimiters = cat(2,config.samples,cat(1,config.groups{:}));
set1 = group_delimiters(group_delimiters(:,3) == config.compare_groups(1),1);
vox1 = imgs(ismember(config.samples,set1));
set2 = group_delimiters(group_delimiters(:,3) == config.compare_groups(2),1);
vox2 = imgs(ismember(config.samples,set2));

img_sizes = cellfun(@(s) size(s),vox1{1},'UniformOutput',false);

vox1 = arrayfun(@(s) s{:}(:),cat(1,vox1{:}),'UniformOutput',false);
vox2 = arrayfun(@(s) s{:}(:),cat(1,vox2{:}),'UniformOutput',false);

% For each class, generate vox
n_classes = length(imgs{1});

% Get samples
for i = 1:n_classes
    vox_res = ones(prod(img_sizes{i}),3);
    fvox_res = ones(1360*1360,3);
    
    % Seperate groups
    set1 = cat(2,vox1{:,i});
    set2 = cat(2,vox2{:,i});
    
    % Generate flatmaps
    fset1 = cell(1,size(set1,2));
    for j = 1:size(set1,2)
        s = convert_to_flatmap(reshape(set1(:,j),img_sizes{i}));
        fset1{j} = s(:);
    end
    fset1 = cat(2,fset1{:});
    
    fset2 = cell(1,size(set2,2));
    for j = 1:size(set2,2)
        s = convert_to_flatmap(reshape(set2(:,j),img_sizes{i}));
        fset2{j} = s(:);
    end
    fset2 = cat(2,fset2{:});
    
    % Calculate means
    m_set1 = mean(set1,2);
    m_set2 = mean(set2,2);
    fm_set1 = mean(fset1,2);
    fm_set2 = mean(fset2,2);
    
    % Remove voxels with low cells
    idx = find(max(m_set1,[],2) & max(m_set2,[],2) >config.minimum_cell_number);
    fidx = find(max(fm_set1,[],2) & max(fm_set2,[],2) >config.minimum_cell_number);
    
    % Fold change
    vox_res(idx,1) = m_set2(idx)./m_set1(idx);
    fvox_res(fidx,1) = fm_set2(fidx)./fm_set1(fidx);
    
    % p value
    [~,p] = ttest2(set1(idx,:)',set2(idx,:)');
    vox_res(idx,2) = p;

    [~,pf] = ttest2(fset1(fidx,:)',fset2(fidx,:)');
    fvox_res(fidx,2) = pf;

    %p = zeros(1,length(idx));
    %for j = 1:length(idx)
    %    [~,p(j)] = ttest2(set1(idx(j),:),set2(idx(j),:)); 
    %end
    %vox_res(idx,2) = p;
    
    %for j = 1:length(fidx)
    %    [~,pf(j)] = ttest2(fset1(fidx(j),:),fset2(fidx(j),:)); 
    %end
    %fvox_res(fidx,2) = pf;

    % Adjusted p value. Apply threshold to adjusted p value
    [~,~,~, adj_p] = fdr_bh(p);
    adj_p(adj_p>0.05) = 1;
    vox_res(idx,3) = adj_p;
    
    [~,~,~, adj_p] = fdr_bh(pf);
    adj_p(adj_p>0.05) = 1;
    fvox_res(fidx,3) = adj_p;
    
    % Reshape and save to file
    vox_res = reshape(vox_res,[img_sizes{i},3]);
    fvox_res = reshape(fvox_res,[1360,1360,3]);    
    
    % Save file
    fname = fullfile(config.results_directory,sprintf("%s_%s_voxels.nii",...
        config.prefix,config.class_names(i)));
    niftiwrite(vox_res,fname)
        
    fname = fullfile(config.results_directory,sprintf("%s_%s_flatmap.nii",...
        config.prefix,config.class_names(i)));
    niftiwrite(fvox_res,fname)
end

end