function vs = create_stats_volume2(df_stats, markers, structure_depth, compare_structures_by, orientation)
% Create a new annotation volume where each structure index is colored with
% a new field from the summary stats table
av = load('./supplementary_data/annotationData.mat');
downsample_factor = 0.1;

% Permute the annotation volume based on desired brain orientation
if isequal(orientation,"sagittal")
    av.annotationVolume = permute(av.annotationVolume,[1,3,2]);
elseif isequal(orientation, "axial")
    av.annotationVolume = permute(av.annotationVolume,[3,2,1]);
end

% Subset slices to reduce data sizes and make things smoother
[nrows,ncols,nslices] = size(av.annotationVolume);
nslices = round(nslices*downsample_factor);
av.annotationVolume = imresize3(av.annotationVolume,[nrows,ncols,nslices],'Method','nearest');

% Calculate number of markers
nmarkers = length(markers);

% Begin building visualization structure
% Add structure info + initial slice positions
vs.info = table2struct(df_stats(:,1:9));
vs.depth = structure_depth;
vs.z = round(size(av.annotationVolume,3)/2);
vs.bregma = (540/1320)*downsample_factor;
vs.marker_names = markers;

% Add colors with indexes for major structure divisions in Allen reference
c = load(fullfile('./visualization','allen_ccf_colormap_2017.mat'));
vs.av_colors = c.cmap;
vs.av_cmap_labels = ["Isocortex","OLF","HPF","CTXsp","STR","PAL","TH","HY","MB","HB","CB","FT"];
vs.av_cmap_pos = [5,379,454,555,573,608,641,715,806,882,1014,1101];

% Find structures to keep based on depth limit
% Reword this, slow
if isequal(compare_structures_by,'depth') && ~isempty(structure_depth)
    rm_idx = df_stats.depth > structure_depth;
    deep = df_stats(rm_idx,:);
    a = cellfun(@(s) strsplit(s,'/'),deep.structure_id_path,'UniformOutput',false);
    new_id = cellfun(@(s) str2double(s{structure_depth+2}),a);
    for i = 1:length(new_id)
        new_index(i) = df_stats.index(df_stats.id == new_id(i)); 
    end
    old_index = deep.index;
    for i = 1:length(old_index)
        av.annotationVolume(av.annotationVolume == old_index(i)) = new_index(i);
    end
    % Potentially edit
    vs.av = av.annotationVolume+1;
    
elseif isequal(compare_structures_by,'csv')
    [av.annotationVolume,av.annotationIndexes] = ...
        bin_annotation_structures(av.annotationVolume,'cortex','true');
    
    binnedMask = ismember(av.annotationVolume,av.annotationIndexes);
    
    % Potentially edit
    vs.av = av.annotationVolume+1;
    vs.binMask = binnedMask;
end 

% Precompute annotation indexes for each slice
% Add boundary image for each slice
% Add AP position of Bregma
vs.boundaries = zeros(size(av.annotationVolume),'logical');
for i = 1:nslices
    vs.av_idx(i) = {unique(vs.av(:,:,i))};
    vs.ap(i) = ((540*downsample_factor) - i)/(100*downsample_factor);
    vs.boundaries(:,:,i) = boundarymask(av.annotationVolume(:,:,i),8);
end

vs.bin_idx = cellfun(@(s) ismember(s-1,av.annotationIndexes),vs.av_idx,...
    'UniformOutput',false);

% Add calculated statistics for each marker
stats_columns = table2array(df_stats(:,10:end));
stats_headers = df_stats.Properties.VariableNames(10:end);

vs.indexes = df_stats(stats_columns(:,1)>0,:).index+1;

% Now for each marker, add counts and densities
for i = 1:nmarkers    
    % First add volume only to the first channel
    if i == 1
        if any(contains(stats_headers,"Volume"))
            volume_headers = stats_headers(contains(stats_headers,"Volume"));
            vs.markers.volumes.title = string(strrep(volume_headers,'_',' '));
            vs.markers(i).volumes.name = "Volume";
            
            values = stats_columns(:,contains(stats_headers,"Volume"));
            
            vs.markers(i).volumes.values(:,1:5) = values(:,1:5);
            vs.markers(i).volumes.values(:,6:7) = -log10(values(:,6:7));
            vs.markers(i).volumes.values(:,8) = values(:,8);
            
            values = values(vs.indexes,:);

            vs.markers(i).volumes.cmin = min(values);
            vs.markers(i).volumes.cmin(6:8) = 0;
            vs.markers(i).volumes.cmax = max(values);
            vs.markers(i).volumes.cmax(6:7) = max(-log10(values(:,6:7)));
            vs.markers(i).volumes.cmax(8) = 4;  
            
            pchange_max = max(abs(vs.markers(i).volumes.cmin(5)),vs.markers(i).volumes.cmax(5))*1.1;
            vs.markers(i).volumes.cmin(5) = -pchange_max;
            vs.markers(i).volumes.cmax(5) = pchange_max;
            
        end
    end

    if any(contains(stats_headers,"Counts"))
        count_headers = stats_headers(contains(stats_headers,"Counts")...
            & contains(stats_headers,markers(i)));
        vs.markers(i).counts.title = string(strrep(count_headers,'_',' '));
        vs.markers(i).counts.name = markers(i);

        values = stats_columns(:,contains(stats_headers,"Counts")...
            & contains(stats_headers,markers(i)));
        
        vs.markers(i).counts.values(:,1:5) = values(:,1:5);
        vs.markers(i).counts.values(:,6:7) = -log10(values(:,6:7));
        vs.markers(i).counts.values(:,8) = values(:,8);
            
        values = values(vs.indexes,:);

        vs.markers(i).counts.cmin = min(values);
        vs.markers(i).counts.cmin(6:8) = 0;
        vs.markers(i).counts.cmax = max(values);
        vs.markers(i).counts.cmax(6:7) = max(-log10(values(:,6:7)));
        vs.markers(i).counts.cmax(8) = 4;
        
        pchange_max = max(abs(vs.markers(i).counts.cmin(5)),vs.markers(i).counts.cmax(5))*1.1;
        vs.markers(i).counts.cmin(5) = -pchange_max;
        vs.markers(i).counts.cmax(5) = pchange_max;

    end
    
    if any(contains(stats_headers,"Density"))
                vs.markers(i).density.name = markers(i);
        count_headers = stats_headers(contains(stats_headers,"Density")...
            & contains(stats_headers,markers(i)));
        vs.markers(i).density.title = string(strrep(count_headers,'_',' '));
        vs.markers(i).density.name = markers(i);

        values = stats_columns(:,contains(stats_headers,"Density")...
            & contains(stats_headers,markers(i)));
        
        vs.markers(i).density.values(:,1:5) = values(:,1:5);
        vs.markers(i).density.values(:,6:7) = -log10(values(:,6:7));
        vs.markers(i).density.values(:,8) = values(:,8);
        
        values = values(vs.indexes,:);
        
        vs.markers(i).density.cmin = min(values);
        vs.markers(i).density.cmin(6:8) = 0;
        vs.markers(i).density.cmax(:,6:7) = max(-log10(values(:,6:7)));
        vs.markers(i).density.cmax(8) = 4;  

        pchange_max = max(abs(vs.markers(i).density.cmin(5)),vs.markers(i).density.cmax(5))*1.1;
        vs.markers(i).density.cmin(5) = -pchange_max;
        vs.markers(i).density.cmax(5) = pchange_max;
    end
    
    % Set colors
    brew1 = brewermap(255,'YlOrBr');
    brew1 = [1,1,1;brew1(2:end,:)];
    vs.markers(i).colors(1:4) = {brew1};
    
    brew2 = brewermap(255,'*RdBu');
    brew2(123,:) = [1 1 1];
    vs.markers(i).colors(5) = {brew2};
    
    brew3 = brewermap(255,'*RdBu');
    brew3(123,:) = [1 1 1];
    vs.markers(i).colors(6:7) = {brew3};
    
    brew4 = brewermap(5,'Purples');
    brew4 = [1,1,1;brew4(2:end,:)];
    vs.markers(i).colors(8) = {brew4};
    
    % 1.30103,2,3,5  
    
    vs.markers(i).log(1:4) = true;
    vs.markers(i).log(5:8) = false;
end


end