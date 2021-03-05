function vs = create_stats_volume2(df_stats, markers, orientation, custom_comp)
% Create a new annotation volume where each structure index is colored with
% a new field from the summary stats table
home_path = fileparts(which('NM_config.m'));
av = load(fullfile(home_path,'data','annotationData.mat'));
downsample_factor = 0.1;

% Permute the annotation volume based on desired brain orientation
if isequal(orientation,"sagittal")
    av.annotationVolume = permute(av.annotationVolume,[1,3,2]);
elseif isequal(orientation, "axial")
    av.annotationVolume = permute(av.annotationVolume,[3,2,1]);
end

% Add annotation info for all structures
temp_tbl = readtable(fullfile(home_path,'annotations','structure_template.csv'));
vs.info = table2struct(temp_tbl(:,1:9));

% Subset slices to reduce data sizes and make things smoother
[nrows,ncols,nslices] = size(av.annotationVolume);
nslices = round(nslices*downsample_factor);
av.annotationVolume = imresize3(av.annotationVolume,[nrows,ncols,nslices],'Method','nearest');

% Calculate number of markers
nmarkers = length(markers);

% Begin building visualization structure
% Add structure info + initial slice positions
vs.z = round(size(av.annotationVolume,3)/2);
vs.bregma = (540/1320)*downsample_factor;
vs.marker_names = markers;

% Add colors with indexes for major structure divisions in Allen reference
c = load(fullfile(home_path,'src','utils','allen_ccf_colormap_2017.mat'));
vs.av_colors = c.cmap;
vs.av_cmap_labels = ["Isocortex","OLF","HPF","CTXsp","STR","PAL","TH","HY","MB","HB","CB","FT"];
vs.av_cmap_pos = [5,379,454,555,573,608,641,715,806,882,1014,1101];

% Bin to correct structures
df_stats = trim_to_deepest_structure(df_stats);
[aV,aI] = bin_annotation_structures(av.annotationVolume,df_stats,'true');
binnedMask = ismember(aV,aI);
% Potentially edit
vs.av = aV+1;
vs.indexes = aI;
vs.binMask = binnedMask;

% Precompute annotation indexes for each slice
% Add boundary image for each slice
% Add AP position of Bregma
vs.boundaries = zeros(size(av.annotationVolume),'logical');
vs.av_stat_idx = cell(1,nslices);
for i = 1:nslices
    vs.av_idx(i) = {unique(vs.av(:,:,i))};
    vs.av_stat_idx{i} = find(ismember(aI+1,vs.av_idx{i}));
    vs.boundaries(:,:,i) = boundarymask(aV(:,:,i),8);
    if isequal(orientation,"coronal")
        vs.ap(i) = ((540*downsample_factor) - i)/(100*downsample_factor);
    end
end

vs.bin_idx = cellfun(@(s) ismember(s-1,av.annotationIndexes),vs.av_idx,...
    'UniformOutput',false);

% Now for each marker, add counts and densities
stats_headers = df_stats.Properties.VariableNames(10:end);
for i = 1:nmarkers    
    % First add volume only to the first channel
    if i == 1 && any(contains(stats_headers,"Volume"))
        vs = add_stats_to_vs(vs,df_stats,markers(i),"Volume",i);
    end
    vs = add_stats_to_vs(vs,df_stats,markers(i),"Counts",i);
    vs = add_stats_to_vs(vs,df_stats,markers(i),"Density",i);
end

for i = 1:length(custom_comp)
    vs = add_stats_to_vs(vs,df_stats,custom_comp(i),"Custom",i+nmarkers);
end


% Update indexes to match volume
vs.indexes = vs.indexes+1;

% Save volume
if isequal(orientation,"coronal")
    save_path = fullfile(fileparts(which('NM_config')),'data','coronalSliceVol.mat');
elseif isequal(orientation,"sagittal")
    save_path = fullfile(fileparts(which('NM_config')),'data','sagittalSliceVol.mat');
elseif isequal(orientation,"axial")
    save_path = fullfile(fileparts(which('NM_config')),'data','axialSliceVol.mat');
end
save(save_path,'vs')

end


function vs = add_stats_to_vs(vs,df_stats,marker,colname,idx)

% Add calculated statistics for each marker
stats_columns = table2array(df_stats(:,10:end));
stats_headers = df_stats.Properties.VariableNames(10:end);

header = stats_headers(contains(stats_headers,colname));
vs.markers(idx).(colname).title = string(strrep(header,'_',' '));
vs.markers(idx).(colname).name = colname;

if isequal(colname,"Volume")
    values = stats_columns(:,contains(stats_headers,colname));
elseif isequal(colname,"Custom")
    values = stats_columns(:,contains(stats_headers,colname)...
        & contains(stats_headers,marker));
else
    values = stats_columns(:,contains(stats_headers,colname)...
            & contains(stats_headers,marker));
end

vs.markers(idx).(colname).values(:,1:5) = values(:,1:5);
vs.markers(idx).(colname).values(:,6:7) = -log10(values(:,6:7));
vs.markers(idx).(colname).values(:,8) = values(:,8);

vs.markers(idx).(colname).cmin = min(values);
vs.markers(idx).(colname).cmin(6:8) = 0;
vs.markers(idx).(colname).cmax = max(values);
vs.markers(idx).(colname).cmax(6:7) = max(-log10(values(:,6:7)));
vs.markers(idx).(colname).cmax(8) = 4;  

pchange_max = max(abs(vs.markers(idx).(colname).cmin(5)),vs.markers(idx).(colname).cmax(5))*1.1;
vs.markers(idx).(colname).cmin(5) = -pchange_max;
vs.markers(idx).(colname).cmax(5) = pchange_max;

% Set colors
brew1 = brewermap(201,'YlOrBr');
vs.markers(idx).colors(1:4) = {brew1};
    
brew2 = brewermap(201,'*RdBu');
vs.markers(idx).colors(5) = {brew2};
    
brew3 = brewermap(201,'BuPu');
vs.markers(idx).colors(6:7) = {brew3};
    
brew4 = brewermap(5,'Purples');
vs.markers(idx).colors(8) = {brew4};
    
end
