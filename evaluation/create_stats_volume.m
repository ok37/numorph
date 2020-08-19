function [A,V] = create_stats_volume(df_stats,column)
% Create a new annotation volume where each structure index is colored with
% a new field from the summary stats table
av = load('./supplementary_files/annotationData.mat');

downsample_factor = 0.1;

[nrows,ncols,nslices] = size(av.annotationVolume);
nslices = round(nslices*downsample_factor);

av.annotationVolume = imresize3(av.annotationVolume,[nrows,ncols,nslices],'Method','nearest');

column_data = table2array(df_stats(:,9+column));

% Data scaling to 255
min_val = min(column_data(av.annotationIndexes));
max_val = max(column_data(av.annotationIndexes));

% Automatically convert to log10 if max value is >1000
if log10(max_val - min_val) > 3
    column_data = log10(column_data);
    column_data(column_data == -Inf) = 0;
    max_val = log10(max_val);
end

s = 255/(max_val-min_val);

V = zeros(size(av.annotationVolume),'uint8');
for i = 1:length(av.annotationIndexes)
    V(av.annotationVolume == av.annotationIndexes(i)) = round((column_data(av.annotationIndexes(i))-min_val)*s);
end

A = av.annotationVolume;
end