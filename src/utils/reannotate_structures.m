% Create new table with containing the groups of interest
new_template_path = './supplementary_files/structure_template.csv';
df_stats = readtable(new_template_path);


I = nrrdread('annotation_10.nrrd');


%%
ids = df_stats.id;
indexes = df_stats.index;
indexes2 = unique(I(:));

I2 = zeros(size(I),'uint16');

%%

for i = 2:length(indexes2)
    % Check indexes is in id column
    row_index = find(indexes2(i) == ids);
    new_index = indexes(row_index);

    I2(I == indexes2(i)) = new_index;
end
niftiwrite(I2,'atlas.nii','Compressed',true)
%364