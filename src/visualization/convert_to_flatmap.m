function fm = convert_to_flatmap(img)


home_path = fileparts(which('NM_config'));
fc_path = fullfile(home_path,'data','annotation_data','flatviewCortex.mat');
load(fc_path,'voxelmap')
load(fc_path,'voxelpaths')

if length(size(img)) == 4
    img = squeeze(img);
end

img = cat(3,img,zeros(size(img)));
img = permute(img,[3,2,1]);
img = flip(img,2);
img = flip(img,1);
results = zeros(1,length(voxelpaths));
for j = 1:length(voxelpaths)
    c = voxelpaths(:,j);
    results(j) = mean(img(c(c>1)));
end

fm = zeros(1,length(voxelmap(:)));
for j = 1:length(fm)
    if voxelmap(j)>1
        fm(j) = results(voxelmap(j));
    end
end

fm = reshape(fm,2720,1360);
fm = fm(1361:end,:);
fm = imrotate(fm,-90);
end