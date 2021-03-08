function I_final = create_image_slices(centroids, path_table_stitched, config)
% Create overlay of centroids with cell-type

spacing = 100;
s = 0.1;

z_positions = height(path_table_stitched)/length(config.markers);
z_positions = 100:spacing:z_positions;

centroids(:,1:3) = centroids(:,1:3)-1;
idx = any(centroids(:,3) == z_positions,2);
centroids = centroids(idx,:);

pos = centroids(:,1:3);
val = centroids(:,end);

tempI = imread(path_table_stitched.file{1});
dims = round(size(tempI)*s);

I_final = zeros([dims(1), dims(2), 3, length(z_positions)],'uint8');
se = strel('disk',10);

for z = 1:length(z_positions)
    R = zeros(size(tempI),'uint8');
    G = zeros(size(tempI),'uint8');
    B = zeros(size(tempI),'uint8');
    
    
    pos1 = pos(pos(:,3) == z_positions(z),:);
    val1 = val(pos(:,3) == z_positions(z),:);
    
    for i = 1:size(pos1,1)
        if val1(i) == 2
            G(pos1(i,1),pos1(i,2)) = 254;
        end
        if val1(i) == 3
            R(pos1(i,1),pos1(i,2)) = 254;
        end
        if val1(i) == 4
            B(pos1(i,1),pos1(i,2)) = 254;
        end
    end
    R = imdilate(R,se);
    G = imdilate(G,se);
    B = imdilate(B,se);
    
    R = imresize(R,[dims(1) dims(2)]);
    G = imresize(G,[dims(1) dims(2)]);
    B = imresize(B,[dims(1) dims(2)]);
        
    I_final(:,:,:,z) = cat(3,R,G,B);
end

fprintf('%s\t Writing .tif image \n',datetime('now'))
save_name = sprintf('%s_all_%d.tif',config.sample_name,spacing);
options.color = true;
options.overwrite = true;
options.verbose = false;
saveastiff(I_final,save_name,options)

end