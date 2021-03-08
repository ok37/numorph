function I_final = create_image_slice(centroids, path_table_stitched, config, z_position)
% Create overlay of centroids with cell-type

if nargin < 4
    z_position = 700;
end

z_position = z_position-1:z_position+1;

centroids(:,1:3) = centroids(:,1:3)+1;
idx = any(centroids(:,3) == z_position,2);
centroids = centroids(idx,:);

pos = centroids(:,1:3);
val = centroids(:,end);

tempI = imread(path_table_stitched.file{1});
dims = round(size(tempI));

I_final = zeros([dims(1), dims(2), 3],'uint8');
se = strel([0 1 0;1 1 1;0 1 0]);

for z = 1:3
    R = zeros(size(tempI),'uint8');
    G = zeros(size(tempI),'uint8');
    B = zeros(size(tempI),'uint8');
    
    pos1 = pos(pos(:,3) == z_position(z),:);
    val1 = val(pos(:,3) == z_position(z),:);
    
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
    if z == 2
        R = imdilate(R,se);
        G = imdilate(G,se);
        B = imdilate(B,se);
    end
    
    if i == 1
        I_final = cat(3,R,G,B);
    else
        I_final = I_final + cat(3,R,G,B);
    end
end

fprintf('%s\t Writing .png image \n',datetime('now'))
save_name = sprintf('%s_slice_%d.png',config.sample_name,z_position(2));
imwrite(I_final,save_name)

end