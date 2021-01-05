a = dir(pwd);
pad = 4;


for i = 1:length(a)-2
    sub = load(a(i+2).name);
    centroid = sub.centroid;
    
    z_start = min(centroid(3,:))+pad/2;
    z_end = max(centroid(3,:))-pad/2;
    z_range = z_start:z_end;
    
    
    if i == 1
        centroid = centroid(:,ismember(centroid(3,:),z_range));
        centroids = centroid;

    else
        centroid = centroid(:,ismember(centroid(3,:),z_range));
        centroids = horzcat(centroids,centroid);
    end 
end