points = readmatrix('native_landmarks110_spread.csv');
points2 = readmatrix('native_landmarks11.csv');

%% Check points distance
points = points(:,6:8);

distance = zeros(1,size(points,1));
for i = 1:length(distance)
   p = points(i,:);
   d = sqrt((points(:,1)-p(1)).^2 + (points(:,2)-p(2)).^2 + (points(:,3)-p(3)).^2);
   distance(i) = min(d(d>0));
end


[sorted_distance,idx] = sort(distance);
p_idx = 1:200;

p_idx = p_idx(idx);
%% List of point indexes to replace
miss_indx = zeros(1,200);

for i = 1:size(points,1)
   if all(round(points(i,6:8)) ~= round(points2(i,6:8)))
        miss_indx(i) = 1;
   end
end

disp(sum(miss_indx(:)))






