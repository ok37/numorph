set1 = readmatrix('native_landmarks11.csv');
set1a = readmatrix('native_landmarks11_spread.csv');

ref_set1 = round(set1(:,6:8));
ref_set1a = round(set1a(:,6:8));


p_idx = find(all(ref_set1 ~= ref_set1a,2));

%%
new_points = readtable('native_landmarks11_spread.csv');
new_points = new_points(p_idx,:);
writetable(new_points,'new_points.csv')


%%
set2 = readtable('native_landmarks11.csv');
set3 = readtable('native_landmarks16.csv');
set4 = readtable('native_landmarks110.csv');







