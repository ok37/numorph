function df = trim_to_deepest_structure(df)

ids = df.id;
a = cellfun(@(s) strsplit(s,'/'),df.structure_id_path,'UniformOutput',false);
a = cellfun(@(s) rmmissing(str2double(s)),a,'UniformOutput',false);

% Count id instances in structure tree 
x = zeros(1,length(ids));
for i = 1:length(x)
    x(i) = sum(cellfun(@(s) any(ismember(s,ids(i))),a));
end

% Remove rows id appears more than once
df = df(x==1,:);

end