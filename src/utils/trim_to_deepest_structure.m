function df = trim_to_deepest_structure(df, idxs)

trim = false;
if nargin<2
    trim = true;
end

if trim 
    ids = df.info.id;
    a = cellfun(@(s) strsplit(s,'/'),df.info.structure_id_path,'UniformOutput',false);
    a = cellfun(@(s) rmmissing(str2double(s)),a,'UniformOutput',false);
    
    % Count id instances in structure tree 
    x = zeros(1,length(ids));
    for i = 1:length(x)
        x(i) = sum(cellfun(@(s) any(ismember(s,ids(i))),a));
    end
    
    % Remove rows id appears more than once
    df.info = df.info(x==1,:);
    for i = 1:length(df.stats)
        df.stats(i).stats = df.stats(i).stats(x==1,:);
    end
else
    % Just select annotations that are specified
    idxs_native = df.info.index;
    keep_idx = ismember(idxs_native, idxs);
    for i = 1:length(df.stats)
        df.stats(i).stats = df.stats(i).stats(keep_idx,:);
    end
    
end

end