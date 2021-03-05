function conda_path = add_conda_to_path

conda_path = [];

% First check if conda binary already in PATH
PATH = getenv('PATH');
if contains(PATH,'conda3')
    a = strsplit(PATH,':');
    idx = cellfun(@(s) contains(s,'conda3'),a);    
    conda_path = a{find(idx,1)};
    return
end

% Scan for location of conda binary
if  ismac || isunix       
    userDir = char(java.lang.System.getProperty('user.home'));
    locations = {'.bash_profile','.bashrc','.zshrc','.zprofile'};
    
    for i = 1:length(locations)
        if isfile(fullfile(userDir,locations{i}))
            fid = fopen('~/.bash_profile');
            c = textscan(fid,'%s');
            fclose(fid);
            
            if any(contains(c{:},{'conda3/bin:$PATH'}))
                idx = find(contains(c{:},{'conda3/bin:$PATH'}),1);
                s1 = regexp(c{1}{idx,:},'/');
                s2 = regexp(c{1}{idx,:},':');
                
                conda_path = c{1}{idx,:}(s1(1):s2(1)-1);
            end
        end
    end
end

% Add to PATH
if ~isempty(conda_path)
    PATH = conda_path + ":" + PATH; 
    setenv('PATH',PATH)
end

end