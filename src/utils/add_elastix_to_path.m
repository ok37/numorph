function add_elastix_to_path(elastix_path_bin,elastix_path_lib)

% Add elastix paths if present
if exist('elastix_path_bin','var')
    path1 = getenv('PATH');
    if ~contains(path1,'elastix')
        path1 = [path1, ':', elastix_path_bin];
        setenv('PATH',path1)
    end
end
if exist('elastix_path_lib','var')
    if ismac
        path2 = getenv('DYLD_LIBRARY_PATH');
        if ~contains(path2,'elastix')
            path2 = [path2, ':', elastix_path_lib];
            setenv('DYLD_LIBRARY_PATH',path2)
        end
    else
        path2 = getenv('LD_LIBRARY_PATH');
        if ~contains(path2,'elastix')
            path2 = [path2, ':', elastix_path_lib];
            setenv('LD_LIBRARY_PATH',path2)
        end
    end
end

end