function add_elastix_to_path(elastix_path_bin,elastix_path_lib)

% Add elastix paths if present
path1 = getenv('PATH');
if ~contains(path1,elastix_path_bin)
    disp('hello')
    path1 = [elastix_path_bin, ':', path1];
    setenv('PATH',path1)
end

if ismac
    path2 = getenv('DYLD_LIBRARY_PATH');
    if ~contains(path2,elastix_path_lib)
        path2 = [elastix_path_lib, ':', path2];
        setenv('DYLD_LIBRARY_PATH',path2)
    end
else
    path2 = getenv('LD_LIBRARY_PATH');
    if ~contains(path2,elastix_path_lib)
        path2 = [elastix_path_lib, ':', path2];
        setenv('LD_LIBRARY_PATH',path2)
    end
end

end