function check_file_integrity(input_directory, chunk_size)
%--------------------------------------------------------------------------
% Check for corrupt images in directories by opening .tiffs using MATALB's
% imread function.
%--------------------------------------------------------------------------

if nargin<2
    chunk_size = 10;
end

dir_files = dir(input_directory);

% Find all tif files in directory and subdirectory
file_list = {};
for i = 1:length(dir_files)
   if isequal(dir_files(i).name,'.') || isequal(dir_files(i).name,'..')
       continue
   end
   
   if contains(dir_files(i).name,'.tif') || contains(dir_files(i).name,'.tiff')
       file_list = cat(1,file_list,{fullfile(input_directory,dir_files(i).name)});
       continue
   end
   
   if isfolder(fullfile(input_directory,dir_files(i).name))
      a = 1;
      
   end
   
    
end





end