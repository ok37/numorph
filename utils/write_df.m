function write_df(df,file_path,overwrite_flag,matrix_or_table)
%--------------------------------------------------------------------------
% Save centroid list or data matrix as new file or overwrite existing file.
%--------------------------------------------------------------------------

if nargin<3
    matrix_or_table = 'matrix';
elseif ~isequal(matrix_or_table,'matrix') || ~isequal(matrix_or_table,'table')
    error("Enter whether to save as matrix or table")
end

% If overwriting, just save file
if overwrite_flag
   if isequal(matrix_or_table,'matrix')
       writematrix(df, file_path)
   else
        writetable(df, file_path)
   end
   return
end

% Check if file exists
if exist(file_path,'file') == 2
    exists = true;
else
    exists = false;
end




end