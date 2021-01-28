function save_samples(config, process, path_table)
%--------------------------------------------------------------------------
% Save samples results during image processing steps.
%--------------------------------------------------------------------------

% Defaults
spacing = [3,3,20]; % Downsampling for alignment check

% Read image filename information if not provided
if nargin<3
    path_table = path_to_table(config, [], false);
end
        
% Create samples directory
samples_directory = fullfile(config.output_directory,'samples');
if exist(samples_directory,'dir') ~= 7
    mkdir(samples_directory);
end

switch process
    case 'intensity'
        % Check for adjustment parameters or try loading
        if ~isfield(config,'adj_params') || isempty(config.adj_params)
            try
                load(fullfile(config.output_directory,'variables','adj_params.mat'),'adj_params')
                [adj_params, config] = check_adj_parameters(adj_params,config);
                config.adj_params = adj_params;
            catch
                error("No adjustment parameters found in config structure "+...
                    "or in the output directory")
            end
        end
        
        % Create sub-directory
        sub_directory = fullfile(samples_directory,'intensity_adjustment');
        if exist(sub_directory,'dir') ~= 7
            mkdir(sub_directory);
        end
        
        % Create image
        % Take middle z position
        path_table = path_table(path_table.z == median(path_table.z),:);
        fprintf("%s\t Saving adjusted images for middle slice \n",datetime('now'));
        for i = 1:height(path_table)            
            % Read image
            img = imread(path_table(i,:).file{1});
            r = path_table(i,:).y;
            c = path_table(i,:).x;
            z = path_table(i,:).z;
            c_idx = path_table(i,:).channel_num;
            
            % Apply intensity adjustments
            img_adj = apply_intensity_adjustment(img,'params',...
                config.adj_params.(config.markers(c_idx)),...
                    'r',r,'c',c);

            % Apply any type of post-processing
            img_adj = postprocess_image(config, img_adj, c_idx);
            
            % Concatenate images
            img = horzcat(img,img_adj);
            
            % Save results
            file_name = sprintf('%s_%d_%d_%d.tif',path_table(i,:).markers,r,c,z);
            imwrite(img,fullfile(sub_directory,file_name))
        end
        
        % Check for exporgraphics function in MATLAB 2020
        if exist('exportgraphics','file') == 0
                warning("Could not save flatfield images due to missing "+...
                    "export function. Update to MATLAB 2020.")
            return
        end
        
        % Save shading correction images
        fprintf("%s\t Saving current intensity adjustment visualizations \n",datetime('now'));
        markers = fields(config.adj_params);
        for i = 1:length(markers)
            params = config.adj_params.(markers{i});
            
            % Save flatfield 
            flatfield = params.flatfield;
            fig = figure('visible','off');
            imagesc(flatfield,[0.5,1.5])
            colorbar
            axis image
            exportgraphics(fig,fullfile(sub_directory,...
                sprintf('flatfield_%d.png',i)))
            
            % Save y_adj 
            y_adj = params.y_adj;
            y_adj = repmat(y_adj,1,size(flatfield,2));
            fig = figure('visible','off');
            imagesc(y_adj,[0.5,1.5])
            colorbar
            axis image
            exportgraphics(fig,fullfile(sub_directory,...
                sprintf('y_adj_%d.png',i)))
            
            % Save t_adj
            t_adj = params.t_adj;
            if t_adj~=1
                fig = figure('visible','off');
                heatmap(t_adj)
                caxis([0.5,1.5])
                colormap parula
                exportgraphics(fig,fullfile(sub_directory,...
                    sprintf('tile_adj_%d.png',i)))
            end
        end
    case 'alignment'
        % Run alignment check for each tile position present
        x = unique(path_table.x);
        y = unique(path_table.y);
        n_tiles = length(x)*length(y);
        
        if n_tiles == 1
            check_alignment(config, {y, x}, [], spacing);
        else
            spacing(3) = spacing(3)*5;
            
            I_adj2 = cell(length(y),1);
            for i = 1:length(y)
                I_adj = cell(1,length(x));
                for j = 1:length(x)
                    I_adj{j} = check_alignment(config, {y(i), x(j)}, [], spacing);
                end
                I_adj2{i} = cat(2,I_adj{:});
            end
            I_adj = cat(1,I_adj2{:});
            
            % Save multi-tile
            save_directory = fullfile(config.output_directory,'samples','alignment');
            for i = 1:size(I_adj,4)
                fname = fullfile(save_directory,sprintf('%s_full.tif',config.markers(i)));
                options.overwrite = true;
                options.message = false;
                saveastiff(squeeze(I_adj(:,:,:,i)),char(fname),options);
            end
        end
end

end