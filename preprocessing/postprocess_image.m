function I = postprocess_image(config, I, marker_idx)
%--------------------------------------------------------------------------
% Postprocess an image by filtering, subtracting background, or enhancing
% blobs as specified in the config structure.
%--------------------------------------------------------------------------

% Rescale intensities
if isequal(config.rescale_intensities(marker_idx),"true")
    I = apply_intensity_adjustment(I,'l_thresh',config.adj_params.lowerThresh(marker_idx),...
        'u_thresh',config.adj_params.upperThresh(marker_idx),'gamma',config.adj_params.gamma(marker_idx));
end

% Apply background subtraction
if isequal(config.subtract_background(marker_idx),"true")
    I = smooth_background_subtraction(I, 'false', config.nuc_radius);
end

% Apply Difference-of-Gaussian filter
if isequal(config.DoG_img(marker_idx),"true")
    I = dog_adjust(I,config.nuc_radius,config.DoG_minmax,config.DoG_factor(marker_idx));                
end

% Apply smoothing filter
if ~isequal(config.smooth_img(marker_idx),"false")
    I = smooth_image(I, config.smooth_img(marker_idx),config.smooth_sigma(marker_idx));                
end

end