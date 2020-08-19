function [pc_img,ref_img,tform,type,cc] = calculate_phase_correlation(mov_img,ref_img,peaks,usfac)
warning('off','images:imregcorr:weakPeakCorrelation')

%%%%%%%Important: Any predicted translation above this threshold in either
%%%%%%%dimension will be considered invalid and the transform will be
%%%%%%%thrown out
shift_threshold = 120;


%Set precision as 1 pixel unless specified
if isempty(usfac)
    usfac = 1;
end

%Set number of peaks as 1 unless specified
if isempty(peaks)
    peaks = 1;
end

%Perform subpixel registration using phase correlation. 
%Output is y_translation;x_translation
output = dftregistration(ref_img,mov_img,peaks,usfac);
type = 1;

if isempty(output)
%If pc using dft doesn't give anything useful, try MATLAB's imregcorr. This
%can sometimes give something using useful when SNR in one of the images is
%low. However only top peak is chosen
    tform = imregcorr(mov_img,ref_img,'translation','Window',false);
    output = [tform.T(6); tform.T(3)];
    type = 2; 
    %Remove if contains large shifts
    if output(1)>shift_threshold || output(2)>shift_threshold
        output = [];
    end
    %Try first without windowing. Then try with windowing
    if isempty(output)
        tform = imregcorr(mov_img,ref_img,'translation','Window',true);
        output = [tform.T(6); tform.T(3)];
        type = 3; 
        %Remove if contains large shifts
        if output(1)>shift_threshold || output(2)>shift_threshold
            output = [];
        end
    end
end

%If no result, give empty transform
if isempty(output)
    pc_img = mov_img;
    tform = [];
    
    ref_img(pc_img == 0) = 0;
    pc_img(ref_img == 0) = 0;
    
    cc = corr2(ref_img(ref_img>0),pc_img(pc_img>0));
    type = 0;
else
    %Translate image
    pc_img = imtranslate(mov_img,[output(2) output(1)]);

    %Save translation as affine2d object
    tform = affine2d([1 0 0; 0 1 0; output(2) output(1) 1]);
    
    %Calculate cross-correlation
    cc = corr2(ref_img,pc_img);
end
end