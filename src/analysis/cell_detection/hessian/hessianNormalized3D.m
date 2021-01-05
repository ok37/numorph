function H = hessianNormalized3D(img, s1, s2, s3)
% hessianNormalized2D gives the normalized hessian of an image using
% DIPimage
    img = imgaussfilt3(img,0.3);
    H = Hessian3D(img);

    % Normalize the hessian
    H{1,1} = s1^2 * H{1,1};
    H{1,2} = s1*s2 * H{1,2};
    H{2,1} = s1*s2 * H{2,1};
    H{2,2} = s2^2 * H{2,2};
    
    H{1,3} = s1^s3 * H{1,3};
    H{2,3} = s2*s3 * H{2,3};
    H{3,1} = s3*s1 * H{3,1};
    H{3,2} = s3*s2 * H{3,2};
    H{3,3} = s3^2 * H{3,3};

end