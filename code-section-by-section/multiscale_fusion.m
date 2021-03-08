function pyramid = multiscale_fusion(laplacian_pyramid_gc, laplacian_pyramid_sh, gaussian_pyramid_gc_weights, gaussian_pyramid_sh_weights)
%{
    Description: Given the Laplacian pyramids of the gamma corrected and the
    sharpened images along with the Gaussian pyramids of the normalized
    weights of them, the function applies the multiscale fusion.

    Input:
        - laplacian_pyramid_gc: the Laplacian pyramid of the gamma
        corrected image
        - laplacian_pyramid_sh: the Laplacian pyramid of the sharpened
        image
        - gaussian_pyramid_gc_weights: the Gaussian pyramid of the
        normalized weights of the gamma corrected image
        - gaussian_pyramid_sh_weights: the Gaussian pyramid of the
        normalized weights of the sharpened image

    Output:
        - pyramid: The pyramid resulted from the fusion process
%}
    
    % measuring the number of levels in the pyramid
    num_levels = size(laplacian_pyramid_gc, 1);
    
    % defining the resulting pyramid
    pyramid = cell(num_levels,1);
    
    % fusion
    for i=1:num_levels
        sharpened_sum = double(gaussian_pyramid_sh_weights{i}).*double(laplacian_pyramid_sh{i});
        gamma_corrected_sum = double(gaussian_pyramid_gc_weights{i}).*double(laplacian_pyramid_gc{i});
        pyramid{i} = sharpened_sum + gamma_corrected_sum;
    end
   
end