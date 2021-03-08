function laplacian_pyramid = generate_laplacian_pyramid(gaussian_pyramid)
%{
    Description: Given the Gaussian pyramid of an image, the function
    generates the corresponding Laplacian pyramid by subtracting each level
    of the Gaussian pyramid from the upsampled version of the upper Gaussian level.

    Input: 
        - gaussian_pyramid: a cell array containing the Gaussian pyramid of
        an image

    Output:
        - laplacian_pyramid: the generated Laplacian pyramid
%}
    % calculating the number of levels
    num_levels = size(gaussian_pyramid, 1);
    
    % defining the Laplacian pyramid
    laplacian_pyramid = cell(num_levels, 1);
    laplacian_pyramid{num_levels} = gaussian_pyramid{num_levels};
    
    for i=1:num_levels-1
        upsampled_upper_level = impyramid(gaussian_pyramid{i+1},'expand');
        desired_height = height(gaussian_pyramid{i});
        desired_width = width(gaussian_pyramid{i});
        upsampled_upper_level = imresize(upsampled_upper_level, [desired_height, desired_width]);
        laplacian_pyramid{i} = gaussian_pyramid{i}-upsampled_upper_level;
    end
    
end
