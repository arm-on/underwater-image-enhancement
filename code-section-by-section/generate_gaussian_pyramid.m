function gaussian_pyramid = generate_gaussian_pyramid(img, num_levels)
%{
    Description: Given an input image, the function generates a gaussian
    pyramid of the desired levels.

    Input:
        - img: an image
        - num_levels: the desired number of levels present in the pyramid

    Output:
        - gaussian_pyramid: the generated Gaussian pyramid
%}
    % defining the pyramid
    gaussian_pyramid = cell(num_levels,1);
    
    for i=1:num_levels
        if i==1
            gaussian_pyramid{i} = img;
        else
            gaussian_pyramid{i} = impyramid(gaussian_pyramid{i-1},'reduce');
        end
    end
end