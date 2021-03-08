function white_balanced_img = apply_gray_world(img, percentiles)
%{
    Description: The function applies the Gray World algorithm on a given
    image, and returns the result.

    Input:
        - img: an RGB image (a mxnx3 tensor)
        - percentiles: the percentage of the pixel values from bottom and
        top not to be considered in the algorithm

    Output:
        - white_balanced_img: The result of performing the Gray World
        algorithm on the given image
%}
    linear_img = rgb2lin(img);
    scene_illumination_estimate = illumgray(linear_img,percentiles);
    linear_white_balanced_img = chromadapt(linear_img, scene_illumination_estimate, 'ColorSpace', 'linear-rgb');
    white_balanced_img = lin2rgb(linear_white_balanced_img);
end