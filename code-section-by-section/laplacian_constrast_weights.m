function weights = laplacian_constrast_weights(img)
%{
    Description: The function estimates the global constrast by computing
    the absolute value of a Laplacian filter applied on each input
    luminance channel. Using this procedure, a value will be associated
    with each pixel position which will be called the Laplacian contrast
    weight.

    Input:
        - img: an RGB image (a mxnx3 tensor)
    
    Output:
        - weights: the weights defined by the Laplacian filter (a mxnx3
        tensor)
%}
    laplacian_filter = [0,1,0;1,-4,1;0,1,0];
    lab_img = rgb2lab(img);
    luminance_channel = lab_img(:,:,1);
    weights = imfilter(luminance_channel, laplacian_filter);
end