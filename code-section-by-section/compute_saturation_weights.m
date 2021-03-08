function saturation_weights = compute_saturation_weights(img)
%{
    Description: The function computes the saturation weights so that the
    fusion algorithm is able to adapt to chromatic information by
    advantaging highly saturated regions. The weight map is simply computed
    as the deviation (for every pixel location) between the R, G, and B
    color channels and the luminance of the input.

    The formula:
        W_sat = sqrt(  (1/3) * ((R-L)^2 + (G-L)^2 + (B-L)^2) )

    Input:
        - img: an RGB image (a mxnx3 tensor)

    Output:
        - saturation_weights: the computed weights
%}
    % converting the image to its LAB equivalent to get the luminance
    % channel
    lab_img = rgb2lab(img);
    luminance = lab_img(:,:,1);
    
    % computing the saturation weights
    saturation_weights = zeros(height(img), width(img));
    for h=1:height(img)
        for w=1:width(img)
            R = double(img(h,w,1));
            G = double(img(h,w,2));
            B = double(img(h,w,3));
            
            L = luminance(h,w);
            
            sum_of_differences = (R-L)^2 + (G-L)^2 + (B-L)^2;
            saturation_weights(h,w) = sqrt((1/3)*sum_of_differences);
        end
    end
    
end