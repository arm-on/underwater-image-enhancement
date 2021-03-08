function saliency_weights = compute_saliency_weights(img)
%{
    Description: The function computes the saliency weights for a given
    image. These weights aim at emphasizing the salient objects that lose
    their prominence in the underwater scene. The formula is written below:
    
    W_S(w,h) = mean_pixel_value - gaussian_filter_img(x,y)

    Input: 
        - img: an RGB image (a mxnx3 tensor)

    Output:
        - saliency_weights: the computed weights for each pixel
%}
    saliency_weights = zeros(height(img), width(img));
    
    lab_img = rgb2lab(img);
    
    % decomposing the LAB image to its three channels
    lab_l = lab_img(:,:,1);
    lab_a = lab_img(:,:,2);
    lab_b = lab_img(:,:,3);
    
    % measuring the mean of each channel of the LAB image
    mean_l = mean(lab_l(:));
    mean_a = mean(lab_a(:));
    mean_b = mean(lab_b(:));
    
    % applying a gaussian filter on the image, and calculting its LAB
    % version
    gaussian_filtered_img = imgaussfilt(img);
    lab_gaussian_img = rgb2lab(gaussian_filtered_img);
    
    % calculating the weight of each position using the formula written in
    % the description of the function
    for h=1:height(img)
        for w=1:width(img)
            
            curr_gauss_l = lab_gaussian_img(h,w,1);
            curr_gauss_a = lab_gaussian_img(h,w,2);
            curr_gauss_b = lab_gaussian_img(h,w,3);
            
            curr_weight = sqrt((mean_l-curr_gauss_l)^2+...
                (mean_a-curr_gauss_a)^2+...
                (mean_b-curr_gauss_b)^2);
            
            saliency_weights(h,w) = curr_weight;
        end
    end
end

