function sharpened_img = sharpen(img)
%{
    Description: The function first filters a given image using a Gaussian
    filter. Then, it computes a mask which by subtracting the filtered
    image from its original version . The result is called an unsharp mask.
    Afterwards the function adds the mask to the original image, and
    divides the result by two. The result of this step is returned as the
    sharpened version of the original image.

    Input:
        - img: an RGB image (a mxnx3 tensor)
    
    Output:
        - sharpened_img: the sharpened version of the given image (a mxnx3
        tensor)
%}
    gaussian_filtered_img = imgaussfilt(img);
    mask = img-gaussian_filtered_img;
    mask = histeq(mask);
    % mask = hist_stretch(mask, 0, 255); % this line makes some red artifacts
    % sharpened_img = (img + mask)/2;
    sharpened_img = (im2double(img) + im2double(mask))/2;
    sharpened_img = im2uint8(sharpened_img);
end

function hist_stretched_img = hist_stretch(img, min_val, max_val)
%{
    Description: The function adjusts an image by stretching its histogram
    and making it having more contrast.

    Input:
        - img: an RGB image
        - min_val: the desired minimum intensity after stretching
        - max_val: the desired maximum intensity after stretching

    Output:
        - hist_stretched_img: the result of stretching
%}
    curr_min_val = min(img(:));
    curr_max_val = max(img(:));
    
    hist_stretched_img = repmat(img, 1);
    for h=1:height(img)
        for w=1:width(img)
            hist_stretched_img(h, w) = (img(h, w) - curr_min_val)*((max_val - min_val)/(curr_max_val - curr_min_val)) + min_val;
        end
    end
end