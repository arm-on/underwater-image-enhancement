function uciqe_val = UCIQE(img)
%{
    Description: The function calculates the UCIQE metric for evaluating an
    underwater image.

    Input:
        - img: an RGB image
    Output:
        - uciqe_val: the UCIQE value assigned to the given image
%}
    img = im2double(img);
    
    c1 = 0.4680;
    c2 = 0.2745;
    c3 = 0.2576;
    
    lab_img = rgb2lab(img);
    hsv_img = rgb2hsv(img);
    
    saturation = hsv_img(:,:,2);
    luminance = lab_img(:,:,1);
    luminance = mat2gray(luminance);
    
    lab2lch_cform = makecform('lab2lch');
    lch_img = applycform(img, lab2lch_cform);
    chroma = lch_img(:,:,2);
    
    chroma_std = std(chroma(:));
    saturation_avg = mean(saturation(:));
    luminance_contrast = max(luminance(:))-min(luminance(:));
    
    uciqe_val = c1*chroma_std+c2*luminance_contrast+c3*saturation_avg;
   
end
