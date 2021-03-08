function eval = evaluate(img, ref_img)
%{
    Description: The function evaluates an image based on three measures:
    UCIQE, NIQE, and MEAN-CIE2000.

    Input: 
        - img: an RGB image
        - ref_img: another RGB image for comparison (used for CIE2000 evaluation)
    Output: 
        - eval: a vector containing the desired metrics
%}
    uciqe_val = UCIQE(img);
    niqe_val = niqe(img);
    cie2000_val = imcolordiff(ref_img,img,'Standard','CIEDE2000','kL',2,'K1',0.048,'K2',0.014);
    cie2000_val = double(mean(cie2000_val(:)));
    eval = cell2mat({uciqe_val, niqe_val, cie2000_val});
end




