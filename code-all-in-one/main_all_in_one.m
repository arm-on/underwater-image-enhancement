img_path = input('Please type the path of an image: ','s');
plot_output = 1;
num_levels = input('Please specify the number of pyramid levels: ');
img = imread(img_path);
%% red (and blue) channel compensation
compensated_img = compensate_channel(img, 1, 'red');
compensated_img = im2uint8(compensated_img); % reverse of im2double
compensated_img = compensate_channel(compensated_img, 1, 'blue');
compensated_img = im2uint8(compensated_img); % reverse of im2double

%% white balancing (using the gray world algorithm)
white_balanced_img = apply_gray_world(compensated_img, 20);

%% gamma correction (the first input of the multiscale fusion)
gamma_value = 1.2;
gamma_corrector = vision.GammaCorrector(gamma_value, 'Correction', 'gamma');
gamma_corrected_img = gamma_corrector(white_balanced_img);

%% sharpening (the second input of the multiscale fusion)
sharpened_img = sharpen(gamma_corrected_img);

%% Computing the Weight Map

% laplacian contrast weight
lc_weights_gc = laplacian_constrast_weights(gamma_corrected_img);
lc_weights_sh = laplacian_constrast_weights(sharpened_img);

% saliency weights computation
saliency_weights_gc = compute_saliency_weights(gamma_corrected_img);
saliency_weights_sh = compute_saliency_weights(sharpened_img);

% measuring the saturation weights
saturation_weights_gc = compute_saturation_weights(gamma_corrected_img);
saturation_weights_sh = compute_saturation_weights(sharpened_img);

% calculating the aggregated weight maps
%_Gamma Corrected
aggregated_weights_gc = lc_weights_gc + saliency_weights_gc + saturation_weights_gc;
%_Sharpened
aggregated_weights_sh = lc_weights_sh + saliency_weights_sh + saturation_weights_sh;

% Normalizing the aggregated weights
regularization = 0.1;
[normalized_weights_gc, normalized_weights_sh] = normalize_weights(aggregated_weights_gc, aggregated_weights_sh, regularization);

%% Pyramid Generation

% decomposing each source input is decomposed into a Laplacian pyramid
gaussian_pyramid_gc = generate_gaussian_pyramid(gamma_corrected_img, num_levels);
laplacian_pyramid_gc = generate_laplacian_pyramid(gaussian_pyramid_gc); % *

gaussian_pyramid_sh = generate_gaussian_pyramid(sharpened_img, num_levels);
laplacian_pyramid_sh = generate_laplacian_pyramid(gaussian_pyramid_sh); % *

% decomposing the normalized weight maps into Gaussian pyramids
gaussian_pyramid_gc_weights = generate_gaussian_pyramid(normalized_weights_gc, num_levels); % *
gaussian_pyramid_sh_weights = generate_gaussian_pyramid(normalized_weights_sh, num_levels); % *


% now, we need 'laplacian_pyramid_gc', 'laplacian_pyramid_sh',
% 'gaussian_pyramid_gc_weights', 'gaussian_pyramid_sh_weights'
pyramid = multiscale_fusion(laplacian_pyramid_gc, laplacian_pyramid_sh, gaussian_pyramid_gc_weights, gaussian_pyramid_sh_weights);

fused_pyramid = uint8(pyramid_fusion(pyramid));

images = {img, compensated_img, white_balanced_img, gamma_corrected_img, ...
    sharpened_img, fused_pyramid};

titles = {'Original Image', 'Channel Compensated', 'White Balanced', 'Gamma Corrected', ...
    'Sharpened', 'Fused Pyramid'};

if plot_output == 1
    plot_results(images, titles, [2,3]);
end

%% Evaluation
fused_eval = evaluate(fused_pyramid, img);
wb_eval = evaluate(white_balanced_img, img);
rc_eval = evaluate(compensated_img, img);
gc_eval = evaluate(gamma_corrected_img, img);
sh_eval = evaluate(sharpened_img, img);

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

function compensated_img = compensate_channel(img, alpha, channel_name)
%{
    Description: The function receives an input image along with a
    parameter (named alpha) and the name of the channel to perform the
    compensation. Then, it compensates for the given channel according to
    the formula written below.

    red(x,y) = red(x,y) + alpha*(g_mean-r_mean)*(1-red(x,y))*green(x,y)

    blue(x,y) = blue(x,y) + alpha*(g_mean-b_mean)*(1-blue(x,y))*green(x,y)

    Input:
        - img: an RGB image (a mxnx3 tensor)
        - alpha: a parameter of the formula
        - channel_name: the channel on which the compensation should be
        performed

    Output:
        - compensated_img: the result of performing the compensation
%}
    
    % re-scaling the image so that each of its pixels have an intensity
    % value in [0,1].
    scaled_img = im2double(img);

    % decomposing the image to its three channels
    red_channel = scaled_img(:,:,1);
    green_channel = scaled_img(:,:,2);
    blue_channel = scaled_img(:,:,3);
    
    % calculating the mean of the required channels by converting each
    % required channel to a vector and then calculating the mean
    green_channel_mean = mean(green_channel(:));
    if strcmp(channel_name, 'red')
        second_channel_mean = mean(red_channel(:));
    elseif strcmp(channel_name, 'blue')
        second_channel_mean = mean(blue_channel(:));
    end
    
    % performing the compensation
    compensated_channel = repmat(green_channel, 1);
    for h=1:height(compensated_channel)
        for w=1:width(compensated_channel)
            if strcmp(channel_name, 'red')
                compensation_channel = red_channel;
            elseif strcmp(channel_name, 'blue')
                compensation_channel = blue_channel;
            end
            curr_pixel_intensity = compensation_channel(h, w);
            compensated_channel(h, w) = curr_pixel_intensity+ alpha*(green_channel_mean-second_channel_mean)*(1-curr_pixel_intensity)*green_channel(h,w);
        end
    end
    
    % combining the channels
    compensated_img = repmat(scaled_img, 1);
    if strcmp(channel_name, 'red')
        compensated_img(:,:,1) = compensated_channel;
    elseif strcmp(channel_name, 'blue')
        compensated_img(:,:,3) = compensated_channel;
    end
    
end

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

function magnify_on_click(handler, img_title)
%{
    Description: The function catches a handler object and extracts the
    data from it. Then, it displays a figure containing that piece of data.

    Input:
        - handler: a handler object
        - img_title: the title of the figure

    Output: None
%}
    data = get(handler,'CData');    
    figure('Name',img_title);
    imshow(data);
    title(img_title);
end

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

function [normalized_weights_gc, normalized_weights_sh] = normalize_weights(agg_gc, agg_sh, regularization)
%{
    Description: The function first takes the aggregated weights for each
    pixel position of two pictures (gamma corrected and sharpened version
    of an original image), and normalizes the weights with respect to a
    regularization term that ensures that each input contributes to the
    output.

    Input:
        - agg_gc: the aggregated weight map of the gamma corrected image
        - agg_sh: the aggregated weight map of the sharpened image
        - regularization: the value of the regularization term
    Output:
        - normalized_weights_gc: the normalized weights with respect to the
        regularization term (for the gamma corrected image)
        - normalized_weights_sh: the normalized weights with respect to the
        regularization term (for the sharpened image)
%}
    normalized_weights_gc = zeros(height(agg_gc), width(agg_gc));
    normalized_weights_sh = zeros(height(agg_gc), width(agg_gc));
    
    for h=1:height(normalized_weights_gc)
        for w=1:width(normalized_weights_gc)
            normalized_weights_gc(h,w) = ...
                (agg_gc(h,w)+regularization)/(agg_gc(h,w)+agg_sh(h,w)+2*regularization);
            normalized_weights_sh(h,w) = ...
                (agg_sh(h,w)+regularization)/(agg_gc(h,w)+agg_sh(h,w)+2*regularization);
        end
    end
    
end

function [mean_fused, mean_rc, mean_wb, mean_gc, mean_sh] = overall_evaluation(folder_path, img_format, max_pics)
%{
    Description: The function receives a folder's path in which there are a
    lot of pictures in a specified format, and calculates the mean of each
    metric on those images.

    Input:
        - folder_path: the path of the folder containing the images
        - img_format: the format of the pictures (e.g., png, jpg)
        - max_pics: maximum number of pictures to consider in the
        evaluation
    
    Output:
        - mean_fused: mean NIQE, UCIQE, and CIE2000 evaluation metrics for
        the fused pyramids
        - mean_rc: mean NIQE, UCIQE, and CIE2000 evaluation metrics for
        the channel compensated image
        - mean_wb: mean NIQE, UCIQE, and CIE2000 evaluation metrics for
        the white-balanced image
        - mean_sh: mean NIQE, UCIQE, and CIE2000 evaluation metrics for
        the sharpened image
        - mean_gc: mean NIQE, UCIQE, and CIE2000 evaluation metrics for
        the gamma corrected_image
%}
    files = dir(fullfile(folder_path, strcat('*.',img_format)));
    files = files(1:max_pics);
    
    file_paths = cell(length(files), 1);
    for i=1:length(files)
        file_paths{i} = strcat(folder_path,'/',files(i).name);
    end
    
    fused_all = cell(length(files), 1);
    rc_all = cell(length(files), 1);
    wb_all = cell(length(files), 1);
    gc_all = cell(length(files), 1);
    sh_all = cell(length(files), 1);
    
    for i=1:length(file_paths)
        disp(i);
        [fused_eval, rc_eval, wb_eval, gc_eval, sh_eval] = main(file_paths{i}, 0);
        fused_all{i} = fused_eval;
        rc_all{i} = rc_eval;
        wb_all{i} = wb_eval;
        gc_all{i} = gc_eval;
        sh_all{i} = sh_eval;
    end
    
    mean_fused = sum(cell2mat(fused_all))/length(fused_all);
    mean_gc = sum(cell2mat(gc_all))/length(gc_all);
    mean_wb = sum(cell2mat(wb_all))/length(wb_all);
    mean_sh = sum(cell2mat(sh_all))/length(sh_all);
    mean_rc = sum(cell2mat(rc_all))/length(rc_all);
    
end

function plot_results(images, titles, plot_rows_cols)
%{
    Description: The function takes a cell array of images along with
    another array corresponding to their titles as input, and plots the
    images in a way that they can be magnified by clicking on them.
    
    Input:
        - images: a cell array of images
        - titles: another cell array of titles corresponding to the images
        - plot_rows_cols: a 1D vector containing the number of rows and
        columns of the plot

    Output: None
%}
    num_images = size(images, 2);
    
    plot_rows = plot_rows_cols(1);
    plot_cols = plot_rows_cols(2);
    
    figure('Name', 'Results');
    axes = cell(num_images,1);
    for i=1:num_images
        subplot(plot_rows, plot_cols, i);
        axes{i} = imshow(images{i});
        title(titles{i});
    end
    
    for i=1:num_images
        set(axes{i},'ButtonDownFcn',@(h,e)magnify_on_click(h, titles{i}));
    end
end


function fused_pyramid = pyramid_fusion(pyramid)
%{
    Description: The function computes the summation of the fused
    contribution of all levels of a given pyramid (after appropriate
    upsampling) and produces a single output image.

    Input:
        - pyramid: an image pyramid (a cell array)

    Output:
        - fused_pyramid: the result of summing all of the levels within the
        pyramid considering the appropriate upsampling whenever needed
%}
    % measuring the number of levels existing in the pyramid
    num_levels = size(pyramid, 1);
    
    % doing the summation
    for i=num_levels:-1:2
        curr_pyramid = impyramid(pyramid{i}, 'expand');
        curr_pyramid = imresize(curr_pyramid, [height(pyramid{i-1}), width(pyramid{i-1})]);
        pyramid{i-1} = pyramid{i-1} + curr_pyramid;
    end
    
    fused_pyramid = pyramid{1};
end


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
