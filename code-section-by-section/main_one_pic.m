img_path = input('Please type the path of an image: ','s');
plot_output = 1;
num_levels = input('Please specify the number of pyramids: ');
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