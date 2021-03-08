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
