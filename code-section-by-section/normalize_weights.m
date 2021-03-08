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