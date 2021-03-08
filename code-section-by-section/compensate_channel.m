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