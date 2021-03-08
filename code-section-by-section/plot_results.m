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

