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