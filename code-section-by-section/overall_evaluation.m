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