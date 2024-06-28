% calc_min_sp_ppw_dist.m
% Script to calculate the minimum distance between the soft palate and the
% posterior pharyngeal wall in segmentations

% Author: Matthieu Ruthven (matthieuruthven@nhs.net)
% Last modified: 13th February 2023

close all
clear all

% Path to segmentations
%file_path = '/Users/andreia/Desktop/matthieu/segmentation/subject_6/aa_mask_7_classes';
%file_path = '/Users/andreia/Desktop/matthieu/segmentation/subject_5/ah_mask_7_classes';
%file_path = '/Users/andreia/Desktop/matthieu/segmentation/subject_10/br_mask_7_classes';
%file_path = '/Users/andreia/Desktop/matthieu/segmentation/subject_12/gc_mask_7_classes';
%file_path = '/Users/andreia/Desktop/matthieu/segmentation/subject_2/mr_mask_7_classes';
% file_path = '/Users/andreia/Desktop/matthieu/segmentation/2d_seg_results_thesis/bwce/test_aa_lr_0.0003_mbs_8_epochs_200_aug_4_bwce_eval_predict_cca.mat';
% file_path = '/Users/andreia/Desktop/matthieu/segmentation/2d_seg_results_thesis/bwce/test_ah_lr_3e-05_mbs_8_epochs_200_aug_4_bwce_eval_predict_cca.mat';
% file_path = '/Users/andreia/Desktop/matthieu/segmentation/2d_seg_results_thesis/bwce/test_br_lr_3e-05_mbs_4_epochs_200_aug_4_bwce_eval_predict_cca.mat';
% file_path = '/Users/andreia/Desktop/matthieu/segmentation/2d_seg_results_thesis/bwce/test_gc_lr_0.0003_mbs_8_epochs_200_aug_4_bwce_eval_predict_cca.mat';
% file_path = '/Users/andreia/Desktop/matthieu/segmentation/2d_seg_results_thesis/bwce/test_mr_lr_0.0003_mbs_4_epochs_200_aug_4_bwce_eval_predict_cca.mat';
file_path = '/Users/andreia/Desktop/matthieu/segmentation/2d_seg_results_thesis/wce/test_aa_lr_0.0003_mbs_4_epochs_200_aug_4_wce_eval_predict_cca.mat';
%file_path = '/Users/andreia/Desktop/matthieu/segmentation/2d_seg_results_thesis/wce/test_ah_lr_0.0003_mbs_4_epochs_200_aug_4_wce_eval_predict_cca.mat';
%file_path = '/Users/andreia/Desktop/matthieu/segmentation/2d_seg_results_thesis/wce/test_br_lr_0.0003_mbs_4_epochs_200_aug_4_wce_eval_predict_cca.mat';
%file_path = '/Users/andreia/Desktop/matthieu/segmentation/2d_seg_results_thesis/wce/test_gc_lr_0.0003_mbs_4_epochs_200_aug_4_wce_eval_predict_cca.mat';
%file_path = '/Users/andreia/Desktop/matthieu/segmentation/2d_seg_results_thesis/wce/test_mr_lr_0.003_mbs_8_epochs_200_aug_4_wce_eval_predict_cca.mat';

% Segmentation type
seg_type = 'wce';

% Load segmentations
if strcmp(seg_type, 'gt')
    
    % Number of frames
    n_frames = 105;
    
    % Preallocate array for segmentations
    segs = zeros(256, 256, n_frames);
    
    % For each frame
    for frame_num = 1:n_frames
        
        % Path to segmentation
        seg_path = sprintf('%s/10fps_aa_mask_%d.mat', file_path, frame_num);
        
        % Load segmentation
        load(seg_path)
        
        % Populate segs
        segs(:, :, frame_num) = mask_frame;
        
    end
    
else
    
    % Load segmentations
    load(file_path)
    
    % Convert segmentations to double
    segs = double(cca_predict);

end

% Value to assign pixels not is posterior head class
pix_val = size(segs, 1);

% Calculate number of frames
n_frames = size(segs, 3);

% Preallocate array for minimum distances
min_dist_array = zeros(n_frames, 1);

% For each frame
for frame_num = 1:n_frames
    
    % Extract the soft palate
    sp = (segs(:, :, frame_num) == 2);
    
    % Calculate the centroid of the soft palate
    stats = regionprops(sp);
    centroid = stats.Centroid;
    
%     % Create another colormap
%     other_cmap = [0 0 0
%                   0 0 1
%                   0 1 1
%                   0 1 0
%                   1 1 0
%                   1 0 1
%                   1 0 0];
%     % Black, blue, cyan, green, yellow, magenta, red
%      
%     % Show segmentations
%     imshow(segs(:, :, frame_num) + 1, other_cmap)
%     hold on
%     plot(round(centroid(1)), round(centroid(2)), 'r*')
%     
%     % Assemble path to file
%     file_path = sprintf('/Users/andreia/Desktop/matthieu/segmentation/2d_seg_results_thesis/10fps_aa_seg_%d_post_pp_centroid.png', frame_num);
% 
%     % Export segmentations to .png file
%     exportgraphics(gcf, file_path)

    % Remove anterior head pixels and soft palate pixels
    head = (segs(:, :, frame_num) == 1);
    head(:, 1:round(centroid(1) + 5)) = 0;
    sp(:, 1:round(centroid(1))) = 0;
    
    % Connected component analysis
    head = bwareafilt(head, 1);
    sp = bwareafilt(sp, 1);
    
%     % Show segmentations
%     imshow((head + 2 * sp) + 1, other_cmap)
%     
%     % Assemble path to file
%     file_path = sprintf('/Users/andreia/Desktop/matthieu/segmentation/2d_seg_results_thesis/10fps_aa_seg_%d_post_pp_sp_ppw.png', frame_num);
% 
%     % Export segmentations to .png file
%     exportgraphics(gcf, file_path)

    % Calculate the distance of each pixel from a soft palate pixel
    sp_dist = bwdist(sp == 1);
    
    % Assign large distances to pixels not in posterior head class
    sp_dist = sp_dist .* (head == 1) + pix_val * (head == 0);
    
    % Calculate minimum distance between soft palate and posterior
    % pharyngeal wall and update min_dist_array
    min_dist_array(frame_num) = (min(sp_dist(:)) - 1) * 1.17;
    
    % Show segmentations
    f = figure;
    imagesc(head + 2 * sp)
    axis off image
    title((min(sp_dist(:)) - 1) * 1.17)
    waitfor(f)
    
end

% Save distances as CSV file
writematrix(min_dist_array, sprintf('mr_%s.csv', seg_type))