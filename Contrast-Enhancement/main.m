clear all;
close all;
clc;
% Please change the path according to your save path
path = '..../...../';
cd (path);
%% Confirm accessing the sequence of images and output path of results
access_path  = ('test_image/');
im_index = dir([access_path '*.png']);
im_num = length(im_index);
% The DoGF filter with multi-scale windows
window = [1 3 5 7];
numwin = length(window(:));
DoGauss = cell(1,numwin);
for k=1:im_num
    % input image
    I_original = imread([access_path,im_index(k).name]);
    if size(I_original, 3) > 1
    I_original = rgb2gray(I_original);
    end
    I_original = double(I_original);
    [row col] = size(I_original);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Phase 1  calculate the structure_tensor!
     [I_structure_tensor] = structure_tensor_calculate(I_original,window,numwin);
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %% Phase 2 Small target enhancement and candidate areas segmentation !
       sigma1 = 0;%variable
       sigma2 = 0;%variable
     for num = 1:numwin
      G1 = fspecial('gaussian', window(num), sigma1);
      G2 = fspecial('gaussian', window(num), sigma2);
      DoGauss{num} = G1-G2;
     end
    Io_gauss = zeros(row,col,numwin);
      for imfG = 1:numwin
       Dog_result = imfilter(I_original,DoGauss{imfG},'replicate'); 
       Io_gauss(:,:,imfG) = Dog_result;
      end
     Dogf_mean = mean(Io_gauss,3);
     I_CE=  Dogf_mean ./(sqrt(I_structure_tensor));
     img_s = std(I_CE(:));
     img_m = mean(I_CE(:)>0);
     detal = 0 * img_s ./img_m ;%variable
     connected_result = I_CE > detal;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Phase 3 Small target detection based on Gaussian gradient contrast! 
     [target_result] = target_deal_score(connected_result,I_CE,I_original);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   figure;  
   subplot(2,3,1);imshow(I_original,[0 255]); title('INPUT-image');
   subplot(2,3,2);imshow(I_structure_tensor,[]); title('structure_tensor');
   subplot(2,3,3); imshow(Dogf_mean,[]); title('Dogf-mean-image');
   subplot(2,3,4); imshow(I_CE,[]); title('mean-Dog-struct-image');
   subplot(2,3,5); imshow(connected_result,[]); title('connected-resultn-image');
   subplot(2,3,6); imshow(target_result,[]); title('finall-image');
end
