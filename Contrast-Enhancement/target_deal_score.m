function [target_result] = target_deal_score(connected_result,I_CE,I_original)
  gs_k= [3, 6, 3;
         6, 12,6 ;
         3, 6, 3];
 I_original_o = double(I_original);
 I_original_g = imfilter(I_original_o,gs_k, 'symmetric'); 
 [row col] = size(I_CE);
 %% connected region
 connected_struct = bwconncomp(connected_result); 
 connected_num =  connected_struct.NumObjects;
 connected_result_final =I_CE;
 for jj = 1 : connected_num
     connected_tar = connected_struct.PixelIdxList(jj);
     connected_size = [connected_tar{1}]; 
     [row1 col1] = size(connected_size);
     clear connected_position ; 
      for j = 1: row1
         connected_position(j,1) =  mod(connected_size(j,1),row);
          if connected_position(j,1)==0
             connected_position(j,1) = row;
          end
          if connected_position(j,1)==row
             connected_position(j,2) = connected_size(j,1)/row;
          else
              connected_position(j,2) = (floor(connected_size(j,1)/row) )+1;
          end
      end   
      if  row1 > 1
      min_rowcol = min(connected_position);
      max_rowcol = max(connected_position);
      min_row = min_rowcol(1);
      max_row = max_rowcol(1);
      min_col = min_rowcol(2);
      max_col = max_rowcol(2);
      range_row = max_row - min_row +1;
      range_col = max_col - min_col +1 ;
      path_struct = I_original_g(min_row:max_row ,min_col:max_col);
      hig_point = max(path_struct(:));
      [x y]=find(path_struct == hig_point,1);
      highest_row = min_row + x -1;
      highest_col = min_col + y -1;
      else
      highest_row = connected_position(1,1);
      highest_col = connected_position(1,2);
      range_row =1;
      range_col =1;
      max_row =highest_row;
      min_row =highest_row;
      max_col =highest_col;
      min_col =highest_col;
      end
      %% Gaussian template and edge processing!
      Gauss_mode = zeros(3);
      if (highest_row == 1 && highest_col ~= 1 && highest_col ~= col)
          Gauss_mode(1,:)=1;
          Gauss_mode(2:end,:)= I_original_g((highest_row+1):(highest_row+2), (highest_col- 1):(highest_col+1));
      elseif (highest_col == 1 && highest_row ~= 1 && highest_row ~= row)
          Gauss_mode(:,1)= 1;
          Gauss_mode(:,2:end) = I_original_g((highest_row-1):(highest_row+1), (highest_col+1):(highest_col+2));
      elseif (highest_col == 1&&highest_row == 1)
          Gauss_mode(:,1)= 1;
          Gauss_mode(1,:)=1;
          Gauss_mode(2:end,2:end) = I_original_g((highest_row+1):(highest_row+2), (highest_col+1):(highest_col+2));
      elseif (highest_row == row && highest_col ~= 1 && highest_col ~= col)
          Gauss_mode(3,:)=1;
          Gauss_mode(1:2,:) = I_original_g((highest_row-1):(highest_row), (highest_col- 1):(highest_col+1));
      elseif (highest_col == col&& highest_row ~= 1 && highest_row ~= row) 
          Gauss_mode(:,3)=1;
          Gauss_mode(:,1:2) = I_original_g((highest_row-1):(highest_row+1), (highest_col- 1):(highest_col));
      elseif (highest_col == col && highest_row == row)
          Gauss_mode(:,3)= 1;
          Gauss_mode(3,:)=1;
          Gauss_mode(1:2,1:2) = I_original_g((highest_row-1):(highest_row), (highest_col- 1):(highest_col));
       elseif (highest_col == 1 && highest_row == row)
          Gauss_mode(:,1)= 1;
          Gauss_mode(3,:)=1;
          Gauss_mode(1:2,2:end) = I_original_g((highest_row-1):(highest_row), (highest_col+1):(highest_col+2));
       elseif (highest_col == col && highest_row == 1)
          Gauss_mode(1,:)=1;
          Gauss_mode(:,3)= 1;
          Gauss_mode(2:end,1:2)= I_original_g((highest_row+1):(highest_row+2), (highest_col- 1):(highest_col));
       else 
          Gauss_mode = I_original_g((highest_row-1):(highest_row+1), (highest_col- 1):(highest_col+1));
      end 
    %% Gaussian gradient contrast!
       dx11 = Gauss_mode(1,2) - Gauss_mode(1,1);
       dx13 = Gauss_mode(1,2) - Gauss_mode(1,3);
       dx31 = Gauss_mode(3,2) - Gauss_mode(3,1);
       dx33 = Gauss_mode(3,2) - Gauss_mode(3,3);
       dy11 = Gauss_mode(2,1) - Gauss_mode(1,1);
       dy13 = Gauss_mode(2,3) - Gauss_mode(1,3);
       dy31 = Gauss_mode(2,1) - Gauss_mode(3,1);
       dy33 = Gauss_mode(2,3) - Gauss_mode(3,3);
    score_gauss = (dx11 <0) + (dx13 <0) + (dx31 <0) + (dx33<0) + (dy11 <0) + (dy13 <0) + (dy31 <0) + (dy33<0);
      if ( (range_row == 1 &&score_gauss>0)|| (range_col == 1&&score_gauss>0))
          connected_result(min_row:max_row ,min_col:max_col)= 0;
      elseif  ((range_row /range_col)>3|| (range_col /range_row)>3)
          connected_result(min_row:max_row ,min_col:max_col)= 0;
      elseif   score_gauss > 1
          connected_result(min_row:max_row ,min_col:max_col)= 0;
      end
      connected_result_final = connected_result_final.* connected_result;
 end
   target_result = (connected_result_final.^2);
 