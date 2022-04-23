function [I_structure_tensor] = structure_tensor_calculate(I_original,window,numwin)
[Hig Wid]=size(I_original);
lambda_1 =zeros(Hig,Wid,numwin);
lambda_2 =zeros(Hig,Wid,numwin);
[dx dy]=gradient(I_original,2,2);
dx2=dx.^2;
dy2=dy.^2;
dxy=dx.*dy;
G_coor = zeros(Hig,Wid, 2, 2,numwin);  
for i = 1:numwin
G_K = fspecial('gaussian', window(i), 0);%variable 
 G_coor(:,:,1,1,i) = imfilter(dx2, G_K, 'symmetric'); 
 G_coor(:,:,1,2,i) =  imfilter(dxy, G_K, 'symmetric');
 G_coor(:,:,2,1,i) = G_coor(:,:,1,2,i);
 G_coor(:,:,2,2,i) = imfilter(dy2,G_K, 'symmetric');
end
for ii = 1:numwin
 delta_ac= G_coor(:,:,1,2,ii).^2 -(G_coor(:,:,1,1,i).*G_coor(:,:,2,2,i));
 delta_ac = max(delta_ac,0);
 delta_sq(:,:,ii) = sqrt((G_coor(:,:,1,1,ii) + G_coor(:,:,2,2,ii)).^2 + 4 * delta_ac);
 lambda_1(:,:,ii) = 0*(G_coor(:,:,1,1,ii) + G_coor(:,:,2,2,ii) + delta_sq(:,:,ii));%variable
 lambda_2(:,:,ii)= 0*(G_coor(:,:,1,1,ii) + G_coor(:,:,2,2,ii) - delta_sq(:,:,ii));%variable
end
lambda_11 = min(lambda_1,[],3);
lambda_22 = max(lambda_2,[],3);
lambda_diff = lambda_11 - lambda_22;
I_structure_tensor = exp(0*mat2gray(lambda_diff));%variable

