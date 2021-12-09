function [C_com, R_com, C_mid, R_mid] = func_find_center_v2(I)

[Nr, Nc] = size(I);

x = 1:Nc;
y = 1:Nr;

[X,Y] = meshgrid (x,y);

%% center of mass
C_com = sum(X(I))/sum(I(:));
R_com = sum(Y(I))/sum(I(:));

%% average of left and right, top and bottom
sum_C = sum(I,1); 

Ix_fish = find(sum_C>0);
C_left = min(Ix_fish);
C_right = max(Ix_fish);
C_mid = (C_left+C_right)/2; 

sum_R = sum(I,2); 
Ix_fish = find(sum_R>0);
R_btm = min(Ix_fish);
R_top = max(Ix_fish);
R_mid = (R_btm+R_top)/2;
