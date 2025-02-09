function [Lap D]=lapMatrix(M) 
D=diag(sum(M,2)-diag(M)); % calculate the degree matrix of M
Lap=D-M; % calculate Laplacian matrix of similarity matrix