function [Lap, D] = lapMatrix(M)
% 标准拉普拉斯矩阵 L = D - M
D = diag(sum(M, 2));
Lap = D - M;
end


% % function [Lap D]=lapMatrix(M) 
% % D=diag(sum(M,2)-diag(M)); % calculate the degree matrix of M
% % Lap=D-M; % calculate Laplacian matrix of similarity matrix
% 
% function [Lap D]=lapMatrix(M) 
% % 标准拉普拉斯矩阵 L = D - M, D = diag(sum(M,2))
% D = diag(sum(M,2));
% Lap = D - M;
% end