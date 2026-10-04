function S = sinkhorn_proj(S, max_iter, tol)
% Sinkhorn-Knopp 投影到双随机矩阵集合
if nargin < 2, max_iter = 50; end
if nargin < 3, tol = 1e-6; end
n = size(S, 1);
S(S < 0) = 0;
S = S + eps;
for iter = 1:max_iter
    r = sum(S, 2);
    S = S ./ r;
    c = sum(S, 1);
    S = S ./ c;
    if norm(sum(S,2)-1,1) < tol && norm(sum(S,1)-1,1) < tol
        break;
    end
end
S(S < 0) = 0;
S = S ./ sum(S,2);
S = S ./ sum(S,1);
end