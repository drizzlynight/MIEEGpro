function [L_S D_S S]=GauSimilarity(E,sigma)
[mE nE]=size(E);
D_S=zeros(mE,mE);
for i=1:mE
    for j=i:mE
%         g=max(xcorr(E(i,:)-E(j,:)));
        g=norm(E(i,:)-E(j,:));
        S(i,j)=exp(- g^2/(2*sigma^2)); %Gaussian similarity between the i-th time series and the j-th time series
        S(j,i)=S(i,j);
    end
    D_S(i,i)=sum(S(i,:));
end
L_S=D_S-S;  %Laplacian matrix of S