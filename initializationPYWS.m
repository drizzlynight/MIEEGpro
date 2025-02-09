
function [P_0, Y_0, W_0, S_0, M]=initializationPYWS(data,parameter)

data_label=data(1:parameter.nlabel,:);
data_unlabel=data(parameter.nlabel+1:end,2:end);

L_label=data_label(1:parameter.nlabel,1);

a=0.99;
b=0.01;
[mD nD]=size(data(:,2:end));

[L_M D_M M]=GauSimilarity(data(:,2:end),parameter.sigma); % calculate Gaussian similarity of EEG signals


%---------------------------------initialize S0----------------------------------------------
numL=parameter.nlabel;

S_0=1/mD*ones(mD);
for i=1:numL
    for j=1:numL
        if L_label(i)==L_label(j)
            S_0(i,j)=1;
            S_0(j,i)=1;
        else if L_label(i)~=L_label(j)
            S_0(i,j)=0;
            S_0(j,i)=0;
            end
        end
    end
end

% S_0=M;
for i=1:numL
    S_0(i,:)=S_0(i,:)/sum(S_0(i,:)); % S*1^T=1;
end

%---------------------------------initialize Y0----------------------------------------------
[L_unlabel Centroid]=kmeans(data_unlabel,parameter.c); 
L_unlabel(find(isnan(L_unlabel)))=randperm(parameter.c,1);
L_ini=[L_label;L_unlabel];
Y_0=reshape_Y_0(L_ini,parameter.c);

%---------------------------------initialize W0----------------------------------------------
W_0=(Y_0*inv(M))';

%---------------------------------initialize P0----------------------------------------
P_0=b*Y_0*inv(eye(mD)-a*S_0);
