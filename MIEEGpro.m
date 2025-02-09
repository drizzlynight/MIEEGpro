function [P_opt, Y_opt, W_opt, S_opt, P_0, Y_0, W_0, S_0, F_t,F_tp1,wh_time]=MIEEGpro(data,parameter)

L_true=data(:,1);
unregu_data=data(:,2:end);

data=z_regularization(unregu_data);% regularize EEG signals;

data=[L_true data];

[P_0, Y_0, W_0, S_0, M]=initializationPYWS(data,parameter);

P_tp1=P_0;
Y_tp1=Y_0;
W_tp1=W_0;
S_tp1=S_0;
gap=100;
F_tp1=10^15;
F_t=F_tp1+10^10;
L_label=data(1:parameter.nlabel,1);
[L_S_tp1 D_S_tp1]=lapMatrix(S_tp1);
wh_time=0;
count=0; 
while gap>parameter.epsilon 
   count=count+1;
    F_tp1=obj_function(M,P_tp1,Y_tp1,W_tp1,S_tp1,L_S_tp1,parameter); % calculate the value of objective function;
    gap=F_t-F_tp1;
    if isnan(F_tp1)
        break;
    end

    P_tp1=update_P(Y_tp1,S_tp1,parameter); %update R_tp1;
    W_tp1=update_W(S_tp1,Y_tp1,parameter); %update W_tp1;
    [D_S_tp1 S_tp1]=update_S(M,S_0,W_tp1,P_tp1,Y_tp1,L_label,parameter); % update S_tp1;
    Y_tp1=update_Y(W_tp1,S_tp1,parameter); %update L_tp1;
    
    F_t=F_tp1;
    wh_time=wh_time+1;
    if wh_time==50
        break;
    end
end

P_opt=P_tp1; % optimal U_opt;
W_opt=W_tp1; % optimal W_opt;
S_opt=S_tp1; % optimal S_opt;
% Y_tp1

[mL nL]=size(Y_tp1);
for i=1:nL
    L_max=max(Y_tp1(:,i));
    for j=1:mL
        if Y_tp1(j,i)==L_max
            Y_tp1(j,i)=1;
        else
            Y_tp1(j,i)=0;
        end
    end
end

Y_opt=Y_tp1;

% [mL nL]=size(L_tp1);
% L_opt=zeros(mL,nL);
% l_max=max(L_tp1);
% 
% for j=1:mL
%     l_index=find(L_tp1(:,j)==l_max(j));
%     L_opt(l_index,j)=1; % optimal L_opt;
% end

