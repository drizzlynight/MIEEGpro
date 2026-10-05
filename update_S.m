function [D_S_tp1, S_tp1, L_S_tp1] = update_S(M, S_current, W_t, P_t, Y_t, L_label, parameter)

S_t = S_current;
[mS, nS] = size(S_t);
t = 1;
while t < parameter.iter + 1
    [Lap, D] = lapMatrix(S_t);
    %dP = parameter.gamma1 * Y_t * inv(Lap + parameter.gamma1 * eye(mS)) * (inv(Lap + parameter.gamma1 * eye(mS)))';
    dL = -eye(mS);
    %P1 = (P_t * Lap + parameter.gamma1 * (P_t - Y_t))' * dP;
    P2 = 0.5 * (P_t' * P_t + parameter.gamma5 * Y_t' * Y_t) * dL;
    P3 = parameter.gamma2 * W_t * (W_t' * S_t - Y_t);
    P4 = parameter.gamma4 * (S_t - M);
    derS = P2 + P3 + P4; %P1 + P2 + P3 + P4;
    S_t = S_t - parameter.eta * derS;
    
    % 非负截断，对角线置0，双随机投影
    S_t(S_t < 0) = 0;
    S_t(1:mS+1:end) = 0;
    S_t = sinkhorn_proj(S_t, 50, 1e-6);
    
    t = t + 1;
end

S_tp1 = S_t;
[D_S_tp1, L_S_tp1] = lapMatrix(S_tp1);
D_S_tp1 = diag(sum(S_tp1, 2));
end



% % // function [D_S_tp1 S_tp1]=update_S(M,S_0,W_t,P_t,Y_t,L_label,parameter)
% % 
% % // S_t=S_0;
% % // [mS nS]=size(S_t);
% % // [Lap D]=lapMatrix(S_t);
% % // [numL colL]=size(L_label);
% % // t=1;
% % // while t<parameter.iter+1
% % //     dP=parameter.gamma1*Y_t*inv(Lap+parameter.gamma1*eye(mS))*(inv(Lap+parameter.gamma1*eye(mS)))';
% % //     dL=zeros(mS)-eye(mS);
% % //     P1=(P_t*Lap+parameter.gamma1*(P_t-Y_t))'*dP;
% % //     P2=0.5*(P_t'*P_t+parameter.gamma5*Y_t'*Y_t)*dL;
% % //     P3=parameter.gamma2*W_t*(W_t'*S_t-Y_t);
% % //     P4=parameter.gamma4*(S_t-M);
% % //     derS=P1+P2+P3+P4; % derivative of S
% % //     S_t=S_t-parameter.eta*derS;
% % // %     for i=1:mS
% % // %         for j=1:mS
% % // %             if S_t(i,j)<=0
% % // %                 S_t(i,j)=0;
% % // %             end
% % // %         end
% % // %     end
% % // % 
% % // %     for i=1:parameter.nlabel
% % // %         for j=1:parameter.nlabel
% % // %             if L_label(i)==L_label(j)
% % // %                 S_t(i,j)=1;
% % // %                 S_t(j,i)=1;
% % // %             else if L_label(i)~=L_label(j)
% % // %                     S_t(i,j)=0;
% % // %                     S_t(j,i)=0;
% % // %                 end
% % // %             end
% % // %         end
% % // %     end
% % //     for i=1:mS
% % //         S_t(i,:)=S_t(i,:)/sum(S_t(i,:)); % S*1^T=1;
% % //     end
% % //     [L_S_t D_S_t]=lapMatrix(S_t);
% % //     t=t+1;
% % // end
% % 
% % // S_tp1=S_t;
% % // D_S_tp1=D_S_t;
% % // L_S_tp1=L_S_t;
% 
% function [D_S_tp1 S_tp1]=update_S(M,S_0,W_t,P_t,Y_t,L_label,parameter)
% 
% S_t=S_0;
% [mS nS]=size(S_t);
% [numL colL]=size(L_label);
% t=1;
% while t<parameter.iter+1
%     % 每次迭代前重新计算 Lap 和 D
%     [Lap D]=lapMatrix(S_t);
%     dP=parameter.gamma1*Y_t*inv(Lap+parameter.gamma1*eye(mS))*(inv(Lap+parameter.gamma1*eye(mS)))';
%     dL=zeros(mS)-eye(mS);
%     P1=(P_t*Lap+parameter.gamma1*(P_t-Y_t))'*dP;
%     P2=0.5*(P_t'*P_t+parameter.gamma5*Y_t'*Y_t)*dL;
%     P3=parameter.gamma2*W_t*(W_t'*S_t-Y_t);
%     P4=parameter.gamma4*(S_t-M);
%     derS=P1+P2+P3+P4;
%     S_t=S_t-parameter.eta*derS;
%     
%     % 非负截断，对角线置0
%     S_t(S_t < 0) = 0;
%     S_t(1:mS+1:end) = 0;
%     % 双随机投影
%     S_t = sinkhorn_proj(S_t, 50, 1e-6);
%     
%     % 重新计算 Lap 和 D 用于下一次迭代
%     [Lap D] = lapMatrix(S_t);
%     t=t+1;
% end
% 
% S_tp1=S_t;
% D_S_tp1=D;
% L_S_tp1=Lap;
