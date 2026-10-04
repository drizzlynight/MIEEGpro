function W_tp1=update_W(S_t,Y_t,parameter) % calculate W
[mS nS]=size(S_t);

W_tp1=parameter.gamma2*inv(parameter.gamma2*S_t*S_t'+parameter.gamma3*eye(mS))*S_t*Y_t';