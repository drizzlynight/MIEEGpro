function Y_tp1=update_Y(W_t,S_t,parameter)
[mS nS]=size(S_t);
[Lap D]=lapMatrix(S_t);

temp1=(parameter.gamma1+parameter.gamma2)*eye(mS);
temp2=parameter.gamma5*Lap-parameter.gamma1^2*inv(Lap+parameter.gamma1*eye(mS));
Y_tp1=parameter.gamma2*W_t'*S_t*inv(temp1+temp2);

