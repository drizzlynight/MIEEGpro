function P_tp1=update_P(Y_t,S_t,parameter)
[mS nS]=size(S_t);
[Lap D]=lapMatrix(S_t);

P_tp1=parameter.gamma1*Y_t*inv(Lap+parameter.gamma1*eye(mS));

