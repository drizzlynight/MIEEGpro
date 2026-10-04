function Obj=obj_function(M,P,Y,W,S,L_S,parameter)
part1=0.5*trace(P*L_S*P');
part2=0.5*parameter.gamma1*trace((P-Y)'*(P-Y));
part3=0.5*parameter.gamma2*trace((W'*S-Y)'*(W'*S-Y));
part4=0.5*parameter.gamma3*trace(W'*W);
part5=0.5*parameter.gamma4*trace((S-M)'*(S-M));
part6=0.5*parameter.gamma5*trace(Y*L_S*Y');
Obj=part1+part2+part3+part4+part5+part6;