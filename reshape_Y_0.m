function Y_0_matrix=reshape_Y_0(Y_0,c)
index=find(Y_0==1);
Y_0(index)=1;
index1=find(Y_0==2);
Y_0(index1)=2;

mL_0=length(Y_0);
Y_0_matrix=zeros(c,mL_0);

for i=1:mL_0
    Y_0_matrix(Y_0(i),i)=1;
end