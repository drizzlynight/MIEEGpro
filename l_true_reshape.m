function true_L_matrix=l_true_reshape(L_true,C)
index=find(L_true==1);
L_true(index)=1;
index1=find(L_true==2);
L_true(index1)=2;
mY=length(L_true);
true_L_matrix=zeros(C,mY);
for j=1:mY
    true_L_matrix(L_true(j),j)=1;
end
end
