close all;
clear all;
clc;

load EEGdataset.mat;
data=EEGdataset;
%data(isnan(data))=0;
L_true=data(:,1);
[mD nD]=size(data);

parameter.c=4; % number of clusters
parameter.nlabel=floor(0.0 *mD); % random setting for number of labeled EEG signals
parameter.iter=50; % the number of internal iterations
parameter.eta=0.01; % learning rate
parameter.epsilon=0.01; % internal convergence parameter

L_true=L_true(parameter.nlabel+1:end,:)';  % true labels without labled EEG    
l_true_matrix=l_true_reshape(L_true,parameter.c);


Result=[];

for parameter_sigma=10.^[-2 -1 0 1 2]
   parameter.sigma=parameter_sigma;
    for parameter_gamma1=10.^[-2 -1 0 1 2]
        parameter.gamma1=parameter_gamma1;
        for parameter_gamma2=10.^[-2 -1 0 1 2]
            parameter.gamma2=parameter_gamma2; 
            for parameter_gamma3=10.^[-2 -1 0 1 2]
                parameter.gamma3=parameter_gamma3;
                for parameter_gamma4=10.^[-2 -1 0 1 2]
                    parameter.gamma4=parameter_gamma4;
                    for parameter_gamma5=10.^[-2 -1 0 1 2]
                        parameter.gamma5=parameter_gamma5;
                        tic;
                        [P_opt, Y_opt, W_opt, S_opt, P_0, Y_0, W_0, S_0, F_t, F_tp1, wh_time]=MIEEGpro(data,parameter);
                        time=toc; 
                        count=0;
                        Y_label=[];
                        Y_opt=Y_opt(:,parameter.nlabel+1:end); 
                        [r cc]=size(Y_opt);
                        for ct=1:cc
                            for rt=1:r
                                if Y_opt(rt,ct)==1
                                    Y_label(ct)=rt;
                                    count=count+1;
                                end
                                if count==0
                                    Y_label(ct)=randperm(parameter.c,1);
                                end
                            end
                        end
                        ARI=rand_index(Y_label,L_true,'adjusted');
                        RI=rand_index(Y_label,L_true);
                        NMI=nmi(Y_label,L_true);
                        [micro, macro]=micro_macro_PR(Y_label,L_true);
                        kap=kappaindex(Y_label,L_true,parameter.c);
                        runtime=time;
                        Result=[Result; RI ARI NMI macro.fscore micro.fscore kap runtime];
                    end
                end
            end
        end
    end
end
sortrows(Result)
