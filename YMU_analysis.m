clear;clf;close;clc

fprintf('loading data...\n')
load('depression_aal2TC_YMU.mat')
load('covariates_YMU.mat')
load('UKB\UKB_bscmodel_aal2.mat')
%% processing data
n_ROI = 94;
n_person = length(RSts_aal2_YMU);
% calculating FC
for person = 1:n_person
    X = corr(RSts_aal2_YMU{person}(:,1:n_ROI));
    xtemp = [];
    for edge = 1:n_ROI
        xtemp = [xtemp,X1(edge,edge+1:n_ROI)];
    end
    FC(person,:) = xtemp;
end
FC = 0.5*log((1+FC)./(1-FC)); % Fisher-z transformation
fprintf('finish preparation.\n')

regresscov = [cov_YMU.age,cov_YMU.age.^2,cov_YMU.age.^3]; % age to its third order as covariates to be regressed out from FC.
fprintf('starting regressing out covariates\n')

for i = 1:(n_ROI*(n_ROI-1)/2)
    glmstruct = fitglm(regresscov,FC(:,i));
    b = glmstruct.Coefficients.Estimate;
    FC_regressed(:,i) = FC(:,i)-[ones(n_person,1),regresscov]*b;
end

% predicting on test set

[~,prob] = predict(svmstruct,FC_regressed); % calculating the classification accuracy    
fprintf('finish calculating\n');

BSC_aal2 = normcdf(prob(:,2)); % calculating the probability (i.e. the brain sex continuum)
output = double(BSC_aal2>0.5);
acc_tot = 1-sum(abs(output-cov_YMU.sex))/n_person;

acc_hc_f = 1-sum(cov_YMU.sex(cov_YMU.diag==0 & cov_YMU.sex==0) ~= double(BSC_aal2(cov_YMU.diag==0 & cov_YMU.sex==0)>0.5))/sum(cov_YMU.diag==0 & cov_YMU.sex==0);
acc_hc_m = 1-sum(cov_YMU.sex(cov_YMU.diag==0 & cov_YMU.sex==1) ~= double(BSC_aal2(cov_YMU.diag==0 & cov_YMU.sex==1)>0.5))/sum(cov_YMU.diag==0 & cov_YMU.sex==1);
acc_MDD_f = 1-sum(cov_YMU.sex(cov_YMU.diag==1 & cov_YMU.sex==0) ~= double(BSC_aal2(cov_YMU.diag==1 & cov_YMU.sex==0)>0.5))/sum(cov_YMU.diag==1 & cov_YMU.sex==0);
acc_MDD_m = 1-sum(cov_YMU.sex(cov_YMU.diag==1 & cov_YMU.sex==1) ~= double(BSC_aal2(cov_YMU.diag==1 & cov_YMU.sex==1)>0.5))/sum(cov_YMU.diag==1 & cov_YMU.sex==1);
%% bootstrap test to measure difference between accuracy
% bootstrap groups
BSC_f_hc = BSC_aal2(cov_YMU.sex==0 & cov_YMU.diag==0);
BSC_f_mdd = BSC_aal2(cov_YMU.sex==0 & cov_YMU.diag==1);
BSC_m_hc = BSC_aal2(cov_YMU.sex==1 & cov_YMU.diag==0);
BSC_m_mdd = BSC_aal2(cov_YMU.sex==1 & cov_YMU.diag==1);
n_bstp = 10000;
clear acc;clear acc_a
for i_bstp = 1:n_bstp
    acc(1,i_bstp) = 1-sum(double(BSC_f_hc(randsample(length(BSC_f_hc),length(BSC_f_hc),1))>0.5)~=0)/length(BSC_f_hc);
    acc(2,i_bstp) = 1-sum(double(BSC_f_mdd(randsample(length(BSC_f_mdd),length(BSC_f_mdd),1))>0.5)~=0)/length(BSC_f_mdd);
    acc(3,i_bstp) = 1-sum(double(BSC_m_hc(randsample(length(BSC_m_hc),length(BSC_m_hc),1))>0.5)~=1)/length(BSC_m_hc);
    acc(4,i_bstp) = 1-sum(double(BSC_m_mdd(randsample(length(BSC_m_mdd),length(BSC_m_mdd),1))>0.5)~=1)/length(BSC_m_mdd);
    %fprintf('boostrapping %.0f \n',i_bstp)
end
% calculating the difference of accuracy in females and males
acc_diff_f = sort(acc(2,:)-acc(1,:),'ascend');
acc_diff_m = sort(acc(4,:)-acc(3,:),'ascend');
fprintf('female: %.4f to %.4f; male: %.4f to %.4f \n',acc_diff_f(n_bstp*0.025+1),acc_diff_f(n_bstp*0.975),acc_diff_m(n_bstp*0.025+1),acc_diff_m(n_bstp*0.975));

%% Figure 2a and b
ind_plot_HC = find(cov_YMU.sex==1&cov_YMU.diag==0);
ind_plot_MDD = find(cov_YMU.sex==1&cov_YMU.diag==1);
boxplot([BSC_aal2(ind_plot_HC);BSC_aal2(ind_plot_MDD)],[ones(length(ind_plot_HC),1);2*ones(length(ind_plot_MDD),1)]);
hold on
plot([ones(length(ind_plot_HC),1)+randn(length(ind_plot_HC),1)*0.1;2*ones(length(ind_plot_MDD),1)+randn(length(ind_plot_MDD),1)*0.1],[BSC_aal2(ind_plot_HC);BSC_aal2(ind_plot_MDD)],'+')
set(gca,'xticklabel',{'HC','MDD'})

mdl = fitglm(cov_YMU{cov_YMU.sex==0,[4,2,3,7]},BSC_aal2(cov_YMU.sex==0)); % 4=diag, 2=age, 3=sex, 7=meanFD
mdl = fitglm(cov_YMU{cov_YMU.sex==1,[4,2,3,7]},BSC_aal2(cov_YMU.sex==1));
%% HAMD plot only YMU (figure 2f)
mdl = fitglm(cov_YMU{:,[2,3,7]},BSC_aal2);
x = BSC_aal2 - [ones(length(BSC_aal2),1),cov_YMU{:,[2,3,7]}]*mdl.Coefficients.Estimate;

mdl = fitglm(cov_YMU{:,[2,3,7]},cov_YMU.HAMD);
y = cov_YMU.HAMD - [ones(length(BSC_aal2),1),cov_YMU{:,[2,3,7]}]*mdl.Coefficients.Estimate;
plot(x,y)
lsline

