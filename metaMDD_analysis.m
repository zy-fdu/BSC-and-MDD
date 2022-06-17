%% metaMDD analysis
clear;clc;clf;close;
fprintf('loading data...\n')
load('RSts_Power_metaMDD.mat')

% calculating FC matrix
n_person = length(RSts_Power_metaMDD);
for i_person = 1:n_person
    X = corr(RSts_Power_metaMDD{i_person});
    xtemp = [];
    for i_ROI = 1:size(RSts_Power_metaMDD{1},2)
        xtemp = [xtemp,X(i_ROI,i_ROI+1:size(RSts_Power_metaMDD{1},2))];
    end
    FC_Power(i_person,:) = xtemp;
end
FC_Power = 0.5*log((1+FC_Power)./(1-FC_Power)); % Fisher-z transformation
fprintf('finish preparation.\n')
%% calculate network-wise FC matrix
load('Power_network_order.mat')
FC_Power_new = FC_Power(power_new_index,power_new_index,:); % turning the ordering to newer version of Power's parcellation
clear FC_net
for i_net = 1:12
    for j_net = 1:12
        for i_subj = 1:size(FC_power_new,3)
            FC_net(i_net,j_net,i_subj) = nanmean(nanmean(FC_power_new(power_network(1,i_net):power_network(2,i_net),power_network(1,j_net):power_network(2,j_net),i_subj)));
            % FC power is a file indicating which network the network belongs to.
            % The first row indicates the starting ROI of the network,
            % while the second row indicates the ending ROI of the network.
        end
    end
end
%% calculating Brain sex continuum
s = unique(cov_MDD.site);
for i_site = 1:length(s)
    dummyvar_site(:,i_site) = double(cov_MDD.site==s(i_site)); % site as one additional covariant when dealing with age.
end
regresscov_metaMDD = [cov_MDD.age,cov_MDD.age.^2,cov_MDD.age.^3,dummyvar_site];
fprintf('starting regressing out covariates\n')

% regressing out age to its third order and site
for i = 1:(264*263/2)
    glmstruct = fitglm(regresscov_metaMDD,FC_Power(:,i));
    b = glmstruct.Coefficients.Estimate;
    FC_Power_regressed(:,i) = FC_Power(:,i)-[ones(n_person,1),regresscov_metaMDD]*b;
end
%

load('UKB\UKB_bscmodel_power.mat') % loading the pretrained BSC model
[~,prob] = predict(svmstruct_power,FC_Power_regressed); % calculating the classification score
clear svmstruct*
fprintf('finish calculating\n');
BSC_metaMDD = normcdf(prob(:,2)); % turning classification score to probability (i.e. the Brain sex continuum)

% calculating test accuracy in each sex*diagnosis group
for i_diag = 1:2
    for i_sex = 1:2
        acc_metaMDD(i_diag,i_sex) = 1-sum(double(BSC_metaMDD(cov_MDD.sex==i_sex-1 & cov_MDD.diag==i_diag-1)>0.5)~=(i_sex-1),'omitnan')/sum(cov_MDD.sex==i_sex-1 & cov_MDD.diag==i_diag-1);
    end
end
%%
s = unique(cov_MDD.site);
for i_site = 1:length(s)
    % subject number of XX in each site, where XX represent:
    % 1.all subjects  2.HCs  3.MDDs  4.females  5.males  6.female HCs
    % 7.female MDDs  8.male HCs  9.male MDDs
    n_metaMDD(1,i_site) = sum(cov_MDD.site==s(i_site));
    n_metaMDD(2,i_site) = sum(cov_MDD.site==s(i_site) & cov_MDD.diag==0);
    n_metaMDD(3,i_site) = sum(cov_MDD.site==s(i_site) & cov_MDD.diag==1);
    n_metaMDD(4,i_site) = sum(cov_MDD.site==s(i_site) & cov_MDD.sex==0);
    n_metaMDD(5,i_site) = sum(cov_MDD.site==s(i_site) & cov_MDD.sex==1);
    n_metaMDD(6,i_site) = sum(cov_MDD.site==s(i_site) & cov_MDD.diag==0 & cov_MDD.sex==0);
    n_metaMDD(7,i_site) = sum(cov_MDD.site==s(i_site) & cov_MDD.diag==1 & cov_MDD.sex==0);
    n_metaMDD(8,i_site) = sum(cov_MDD.site==s(i_site) & cov_MDD.diag==0 & cov_MDD.sex==1);
    n_metaMDD(9,i_site) = sum(cov_MDD.site==s(i_site) & cov_MDD.diag==1 & cov_MDD.sex==1);

    % classification accuracy of XX in each site, where XX represent:
    % 1.all subjects  2.HCs  3.MDDs  4.females  5.males  6.female HCs
    % 7.female MDDs  8.male HCs  9.male MDDs
    acc_metaMDD(1,i_site) = 1-sum(double(BSC_metaMDD(cov_MDD.site==s(i_site))>0.5)~=cov_MDD.sex(cov_MDD.site==s(i_site)),'omitnan')/sum(cov_MDD.site==s(i_site));
    acc_metaMDD(2,i_site) = 1-sum(double(BSC_metaMDD(cov_MDD.site==s(i_site) & cov_MDD.diag==0)>0.5)~=cov_MDD.sex(cov_MDD.site==s(i_site) & cov_MDD.diag==0),'omitnan')/sum(cov_MDD.site==s(i_site) & cov_MDD.diag==0);
    acc_metaMDD(3,i_site) = 1-sum(double(BSC_metaMDD(cov_MDD.site==s(i_site) & cov_MDD.diag==1)>0.5)~=cov_MDD.sex(cov_MDD.site==s(i_site) & cov_MDD.diag==1),'omitnan')/sum(cov_MDD.site==s(i_site) & cov_MDD.diag==1);
    acc_metaMDD(4,i_site) = 1-sum(double(BSC_metaMDD(cov_MDD.site==s(i_site) & cov_MDD.sex==0)>0.5)~=0,'omitnan')/sum(cov_MDD.site==s(i_site) & cov_MDD.sex==0);
    acc_metaMDD(5,i_site) = 1-sum(double(BSC_metaMDD(cov_MDD.site==s(i_site) & cov_MDD.sex==1)>0.5)~=1,'omitnan')/sum(cov_MDD.site==s(i_site) & cov_MDD.sex==1);
    acc_metaMDD(6,i_site) = 1-sum(double(BSC_metaMDD(cov_MDD.site==s(i_site) & cov_MDD.diag==0 & cov_MDD.sex==0)>0.5)~=0,'omitnan')/sum(cov_MDD.site==s(i_site) & cov_MDD.diag==0 & cov_MDD.sex==0);
    acc_metaMDD(7,i_site) = 1-sum(double(BSC_metaMDD(cov_MDD.site==s(i_site) & cov_MDD.diag==1 & cov_MDD.sex==0)>0.5)~=0,'omitnan')/sum(cov_MDD.site==s(i_site) & cov_MDD.diag==1 & cov_MDD.sex==0);
    acc_metaMDD(8,i_site) = 1-sum(double(BSC_metaMDD(cov_MDD.site==s(i_site) & cov_MDD.diag==0 & cov_MDD.sex==1)>0.5)~=1,'omitnan')/sum(cov_MDD.site==s(i_site) & cov_MDD.diag==0 & cov_MDD.sex==1);
    acc_metaMDD(9,i_site) = 1-sum(double(BSC_metaMDD(cov_MDD.site==s(i_site) & cov_MDD.diag==1 & cov_MDD.sex==1)>0.5)~=1,'omitnan')/sum(cov_MDD.site==s(i_site) & cov_MDD.diag==1 & cov_MDD.sex==1);

end
% meta-analysis done by R
%% percentage of subjects in 3 BSC groups
clear perct
S = unique(cov_MDD.site);
for i_site = 1:length(S)
    % subject number of females
    nF(1,i_site) = sum(BSC_metaMDD(cov_MDD.sex==0 & cov_MDD.site==S(i_site))<0.35);
    nF(2,i_site) = sum(BSC_metaMDD(cov_MDD.sex==0 & cov_MDD.site==S(i_site))>0.35 & BSC_metaMDD(cov_MDD.sex==0 & cov_MDD.site==S(i_site))<0.65);
    nF(3,i_site) = sum(BSC_metaMDD(cov_MDD.sex==0 & cov_MDD.site==S(i_site))>0.65);

    % subject number of males
    nM(1,i_site) = sum(BSC_metaMDD(cov_MDD.sex==1 & cov_MDD.site==S(i_site))<0.35);
    nM(2,i_site) = sum(BSC_metaMDD(cov_MDD.sex==1 & cov_MDD.site==S(i_site))>0.35 & BSC_metaMDD(cov_MDD.sex==1 & cov_MDD.site==S(i_site))<0.65);
    nM(3,i_site) = sum(BSC_metaMDD(cov_MDD.sex==1 & cov_MDD.site==S(i_site))>0.65);

    % percentage of female HCs in 3 BSC groups
    perct(1,i_site) = sum(BSC_metaMDD(cov_MDD.sex==0 & cov_MDD.diag==0 & cov_MDD.site==S(i_site))<0.35)/sum(cov_MDD.sex==0 & cov_MDD.diag==0 & cov_MDD.site==S(i_site));
    perct(2,i_site) = sum(BSC_metaMDD(cov_MDD.sex==0 & cov_MDD.diag==0 & cov_MDD.site==S(i_site))>0.35 & BSC_metaMDD(cov_MDD.sex==0 & cov_MDD.diag==0 & cov_MDD.site==S(i_site))<0.65)/sum(cov_MDD.sex==0 & cov_MDD.diag==0 & cov_MDD.site==S(i_site));
    perct(3,i_site) = sum(BSC_metaMDD(cov_MDD.sex==0 & cov_MDD.diag==0 & cov_MDD.site==S(i_site))>0.65)/sum(cov_MDD.sex==0 & cov_MDD.diag==0 & cov_MDD.site==S(i_site));

    % percentage of female MDDs in 3 BSC groups
    perct(4,i_site) = sum(BSC_metaMDD(cov_MDD.sex==0 & cov_MDD.diag==1 & cov_MDD.site==S(i_site))<0.35)/sum(cov_MDD.sex==0 & cov_MDD.diag==1 & cov_MDD.site==S(i_site));
    perct(5,i_site) = sum(BSC_metaMDD(cov_MDD.sex==0 & cov_MDD.diag==1 & cov_MDD.site==S(i_site))>0.35 & BSC_metaMDD(cov_MDD.sex==0 & cov_MDD.diag==1 & cov_MDD.site==S(i_site))<0.65)/sum(cov_MDD.sex==0 & cov_MDD.diag==1 & cov_MDD.site==S(i_site));
    perct(6,i_site) = sum(BSC_metaMDD(cov_MDD.sex==0 & cov_MDD.diag==1 & cov_MDD.site==S(i_site))>0.65)/sum(cov_MDD.sex==0 & cov_MDD.diag==1 & cov_MDD.site==S(i_site));

    % percentage of male HCs in 3 BSC groups
    perct(7,i_site) = sum(BSC_metaMDD(cov_MDD.sex==1 & cov_MDD.diag==0 & cov_MDD.site==S(i_site))<0.35)/sum(cov_MDD.sex==1 & cov_MDD.diag==0 & cov_MDD.site==S(i_site));
    perct(8,i_site) = sum(BSC_metaMDD(cov_MDD.sex==1 & cov_MDD.diag==0 & cov_MDD.site==S(i_site))>0.35 & BSC_metaMDD(cov_MDD.sex==1 & cov_MDD.diag==0 & cov_MDD.site==S(i_site))<0.65)/sum(cov_MDD.sex==1 & cov_MDD.diag==0 & cov_MDD.site==S(i_site));
    perct(9,i_site) = sum(BSC_metaMDD(cov_MDD.sex==1 & cov_MDD.diag==0 & cov_MDD.site==S(i_site))>0.65)/sum(cov_MDD.sex==1 & cov_MDD.diag==0 & cov_MDD.site==S(i_site));

    % percentage of male MDDs in 3 BSC groups
    perct(10,i_site) = sum(BSC_metaMDD(cov_MDD.sex==1 & cov_MDD.diag==1 & cov_MDD.site==S(i_site))<0.35)/sum(cov_MDD.sex==1 & cov_MDD.diag==1 & cov_MDD.site==S(i_site));
    perct(11,i_site) = sum(BSC_metaMDD(cov_MDD.sex==1 & cov_MDD.diag==1 & cov_MDD.site==S(i_site))>0.35 & BSC_metaMDD(cov_MDD.sex==1 & cov_MDD.diag==1 & cov_MDD.site==S(i_site))<0.65)/sum(cov_MDD.sex==1 & cov_MDD.diag==1 & cov_MDD.site==S(i_site));
    perct(12,i_site) = sum(BSC_metaMDD(cov_MDD.sex==1 & cov_MDD.diag==1 & cov_MDD.site==S(i_site))>0.65)/sum(cov_MDD.sex==1 & cov_MDD.diag==1 & cov_MDD.site==S(i_site));
end
%% site-wise HAMD correlation
S = unique(cov_MDD.site);
for i_site = 1:length(S)
    % correlation coefficient between HAM-D and BSC in all subjects
    r(1,i_site) = partialcorr(BSC_metaMDD(cov_MDD.site==S(i_site)),cov_MDD.HAMD(cov_MDD.site==S(i_site)),[cov_MDD{cov_MDD.site==S(i_site),[4,5,36]}],'Rows','pairwise');
    % in female MDDs
    r(2,i_site) = partialcorr(BSC_metaMDD(cov_MDD.site==S(i_site) & cov_MDD.sex==0),cov_MDD.HAMD(cov_MDD.site==S(i_site) & cov_MDD.sex==0),[cov_MDD{cov_MDD.site==S(i_site) & cov_MDD.sex==0,[5,36]}],'Rows','pairwise');
    % in male MDDs
    r(3,i_site) = partialcorr(BSC_metaMDD(cov_MDD.site==S(i_site) & cov_MDD.sex==1),cov_MDD.HAMD(cov_MDD.site==S(i_site) & cov_MDD.sex==1),[cov_MDD{cov_MDD.site==S(i_site) & cov_MDD.sex==1,[5,36]}],'Rows','pairwise');
end
% meta-analysis done by R.
%% ANOVA test
BSC_group = zeros(length(BSC_metaMDD),1);
BSC_group(cov_metaMDD.diag==1 & BSC_metaMDD<0.35) = 1; % MDD in female group
BSC_group(cov_metaMDD.diag==1 & BSC_metaMDD>0.35 & BSC_metaMDD<0.65) = 2; % MDD in androgynous group
BSC_group(cov_metaMDD.diag==1 & BSC_metaMDD>0.65) = 3; % MDD in male group
for i_net = 1:12
    for j_net = 1:12
        [~,tbl] = anovan(squeeze(FC_net(i_net,j_net,:)),[BSC_group,cov_metaMDD{:,[4,5,36,3]}],'continuous',[3,4],'display','off');
        % anovan test, controlling 4:sex, 5:age, 36:mean FD, 3:site
        p_anovan(i_net,j_net) = tbl{2,end};
        F(i_net,j_net) = tbl{2,end-1};
    end
end
[ix,iy] = find(p_anovan<0.05/12/11*2); % find significant edges after Bonferroni correction
ind_del = find(ix>iy);
ix(ind_del) = [];
iy(ind_del) = [];
%% bootstrap test to find the difference between BSC groups
clear FC_sig
i_FC = 1;
for i_net = 1:length(ix)
    FC_sig(:,i_FC) = squeeze(FC_net(ix(i_net),iy(i_net),:)); % turning significant FCs into vectors
    i_FC = i_FC+1;
end
S = unique(cov_metaMDD.site);
for i_s = 1:length(S)
    dummycov(:,i_s) = double(cov_metaMDD.site==S(i_s));% site as dummy variable
end
i12 = find(BSC_group==1 | BSC_group==2); % comparing female group and androgynous group
i13 = find(BSC_group==1 | BSC_group==3); % comparing female group and male group
i23 = find(BSC_group==2 | BSC_group==3); % comparing androgynous group and male group
n_bstp = 10000;
for i_bstp = 1:n_bstp
    for i_FC = 1:size(FC_sig,2)
        i_rand = randsample(length(i12),length(i12),1);
        mdl = fitglm([double(BSC_group(i12(i_rand))==2),cov_metaMDD{i12(i_rand),[4,5,36]},dummycov(i12(i_rand),:)],FC_sig(i12(i_rand),i_FC));
        % fitting a glm model to test the difference of FCs between BSC groups.
        % 4:sex, 5:age, 36:meanFD, dummycov:site were used as covariates.
        t(i_FC,i_bstp,1) = mdl.Coefficients.tStat(2); % t-value from glm test

        i_rand = randsample(length(i13),length(i13),1);
        mdl = fitglm([double(BSC_group(i13(i_rand))==3),cov_metaMDD{i13(i_rand),[4,5,36]},dummycov(i13(i_rand),:)],FC_sig(i13(i_rand),i_FC));
        t(i_FC,i_bstp,2) = mdl.Coefficients.tStat(2);

        i_rand = randsample(length(i23),length(i23),1);
        mdl = fitglm([double(BSC_group(i23(i_rand))==3),cov_metaMDD{i23(i_rand),[4,5,36]},dummycov(i23(i_rand),:)],FC_sig(i23(i_rand),i_FC));
        t(i_FC,i_bstp,3) = mdl.Coefficients.tStat(2);
    end
    fprintf('boostrapping %.f\n',i_bstp)
end

for i_group = 1:3
    for i_FC = 1:size(FC_sig,2)
        t(i_FC,:,i_group) = sort(t(i_FC,:,i_group),'ascend'); % sorting t values, so to check whether the difference is significant
    end
end
%% comparing MDDs in 3 BSC groups with HC
clear FC_sig
i_FC = 1;
for i_net = 1:length(ix)
    FC_sig(:,i_FC) = squeeze(FC_net(ix(i_net),iy(i_net),:));% turning significant FCs into vectors
    i_FC = i_FC+1;
end
S = unique(cov_metaMDD.site);
for i_s = 1:length(S)
    dummycov(:,i_s) = double(cov_metaMDD.site==S(i_s));% site as dummy variable
end
clear t;clear p;
for i_bsc_group = 1:3
    for i_FC = 1:size(FC_sig,2)
        mdl = fitglm([cov_metaMDD{cov_metaMDD.diag==0| BSC_group==i_bsc_group,[2,4,5,36]},dummycov(cov_metaMDD.diag==0| BSC_group==i_bsc_group,:)],FC_sig(cov_metaMDD.diag==0| BSC_group==i_bsc_group,i_FC));
        % fitting a glm model to test the difference of FCs between HC and MDD in respective BSC groups.
        % 4:sex, 5:age, 36:meanFD, dummycov:site were used as covariates.
        t(i_FC,i_bsc_group) = mdl.Coefficients.tStat(2);
        p(i_FC,i_bsc_group) = mdl.Coefficients.pValue(2);
    end
end
%% FC t-value plot (Figure 2 e)
p([5,7,9,11,12,13],:) = []; % deleting edges without group differences from bootstrap
t([5,7,9,11,12,13],:) = [];
plot((-log(p).*sign(t))')
hold on
plot([0.5,3.5],[-log(0.05/21),-log(0.05/21)],'black') % multiple-comparison among 7edges * 3groups
plot([0.5,3.5],[log(0.05/21),log(0.05/21)],'black')
plot([0.5,3.5],[0,0],'black')
xlim([0.5,3.5])
legend('DMN-DMN','DMN-MRN','DMN-VN','SMN(H)-SCN','DMN-SCN','VN-SCN','SN-SCN','threshold')
set(gca,'xtick',[1:3],'xticklabel',{'female-like','androgynous','male-like'})