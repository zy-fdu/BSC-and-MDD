%% PKU remission
BSC_group = double(cov_PKU_0w_8w.diag==1 & BSC<0.35)+2*double(cov_PKU_0w_8w.diag==1 & BSC>0.35 & BSC<0.65)+3*double(cov_PKU_0w_8w.diag==1 & BSC>0.65);
remission = double(cov_PKU_8w.HAMD<=7);
[~,chi2,p] = crosstab(BSC_group,remission); % cross-tab chi2 test

% plot of figure 3f
clear r
for i_bsc = 1:3
    r(i_bsc) = sum(BSC_group==i_bsc & remission==0)/sum(BSC_group==i_bsc);
end
barh([3:-1:1],[1;1;1])
hold on
barh([3:-1:1],[r])
set(gca,'YTickLabel',{})
ylim([0.2,3.8])