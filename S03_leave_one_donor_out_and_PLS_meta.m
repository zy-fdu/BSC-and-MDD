
% Function: do leave-one-donor-out analysis and plot.
% Updata date: 2022.07.05
% Email:luolongcao@163.com
clear;
clc;

cd('F:\work_dir\AHBAenrich/');
work_dir ='F:/work_dir/AHBAenrich/file_meta/leave_one_sub_MHC/';
load('AHBA_Mean_scaled_reanote.mat'); 
load ('AHBA_data_reanote.mat');
load('AHBA_ROI_index_80.mat');

% description of the data
% Coordinatesall       ---- MNI Coordinate of all 3702 samples
% cort_expMeanScaled   ---- all gene expression data of the cortical tissues 
%                           from 6 donors
% cort_l_expMS         ---- 1105 cortical tissues for 
%                           15408 genes
% geneSymbol           ---- gene names  
% genes                ---- entrez gene id
% sub_expMeanScaled    ---- all gene expression data of the subcortical tissues 
%                           from 6 donors 
% sub_l_expMS          ---- 645 subcortical tissues in the left hemisphere for 
%                           15408 genes
% Ncort_l_80           ---- sample number in cortex and correspond properity
% Nsub_l_80            ---- sample number in subcortex and correspond properity
%%
Yraw = readtable('F:\work_dir\AHBAenrich\file_meta\meta corr.csv'); 
Yraw_new(:,5:8) = zscore(table2array(Yraw(:,5:8)));
Yraw.Properties.VariableNames
corr(table2array(Yraw(:,5:8)))

XYZ = [Yraw.MNI_X,Yraw.MNI_Y,Yraw.MNI_Z];
XYZ_lables = Yraw.ROI';
len = size(XYZ_lables,2);

ROI_radius =5;%set a radius of a sphere to coordinata probes and ROIs.
thr = [0,946,1839,2202,2731,3201,3702];
%%
%LOO (leave one donor out) and redo PLS.
clc
t=cputime;
feature_column = 6;%column 5~8 features:{'FHC'}    {'FMDD'}    {'MHC'}    {'MMDD'}

%cortical
B = Coordinatesall(Ncort_l_80(:,1),:);
cort_l_expMS= [zscore(cort_l_expMS)];
cort_l_expMS_new = [cort_l_expMS];
% %subcortical
% B = Coordinatesall(Nsub_l_80(:,1),:);
% sub_l_expMS = [zscore(sub_l_expMS)];
% cort_l_expMS_new = [sub_l_expMS];
% % whole brain
% B = Coordinatesall([Ncort_l_80(:,1);Nsub_l_80(:,1)],:);
% cort_l_expMS = [zscore(cort_l_expMS;sub_l_expMS)];
% cort_l_expMS_new = [cort_l_expMS];

N = size(B,1);
all_roi_exp = [];
in_rois = [];
for i = 1:len
    A = XYZ(i,:);
    distance = sqrt(sum((B-repmat(A,N,1)).^2,2));% calculate distance of sample A and all point in B
    temp = find(distance<= ROI_radius);
    if ~isempty(temp)
        all_roi_exp(i,:) = mean(cort_l_expMS_new(temp,:),1);
        in_rois = [in_rois,i];
    else
        all_roi_exp(i,:) = repelem(0,size(cort_l_expMS_new,2));
    end
end

X = all_roi_exp(in_rois,:);
Y = Yraw_new(in_rois,feature_column);%column 5~8 features:{'FHC'}    {'FMDD'}    {'MHC'}    {'MMDD'}
dim=5;
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim);
component_column = find(PCTVAR(2,:) ==max(PCTVAR(2,:)));
weight_whole = stats.W(:,component_column);
corr_LOO = [1];
variance_LOO = [PCTVAR(2,component_column)];

for leave_id=1:6
    
    leave_id
    
%     cortical
    leave_one = find([Ncort_l_80(:,1)>thr(leave_id)] .* [Ncort_l_80(:,1)<=thr(leave_id+1)]);%selece samples of specific donor.
    Ncort_l_80_new = Ncort_l_80;
    Ncort_l_80_new(leave_one,:) = [];
    B = Coordinatesall(Ncort_l_80_new(:,1),:);
    cort_l_expMS_new = [cort_l_expMS];
    cort_l_expMS_new(leave_one,:) = [];
    
% %     subcortical
%     leave_one = find([Nsub_l_80(:,1)>thr(leave_id)] .* [Nsub_l_80(:,1)<=thr(leave_id+1)]);
%     Ncort_l_80_new = Nsub_l_80;
%     Ncort_l_80_new(leave_one,:) = [];
%     B = Coordinatesall(Nsub_l_80(:,1),:);
%     cort_l_expMS_new = [sub_l_expMS];
%     cort_l_expMS_new(leave_one,:) = [];
    
% %     whole brain
%     Nwhole_l_80 = [Ncort_l_80(:,1);Nsub_l_80(:,1)];
%     leave_one = find([Nwhole_l_80(:,1)>thr(leave_id)] .* [Nwhole_l_80(:,1)<=thr(leave_id+1)]);
%     Ncort_l_80_new = Nwhole_l_80;
%     Ncort_l_80_new(leave_one,:) = [];
%     B = Coordinatesall(Ncort_l_80_new(:,1),:);
%     cort_l_expMS_new = [cort_l_expMS;sub_l_expMS];
%     cort_l_expMS_new(leave_one,:) = [];

    N = size(B,1);
    all_roi_exp = [];
    in_rois = [];
    for i = 1:len
        A = XYZ(i,:);
        distance = sqrt(sum((B-repmat(A,N,1)).^2,2));% calculate distance of sample A and all point in B
        temp = find(distance<= ROI_radius);
        if ~isempty(temp)
            all_roi_exp(i,:) = mean(cort_l_expMS_new(temp,:),1);
            in_rois = [in_rois,i];
        else
            all_roi_exp(i,:) = repelem(0,size(cort_l_expMS_new,2));
        end
    end

    X = all_roi_exp(in_rois,:);
    Y = Yraw_new(in_rois,feature_column);%column 5~8 features:{'FHC'}    {'FMDD'}    {'MHC'}    {'MMDD'}

%     dim=5;
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim);
    
    weight_temp = stats.W(:,component_column);
    corr_LOO = [corr_LOO, corr(weight_whole,weight_temp)];
    variance_LOO = [variance_LOO,PCTVAR(2,component_column)];
    
%     subplot(2,6,leave_id);
% %     fig=figure
%     plot(1:dim,cumsum(100*PCTVAR(2,1:dim)),'-o','LineWidth',1.5,'Color',[140/255,0,0]);
%     set(gca,'Fontsize',14)
%     str_xlab=strcat('Number of PLS components of donor_',num2str(leave_id));
%     xlabel(str_xlab,'FontSize',14);
%     ylabel('Percent Variance Explained in Y','FontSize',14);
%     text(1:dim,cumsum(100*PCTVAR(2,1:dim))-3,num2str(cumsum(100*PCTVAR(2,1:dim)).','%.1f'))
%     grid on
% %     fig_name = strcat(work_dir,'leave_',num2str(leave_id),'_PLS_variance.fig');
% %     savefig(fig,fig_name)
% %     size(in_rois,2)
% %     PCTVAR(2,7)    

%     rep=500;
%     dim = 10;
%     [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim);
%     Rsquared_real_y = cumsum(100*PCTVAR(2,1:dim));
%     Rsq_y = [];
%     Pvalue_perm_y = ones(1,dim);
%     for i = 1:dim
%         parfor j = 1:rep % automatic parallel calculate
%         order=randperm(size(Y,1));
%         Yp=Y(order,:);
%         [XL_p,YL_p,XS_p,YS_p,BETA_p,PCTVAR_p,MSE_p,stats_p]=plsregress(X,Yp,dim);
%         temp_y=cumsum(100*PCTVAR_p(2,1:dim));
%         Rsq_y(i,j) = temp_y(i);
%         end
%         Pvalue_perm_y(i) = length(find(Rsq_y(i,:)>=Rsquared_real_y(i)))/rep;
%     end
% 
%     subplot(2,6,6+leave_id);
% %     fig = figure;
%     plot(1:dim, Pvalue_perm_y,'ok','MarkerSize',8,'MarkerFaceColor','r');
%     xlabel('Number of PLS components','FontSize',14);
%     ylabel('p-value','FontSize',14);
%     text(1:dim,Pvalue_perm_y,num2str(Pvalue_perm_y.','%.3f'))
%     grid on    
% %     fig_name = strcat(work_dir,'leave_',num2str(leave_id),'_PLS_perm.fig');
% %     savefig(fig,fig_name)
%     % img = frame2in(getframe(fig));
%     % imwrite(img,'test.jpg');

end
corr_LOO
variance_LOO
exp_time = cputime-t;

% plot of figure with 2 y-axis.
fig = figure;
yyaxis left
bar(1:7, variance_LOO);
xlabel('Leave one donor out','FontSize',14);
ylabel('Percent Variance Explained in Y','FontSize',14);
text(1:7,variance_LOO,num2str(variance_LOO.','%.3f'))
yyaxis right
plot(1:7, corr_LOO,'-o','MarkerSize',8,'MarkerFaceColor','r');
ylabel('Correlation of PLS component 2 weight','FontSize',14);
% text(1:7,corr_LOO,num2str(corr_LOO.','%.3f'))
grid on 
% fig_name = strcat(work_dir,'leave_',num2str(leave_id),'_PLS_perm.fig');
% savefig(fig,fig_name)
% img = frame2in(getframe(fig));
% imwrite(img,'test.jpg');