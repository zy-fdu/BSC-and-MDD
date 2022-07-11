
% Function: do partial least squares analysis (permutation, bootstrap, get genes' z-score) and plot.
% Updata date: 2022.07.05
% Email:luolongcao@163.com
clear;
clc;

cd('F:\work_dir\AHBAenrich/');
work_dir ='F:/work_dir/AHBAenrich/file_meta';
load('AHBA_Mean_scaled_reanote.mat'); 
load ('AHBA_data_reanote.mat');
load('AHBA_ROI_index_80.mat');

% description of the data
% cort_expMeanScaled   ---- all gene expression data of the cortical tissues 
%                           from 6 donors*42 ROIs 
% cort_l_expMS         ---- 1253 cortical tissues in the 42 ROIs for 
%                           15408 genes
% geneSymbol           ---- gene names  
% genes                ---- entrez gene id
% sub_expMeanScaled    ---- all gene expression data of the subcortical tissues 
%                           from 6 donors 
% sub_l_expMS          ---- 182 subcortical tissues in the left hemisphere for 
%                           15408 genes
%%
Yraw = readtable('F:\work_dir\AHBAenrich\file_meta\meta corr.csv'); 
Yraw.Properties.VariableNames
Yraw_new(:,5:8) = zscore(table2array(Yraw(:,5:8)));
corr(table2array(Yraw(:,5:8)))

XYZ = [Yraw.MNI_X,Yraw.MNI_Y,Yraw.MNI_Z];
XYZ_lables = Yraw.ROI';

ROI_radius =5;%set a radius of a sphere to coordinata probes and ROIs.
len = size(XYZ_lables,2);
% %cortical
% B = Coordinatesall(Ncort_l_80(:,1),:);
% cort_l_expMS_new = [zscore(cort_l_expMS)];
% %subcortical
% B = Coordinatesall(Nsub_l_80(:,1),:);
% cort_l_expMS_new = [zscore(sub_l_expMS)];
% whole brain
B = Coordinatesall([Ncort_l_80(:,1);Nsub_l_80(:,1)],:);
cort_l_expMS_new = zscore([cort_l_expMS;sub_l_expMS]);

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
%%
for feature = 1:4
    Y = Yraw_new(in_rois,4+feature);%column 5~8 features:{'FHC'}    {'FMDD'}    {'MHC'}    {'MMDD'}
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%picture start
    %perform full PLS and plot variance in Y explained by top 15 components
    %typically top 2 or 3 components will explain a large part of the variance
    %(hopefully!)
    dim=5;
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim);
    component_column = find(PCTVAR(2,:) ==max(PCTVAR(2,:)));

    %PCTVAR containing the percentage of variance explained by the model. 
    %The first row of PCTVAR contains the percentage of variance explained in X by each PLS component, 
    %and the second row contains the percentage of variance explained in Y
    subplot(2,4,feature);
%     figure
    plot(1:dim,cumsum(100*PCTVAR(2,1:dim)),'-o','LineWidth',1.5,'Color',[140/255,0,0]);
    set(gca,'Fontsize',14)
    str_xlab=strcat('Number of PLS components of_ ',Yraw.Properties.VariableNames{4+feature});
    xlabel(str_xlab,'FontSize',14);
    ylabel('Percent Variance Explained in Y','FontSize',14);
    text(1:dim,cumsum(100*PCTVAR(2,1:dim))-3,num2str(cumsum(100*PCTVAR(2,1:dim)).','%.1f'))
    grid on

%     figure
%     plot(1:dim,cumsum(100*PCTVAR(1,1:dim)),'-o','LineWidth',1.5,'Color',[140/255,0,0]);
%     set(gca,'Fontsize',14)
%     xlabel('Number of PLS components','FontSize',14);
%     ylabel('Percent Variance Explained in X','FontSize',14);
%     grid on
%     
%     % plot correlation of PLS significant component with Y:
%     figure
%     plot(XS(:,component_column),Y,'r.')% the predictor scores XS, that is, the PLS components that are linear combinations of the variables in X
%     [R,P] = corrcoef(XS(:,component_column),Y)
%     temp_str = strcat('XS scores for PLS component ',num2str(component_column));
%     xlabel(temp_str,'FontSize',14);
%     ylabel('Y','FontSize',14);
%     grid on
%     l = lsline;
%     l.Color = 'k';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%picture end
    %
    % permutation testing to assess significance of PLS result as a function of
    % the number of components (dim) included:
    rep=1000;
    dim = 5;
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim);
    PCTVAR
    Rsquared_real_y = cumsum(100*PCTVAR(2,1:dim));
    Rsq_y = [];
    Pvalue_perm_y = ones(1,dim);
    for i = 1:dim
        parfor j = 1:rep % automatic parallel calculate
        order=randperm(size(Y,1));
        Xp=X(order,:);
        [XL_p,YL_p,XS_p,YS_p,BETA_p,PCTVAR_p,MSE_p,stats_p]=plsregress(Xp,Y,dim);

        temp_y=cumsum(100*PCTVAR_p(2,1:dim));
        Rsq_y(i,j) = temp_y(i);
        end
        Pvalue_perm_y(i) = length(find(Rsq_y(i,:)>=Rsquared_real_y(i)))/rep;
    end

    subplot(2,4,4+feature);
%     fig = figure;
    plot(1:dim, Pvalue_perm_y,'ok','MarkerSize',8,'MarkerFaceColor','r');
    xlabel(str_xlab,'FontSize',14);
    ylabel('p-value','FontSize',14);
    text(1:dim,Pvalue_perm_y,num2str(Pvalue_perm_y.','%.3f'))
    grid on
    % img = frame2in(getframe(fig));
    % imwrite(img,'test.jpg');
end
%%
%Bootstrap to get the gene list:
%in order to estimate the error in estimating each gene's PLS1 weight, then
%to caculate the z score
%number of bootstrap iterations:
clc
bootnum=5000;
Y = Yraw_new(in_rois,6);%column 5~8 features:{'FHC'}    {'FMDD'}    {'MHC'}    {'MMDD'}
% Do PLS in 2 dimensions (with 2 components):
dim=5;
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim);

temp_dim = find(PCTVAR(2,:) ==max(PCTVAR(2,:)));
%store regions' IDs and weights in descending order of weight for both components:
[R1,p1]=corr([XS(:,1),XS(:,2),XS(:,temp_dim)],Y);%XS: predicted X score

%align PLS components with desired direction for interpretability 
if R1(1,1)<0  %this is specific to the data shape we were using - will need ammending
    stats.W(:,1)=-1*stats.W(:,1);
    XS(:,1)=-1*XS(:,1);
end
if R1(2,1)<0 %this is specific to the data shape we were using - will need ammending
    stats.W(:,2)=-1*stats.W(:,2);
    XS(:,2)=-1*XS(:,2);
end
if R1(3,1)<0 %this is specific to the data shape we were using - will need ammending
    stats.W(:,temp_dim)=-1*stats.W(:,temp_dim);
    XS(:,temp_dim)=-1*XS(:,temp_dim);
end
% save PLS_XandY X Y
%real weight and sorted by weight
[PLS1w,x1] = sort(stats.W(:,1),'descend');%W: A p-by-ncomp matrix of PLS weights so that XS = X0*W.
PLS1ids=genes(x1);
[PLS2w,x2] = sort(stats.W(:,2),'descend');
PLS2ids=genes(x2);
[PLS3w,x3] = sort(stats.W(:,temp_dim),'descend');
PLS3ids=genes(x3);

%define variables for storing the (ordered) weights from all bootstrap runs
PLS1weights=[];
PLS2weights=[];
PLS3weights=[];

%start bootstrap
parfor i=1:bootnum
    myresample = randsample(size(X,1),size(X,1),1);%replacement
    res(i,:)=myresample; %store resampling out of interest
    Xr=X(myresample,:); % define X for resampled subjects
    Yr=Y(myresample,:); % define X for resampled subjects
    [XL_b,YL_b,XS_b,YS_b,BETA_b,PCTVAR_b,MSE_b,stats_b]=plsregress(Xr,Yr,dim); %perform PLS for resampled data
      
    temp1=stats_b.W(:,1);%extract PLS1 weights
    newW1=temp1(x1); %order the newly obtained weights the same way as initial PLS 
    if corr(PLS1w,newW1)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
        newW1=-1*newW1;
    end
    PLS1weights=[PLS1weights,newW1];%store (ordered) weights from this bootstrap run
    
    temp2=stats_b.W(:,2);%extract PLS2 weights
    newW2=temp2(x2); %order the newly obtained weights the same way as initial PLS 
    if corr(PLS2w,newW2)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
        newW2=-1*newW2;
    end
    PLS2weights=[PLS2weights,newW2];%store (ordered) weights from this bootstrap run
    
    temp3=stats_b.W(:,temp_dim);%extract PLS2 weights
    newW3=temp3(x3); %order the newly obtained weights the same way as initial PLS 
    if corr(PLS3w,newW3)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
        newW3=-1*newW3;
    end
    PLS3weights=[PLS3weights,newW3]; %store (ordered) weights from this bootstrap run    

end

%get standard deviation of weights from bootstrap runs
PLS1sw=std(PLS1weights');
PLS2sw=std(PLS2weights');
PLS3sw=std(PLS3weights');

%get bootstrap weights (Z)
PLS1Z=PLS1w./PLS1sw';
PLS2Z=PLS2w./PLS2sw';
PLS3Z=PLS3w./PLS3sw';

PLS_out_genes = table(genes, geneSymbol, PLS1Z);
PLS2_out_genes = table(genes, geneSymbol, PLS2Z);
PLS3_out_genes = table(genes, geneSymbol, PLS3Z);

cd (work_dir)
writetable(sortrows(PLS_out_genes,-3),'PLS_out_genes_F_MDD_cort_PLS_variance1.csv');
writetable(sortrows(PLS2_out_genes,-3),'PLS_out_genes_F_MDD_cort_PLS_variance2.csv');
writetable(sortrows(PLS3_out_genes,-3),'PLS_out_genes_F_MDD_cort_PLS_variance_max.csv');
cd '..'


