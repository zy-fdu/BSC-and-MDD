clc;
clear;
close all;
%% 
outdir=['.\result\',date];

if ~exist(outdir,'dir')
    mkdir(outdir);
end

logName=['./log/BSCComparison_log',date];
logFiles=dir([logName,'*.txt']);
if ~isempty(logFiles)
    logName=[logName,num2str(length(logFiles)+1)];
end
%    
diary([logName,'.txt']);

load('RSts_aal2_PKU_aal.mat');%% timecourse, RSts_aal2_PKU
load('cov_PKU.mat'); %% cov_PKU
load('svmstruct_UKB.mat');

loc_8w=strfind(cov_PKU.ID,'8w');
index_8w=~cellfun('isempty',loc_8w);

index_0w=ones(size(cov_PKU.ID,1),1);
index_0w=logical(index_0w-index_8w);

cov_PKU_0w=cov_PKU(index_0w,:);
cov_PKU_8w=cov_PKU(index_8w,:);

RSts_aal2_PKU_0w=RSts_aal2_PKU(index_0w);
RSts_aal2_PKU_8w=RSts_aal2_PKU(index_8w);

%% sex classification and compute BSC
[cov_PKU_0w.BSC_0w,cov_PKU_0w.output_0w]=SVMprediction(RSts_aal2_PKU_0w,cov_PKU_0w,svmstruct,1);
[cov_PKU_8w.BSC_8w,cov_PKU_8w.output_8w]=SVMprediction(RSts_aal2_PKU_8w,cov_PKU_8w,svmstruct,'glmstruct.mat');
%% regress out the covariates: age sex and FD
cov_PKU_0w.BSC_0w_regressed=myregression(cov_PKU_0w.BSC_0w,[cov_PKU_0w.age,cov_PKU_0w.sex,cov_PKU_0w.FD]);
cov_PKU_8w.BSC_8w_regressed=myregression(cov_PKU_8w.BSC_8w,[cov_PKU_8w.age,cov_PKU_8w.sex,cov_PKU_8w.FD]);

%% extract subjects with both 0w and 8w
subID_8w=cellfun(@(x)strrep(x,'_8w',''),cov_PKU_8w.ID,'UniformOutput',false);
subID_0w=cellfun(@(x)strrep(x,'_0w',''),cov_PKU_0w.ID,'UniformOutput',false);
[idx_0w,idx_8w]=ismember(subID_0w,subID_8w);
idx_8w(find(idx_8w==0))=[];
cov_PKU_0w_8w=cov_PKU_0w(idx_0w,:);
cov_PKU_0w_8w_8w=cov_PKU_8w(idx_8w,:);
cov_PKU_0w_8w.FD_8w=cov_PKU_0w_8w_8w.FD;
cov_PKU_0w_8w.age_8w=cov_PKU_0w_8w_8w.age;
cov_PKU_0w_8w.BSC_8w=cov_PKU_0w_8w_8w.BSC_8w;
cov_PKU_0w_8w.BSC_8w_regressed=cov_PKU_0w_8w_8w.BSC_8w_regressed;
cov_PKU_0w_8w(cov_PKU_0w_8w.HAMD_8w==-1 & cov_PKU_0w_8w.groupBinary,:)=[]; %%remove the subjects without HAMD_8w

%% paired T three groups 0.35-0.65
index_35=find(cov_PKU_0w_8w.BSC_0w(cov_PKU_0w_8w.groupBinary)<0.35);
index_35_65=find(cov_PKU_0w_8w.BSC_0w(cov_PKU_0w_8w.groupBinary)>=0.35&cov_PKU_0w_8w.BSC_0w(cov_PKU_0w_8w.groupBinary)<=0.65);
index_65=find(cov_PKU_0w_8w.BSC_0w(cov_PKU_0w_8w.groupBinary)>0.65);

cov_PKU_0w_8w.BSCgroup(index_35)=1;
cov_PKU_0w_8w.BSCgroup(index_35_65)=2;
cov_PKU_0w_8w.BSCgroup(index_65)=3;

deltaHAMD=cov_PKU_0w_8w.HAMD_0w-cov_PKU_0w_8w.HAMD_8w;

[h_35,p_35,CI,STATS]=ttest(cov_PKU_0w_8w.BSC_0w(index_35),cov_PKU_0w_8w.BSC_8w(index_35));
[h_35_65,p_35_65,CI,STATS]=ttest(cov_PKU_0w_8w.BSC_0w(index_35_65),cov_PKU_0w_8w.BSC_8w(index_35_65));
[h_65,p_65,CI,STATS]=ttest(cov_PKU_0w_8w.BSC_0w(index_65),cov_PKU_0w_8w.BSC_8w(index_65));

[h_35_regressed,p_35_regressed]=ttest(cov_PKU_0w_8w.BSC_0w_regressed(index_35),cov_PKU_0w_8w.BSC_8w_regressed(index_35));
[h_35_65_regressed,p_35_65_regressed]=ttest(cov_PKU_0w_8w.BSC_0w_regressed(index_35_65),cov_PKU_0w_8w.BSC_8w_regressed(index_35_65));
[h_65_regressed,p_65_regressed]=ttest(cov_PKU_0w_8w.BSC_0w_regressed(index_65),cov_PKU_0w_8w.BSC_8w_regressed(index_65));

display_boxplot_parallel_BSC(cov_PKU_0w_8w);
display_parallel_HAMD(cov_PKU_0w_8w)
display_deltaHAMD3groups(deltaHAMD,cov_PKU_0w_8w.BSCgroup,outdir);

%% display by different medicine
clear group;
group(:,1)=cov_PKU_0w_8w.BSCgroup;
group(:,2)=cov_PKU_0w_8w.medicine;
m=1;

for i=1:8  %eight medicines
    if ~isempty(find(cov_PKU_0w_8w.medicine==i))
        HAMD_0w=cov_PKU_0w_8w.HAMD_0w(find(cov_PKU_0w_8w.medicine==i));
        HAMD_8w=cov_PKU_0w_8w.HAMD_8w(find(cov_PKU_0w_8w.medicine==i));
        BSC_0w=cov_PKU_0w_8w.BSC_0w(find(cov_PKU_0w_8w.medicine==i));
        BSC_8w=cov_PKU_0w_8w.BSC_8w(find(cov_PKU_0w_8w.medicine==i));
        
        if numel(unique(cov_PKU_0w_8w.BSCgroup(find(cov_PKU_0w_8w.medicine==i))))==3 %% at least on subject in each BSC group
            outname=[outdir,'medicine_',num2str(i)];
            display_0w_8w_mean(HAMD_0w,HAMD_8w,cov_PKU_0w_8w.BSCgroup(find(cov_PKU_0w_8w.medicine==i)),outname,'HAMD'); %%display HAMD by medicine
            [p(m),t{m},stats{m}]=anovan(deltaHAMD(find(cov_PKU_0w_8w.medicine==i)),group(find(cov_PKU_0w_8w.medicine==i),1));
            c_HAMD(m)={multcompare(stats)}; m=m+1;
            
            outname=[outdir,'medicine_',num2str(i)];
            display_0w_8w_mean(BSC_0w,BSC_8w,cov_PKU_0w_8w.BSCgroup(find(cov_PKU_0w_8w.medicine==i)),outname,'BSC');  %%display BSC by medicine         
        else
            disp([outname,', ',num2str(numel(find(cov_PKU_0w_8w.medicine==i)))]);
        end
    end
    save([outdir,'HAMDcomparison_bymedicine.mat'],'p','t','stats','c_HAMD');
end

diary off;