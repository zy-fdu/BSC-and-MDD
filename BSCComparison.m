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

load('dep_TC_aal.mat');%% timecourse, dep_TC
load('cov_dep.mat'); %% cov_dep
load('svmstruct_UKB.mat');

loc_8w=strfind(cov_dep.ID,'8w');
index_8w=~cellfun('isempty',loc_8w);

index_0w=ones(size(cov_dep.ID,1),1);
index_0w=logical(index_0w-index_8w);

cov_dep_0w=cov_dep(index_0w,:);
cov_dep_8w=cov_dep(index_8w,:);

dep_TC_0w=dep_TC(index_0w);
dep_TC_8w=dep_TC(index_8w);

%% predict and comput BSC
[cov_dep_0w.BSC_0w,cov_dep_0w.output_0w]=SVMprediction(dep_TC_0w,cov_dep_0w,svmstruct,1);
[cov_dep_8w.BSC_8w,cov_dep_8w.output_8w]=SVMprediction(dep_TC_8w,cov_dep_8w,svmstruct,'glmstruct.mat');
%% regresse the covariates: age sex and FD
cov_dep_0w.BSC_0w_regressed=myregression(cov_dep_0w.BSC_0w,[cov_dep_0w.age,cov_dep_0w.sex,cov_dep_0w.FD]);
cov_dep_8w.BSC_8w_regressed=myregression(cov_dep_8w.BSC_8w,[cov_dep_8w.age,cov_dep_8w.sex,cov_dep_8w.FD]);

%% extract subjects with both 0w and 8w
subID_8w=cellfun(@(x)strrep(x,'_8w',''),cov_dep_8w.ID,'UniformOutput',false);
subID_0w=cellfun(@(x)strrep(x,'_0w',''),cov_dep_0w.ID,'UniformOutput',false);
[idx_0w,idx_8w]=ismember(subID_0w,subID_8w);
idx_8w(find(idx_8w==0))=[];
cov_dep_0w_8w=cov_dep_0w(idx_0w,:);
cov_dep_0w_8w_8w=cov_dep_8w(idx_8w,:);
cov_dep_0w_8w.FD_8w=cov_dep_0w_8w_8w.FD;
cov_dep_0w_8w.age_8w=cov_dep_0w_8w_8w.age;
cov_dep_0w_8w.BSC_8w=cov_dep_0w_8w_8w.BSC_8w;
cov_dep_0w_8w.BSC_8w_regressed=cov_dep_0w_8w_8w.BSC_8w_regressed;
cov_dep_0w_8w(cov_dep_0w_8w.HAMD_8w==-1 & cov_dep_0w_8w.groupBinary,:)=[]; %%remove the subjects without HAMD_8w

%% paired T three groups 0.35-0.65
index_35=find(cov_dep_0w_8w.BSC_0w(cov_dep_0w_8w.groupBinary)<0.35);
index_35_65=find(cov_dep_0w_8w.BSC_0w(cov_dep_0w_8w.groupBinary)>=0.35&cov_dep_0w_8w.BSC_0w(cov_dep_0w_8w.groupBinary)<=0.65);
index_65=find(cov_dep_0w_8w.BSC_0w(cov_dep_0w_8w.groupBinary)>0.65);

cov_dep_0w_8w.BSCgroup(index_35)=1;
cov_dep_0w_8w.BSCgroup(index_35_65)=2;
cov_dep_0w_8w.BSCgroup(index_65)=3;

deltaHAMD=cov_dep_0w_8w.HAMD_0w-cov_dep_0w_8w.HAMD_8w;

[h_35,p_35,CI,STATS]=ttest(cov_dep_0w_8w.BSC_0w(index_35),cov_dep_0w_8w.BSC_8w(index_35));
[h_35_65,p_35_65,CI,STATS]=ttest(cov_dep_0w_8w.BSC_0w(index_35_65),cov_dep_0w_8w.BSC_8w(index_35_65));
[h_65,p_65,CI,STATS]=ttest(cov_dep_0w_8w.BSC_0w(index_65),cov_dep_0w_8w.BSC_8w(index_65));

[h_35_regressed,p_35_regressed]=ttest(cov_dep_0w_8w.BSC_0w_regressed(index_35),cov_dep_0w_8w.BSC_8w_regressed(index_35));
[h_35_65_regressed,p_35_65_regressed]=ttest(cov_dep_0w_8w.BSC_0w_regressed(index_35_65),cov_dep_0w_8w.BSC_8w_regressed(index_35_65));
[h_65_regressed,p_65_regressed]=ttest(cov_dep_0w_8w.BSC_0w_regressed(index_65),cov_dep_0w_8w.BSC_8w_regressed(index_65));

display_boxplot_parallel_BSC(cov_dep_0w_8w);
display_parallel_HAMD(cov_dep_0w_8w)
display_deltaHAMD3groups(deltaHAMD,cov_dep_0w_8w.BSCgroup,outdir);

%% display by different medicine
clear group;
group(:,1)=cov_dep_0w_8w.BSCgroup;
group(:,2)=cov_dep_0w_8w.medicine;
m=1;

for i=1:8  %eight medicines
    if ~isempty(find(cov_dep_0w_8w.medicine==i))
        HAMD_0w=cov_dep_0w_8w.HAMD_0w(find(cov_dep_0w_8w.medicine==i));
        HAMD_8w=cov_dep_0w_8w.HAMD_8w(find(cov_dep_0w_8w.medicine==i));
        BSC_0w=cov_dep_0w_8w.BSC_0w(find(cov_dep_0w_8w.medicine==i));
        BSC_8w=cov_dep_0w_8w.BSC_8w(find(cov_dep_0w_8w.medicine==i));
        
        if numel(unique(cov_dep_0w_8w.BSCgroup(find(cov_dep_0w_8w.medicine==i))))==3 %% at least on subject in each BSC group
            outname=[outdir,'medicine_',num2str(i)];
            display_0w_8w_mean(HAMD_0w,HAMD_8w,cov_dep_0w_8w.BSCgroup(find(cov_dep_0w_8w.medicine==i)),outname,'HAMD'); %%display HAMD by medicine
            [p(m),t{m},stats{m}]=anovan(deltaHAMD(find(cov_dep_0w_8w.medicine==i)),group(find(cov_dep_0w_8w.medicine==i),1));
            c_HAMD(m)={multcompare(stats)}; m=m+1;
            
            outname=[outdir,'medicine_',num2str(i)];
            display_0w_8w_mean(BSC_0w,BSC_8w,cov_dep_0w_8w.BSCgroup(find(cov_dep_0w_8w.medicine==i)),outname,'BSC');  %%display BSC by medicine         
        else
            disp([outname,', ',num2str(numel(find(cov_dep_0w_8w.medicine==i)))]);
        end
    end
    save([outdir,'HAMDcomparison_bymedicine.mat'],'p','t','stats','c_HAMD');
end

diary off;