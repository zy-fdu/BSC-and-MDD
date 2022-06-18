function [BSC,output]=SVMprediction(RSts_aal2_PKU,cov_PKU,svmstruct,glmstruct_path)
%%
n_ROI = 94;

if strcmp(class(glmstruct_path),'char')
    glmstruct=importdata(glmstruct_path);
    glm_compute=0;
else
    glmstruct=cell((n_ROI*(n_ROI-1)/2),1);
    glm_compute=1;
end

% calculating FC
n_person = length(RSts_aal2_PKU);
for person = 1:n_person
    X = corr(RSts_aal2_PKU{person}(:,1:n_ROI));
    xtemp = [];
    for edge = 1:n_ROI
        xtemp = [xtemp,X(edge,edge+1:n_ROI)];
    end
    FC(person,:) = xtemp;
end
FC = 0.5*log((1+FC)./(1-FC));

fprintf('finish preparation.\n')
regresscov = [cov_PKU.age,cov_PKU.age.^2,cov_PKU.age.^3]; % covariates to be regressed out from FC
fprintf('starting regressing out covariates\n')

for i = 1:(n_ROI*(n_ROI-1)/2)
    if glm_compute==1
        glmstruct{i} = fitglm(regresscov,FC(:,i));
    end
    b = glmstruct{i}.Coefficients.Estimate;
    FC_regressed(:,i) = FC(:,i)-[ones(n_person,1),regresscov]*b;
end
if glm_compute==1
    save('glmstruct.mat','glmstruct');
end

% classification on test set
[~,prob] = predict(svmstruct,FC_regressed); % calculating the classification score
fprintf('finish calculating\n');
%
BSC = normcdf(prob(:,2)); % calculating the probability i.e. the BSC
output = double(BSC>0.5); % binarize the BSC to calculate the classification accuracy