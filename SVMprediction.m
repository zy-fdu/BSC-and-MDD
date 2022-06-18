function [gc,output]=SVMprediction(dep_TC,cov_dep,svmstruct,glmstruct_path)
%%
n_ROI = 94;

if strcmp(class(glmstruct_path),'char')
    glmstruct=importdata(glmstruct_path);
    glm_compute=0;
else
    glmstruct=cell((n_ROI*(n_ROI-1)/2),1);
    glm_compute=1;
end


n_person = length(dep_TC);
for person = 1:n_person
    X = corr(dep_TC{person}(:,1:n_ROI));
    xtemp = [];
    for edge = 1:n_ROI
        xtemp = [xtemp,X(edge,edge+1:n_ROI)];
    end
    FC(person,:) = xtemp;
end
FC = 0.5*log((1+FC)./(1-FC));

fprintf('finish preparation.\n')
regresscov = [cov_dep.age,cov_dep.age.^2,cov_dep.age.^3];
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

% predicting on test set

[~,prob] = predict(svmstruct,FC_regressed);

    
fprintf('finish calculating\n');

%
gc = normcdf(prob(:,2));
output = double(gc>0.5);