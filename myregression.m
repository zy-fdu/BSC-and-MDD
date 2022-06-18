function output=myregression(source_data,cov)
n_person=size(cov,1);
glmstruct=fitglm(cov,source_data);
b=glmstruct.Coefficients.Estimate;
output=source_data-[ones(n_person,1),cov]*b;