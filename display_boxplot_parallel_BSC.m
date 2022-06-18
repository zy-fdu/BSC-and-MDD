function display_boxplot_parallel_BSC(cov_dep_0w_8w)

index_35=numel(find(cov_dep_0w_8w.groupBinary==0))+find(cov_dep_0w_8w.BSC_0w(cov_dep_0w_8w.groupBinary)<0.35);
index_35_65=numel(find(cov_dep_0w_8w.groupBinary==0))+find(cov_dep_0w_8w.BSC_0w(cov_dep_0w_8w.groupBinary)>=0.35&cov_dep_0w_8w.BSC_0w(cov_dep_0w_8w.groupBinary)<=0.65);
index_65=numel(find(cov_dep_0w_8w.groupBinary==0))+find(cov_dep_0w_8w.BSC_0w(cov_dep_0w_8w.groupBinary)>0.65);

%% plot paired T  0w 8w, three groups 0.35-0.65
figure()
coordLineStyle = 'k.';
labels={'baseline','8 week'};
group=cell(size(index_35,1),1);
group(cov_dep_0w_8w.sex(index_35))={'male'};
group(~cov_dep_0w_8w.sex(index_35))={'female'};
figure();
boxplot([cov_dep_0w_8w.BSC_0w_regressed(index_35),cov_dep_0w_8w.BSC_8w_regressed(index_35)], 'Symbol', coordLineStyle,'colors','k'); hold on;
parallelcoords([cov_dep_0w_8w.BSC_0w_regressed(index_35),cov_dep_0w_8w.BSC_8w_regressed(index_35)], 'LineStyle', '-',...
'Marker', '.', 'MarkerSize', 10,'Labels',labels,'group',group);
title('female-like','Interpreter', 'none');
set(gca,'Xticklabel',{'baseline','8 week'});
ylabel('SVM_BSC','Interpreter', 'none');  
saveas(gcf, [outdir,'\BSC_0w_8w_femalelike'], 'png');
saveas(gcf, [outdir,'\BSC_0w_8w_femalelike'], 'fig');


group=cell(size(index_35_65,1),1);
group(cov_dep_0w_8w.sex(index_35_65))={'male'};
group(~cov_dep_0w_8w.sex(index_35_65))={'female'};
figure();
boxplot([cov_dep_0w_8w.BSC_0w_regressed(index_35_65),cov_dep_0w_8w.BSC_8w_regressed(index_35_65)], 'Symbol', coordLineStyle,'colors','k'); hold on;
parallelcoords([cov_dep_0w_8w.BSC_0w_regressed(index_35_65),cov_dep_0w_8w.BSC_8w_regressed(index_35_65)], 'LineStyle', '-',...
'Marker', '.', 'MarkerSize', 10,'Labels',labels,'group',group);
title('androgynous','Interpreter', 'none');
set(gca,'Xticklabel',{'baseline','8 week'});
ylabel('BSC','Interpreter', 'none');  
saveas(gcf, [outdir,'\BSC_0w_8w_androgynous'], 'png');
saveas(gcf, [outdir,'\BSC_0w_8w_androgynous'], 'fig');


group=cell(size(index_65,1),1);
group(cov_dep_0w_8w.sex(index_65))={'male'};
group(~cov_dep_0w_8w.sex(index_65))={'female'};
figure();
boxplot([cov_dep_0w_8w.BSC_0w_regressed(index_65),cov_dep_0w_8w.BSC_8w_regressed(index_65)], 'Symbol', coordLineStyle,'colors','k'); hold on;
parallelcoords([cov_dep_0w_8w.BSC_0w_regressed(index_65),cov_dep_0w_8w.BSC_8w_regressed(index_65)], 'LineStyle', '-',...
'Marker', '.', 'MarkerSize', 10,'Labels',labels,'group',group);
title('male-like','Interpreter', 'none');
set(gca,'Xticklabel',{'baseline','8 week'});
ylabel('SVM_BSC','Interpreter', 'none');  
saveas(gcf, [outdir,'\BSC_0w_8w_malelike'], 'png');
saveas(gcf, [outdir,'\BSC_0w_8w_malelike'], 'fig');

end