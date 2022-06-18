function display_boxplot_parallel_gc(cov_dep_0w_8w)

index_35=numel(find(cov_dep_0w_8w.groupBinary==0))+find(cov_dep_0w_8w.gc_0w(cov_dep_0w_8w.groupBinary)<0.35);
index_35_65=numel(find(cov_dep_0w_8w.groupBinary==0))+find(cov_dep_0w_8w.gc_0w(cov_dep_0w_8w.groupBinary)>=0.35&cov_dep_0w_8w.gc_0w(cov_dep_0w_8w.groupBinary)<=0.65);
index_65=numel(find(cov_dep_0w_8w.groupBinary==0))+find(cov_dep_0w_8w.gc_0w(cov_dep_0w_8w.groupBinary)>0.65);

%% plot paired T  0w 8w, three groups 0.35-0.65
figure()
coordLineStyle = 'k.';
labels={'baseline','8 week'};
group=cell(size(index_35,1),1);
group(cov_dep_0w_8w.sex(index_35))={'male'};
group(~cov_dep_0w_8w.sex(index_35))={'female'};
figure();
boxplot([cov_dep_0w_8w.gc_0w_regressed(index_35),cov_dep_0w_8w.gc_8w_regressed(index_35)], 'Symbol', coordLineStyle,'colors','k'); hold on;
parallelcoords([cov_dep_0w_8w.gc_0w_regressed(index_35),cov_dep_0w_8w.gc_8w_regressed(index_35)], 'LineStyle', '-',...
'Marker', '.', 'MarkerSize', 10,'Labels',labels,'group',group);
title('gc<0.35_with regression','Interpreter', 'none');
set(gca,'Xticklabel',{'baseline','8 week'});
ylabel('SVM_gc','Interpreter', 'none');  
saveas(gcf, [outdir,'\gc_0w_8w_threegroups_group1_with regression'], 'png');
saveas(gcf, [outdir,'\gc_0w_8w_threegroups_group1_with regression'], 'fig');


group=cell(size(index_35_65,1),1);
group(cov_dep_0w_8w.sex(index_35_65))={'male'};
group(~cov_dep_0w_8w.sex(index_35_65))={'female'};
figure();
boxplot([cov_dep_0w_8w.gc_0w_regressed(index_35_65),cov_dep_0w_8w.gc_8w_regressed(index_35_65)], 'Symbol', coordLineStyle,'colors','k'); hold on;
parallelcoords([cov_dep_0w_8w.gc_0w_regressed(index_35_65),cov_dep_0w_8w.gc_8w_regressed(index_35_65)], 'LineStyle', '-',...
'Marker', '.', 'MarkerSize', 10,'Labels',labels,'group',group);
title('0.35<gc<0.65_with regression','Interpreter', 'none');
set(gca,'Xticklabel',{'baseline','8 week'});
ylabel('SVM_gc','Interpreter', 'none');  
saveas(gcf, [outdir,'\gc_0w_8w_threegroups_group2_with regression'], 'png');
saveas(gcf, [outdir,'\gc_0w_8w_threegroups_group2_with regression'], 'fig');


group=cell(size(index_65,1),1);
group(cov_dep_0w_8w.sex(index_65))={'male'};
group(~cov_dep_0w_8w.sex(index_65))={'female'};
figure();
boxplot([cov_dep_0w_8w.gc_0w_regressed(index_65),cov_dep_0w_8w.gc_8w_regressed(index_65)], 'Symbol', coordLineStyle,'colors','k'); hold on;
parallelcoords([cov_dep_0w_8w.gc_0w_regressed(index_65),cov_dep_0w_8w.gc_8w_regressed(index_65)], 'LineStyle', '-',...
'Marker', '.', 'MarkerSize', 10,'Labels',labels,'group',group);
title('gc>0.65_with regression','Interpreter', 'none');
set(gca,'Xticklabel',{'baseline','8 week'});
ylabel('SVM_gc','Interpreter', 'none');  
saveas(gcf, [outdir,'\gc_0w_8w_threegroups_group3_with regression'], 'png');
saveas(gcf, [outdir,'\gc_0w_8w_threegroups_group3_with regression'], 'fig');

end