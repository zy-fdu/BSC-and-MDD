function display_parallel_HAMD(cov_PKU_0w_8w)
%% para plot HAMD three groups & mean value
HAMD_0w_mean=[mean(cov_PKU_0w_8w.HAMD_0w(index_35));mean(cov_PKU_0w_8w.HAMD_0w(index_35_65));mean(cov_PKU_0w_8w.HAMD_0w(index_65))];
HAMD_8w_mean=[mean(cov_PKU_0w_8w.HAMD_8w(index_35));mean(cov_PKU_0w_8w.HAMD_8w(index_35_65));mean(cov_PKU_0w_8w.HAMD_8w(index_65))];
HAMD_0w_std=[std(cov_PKU_0w_8w.HAMD_0w(index_35));std(cov_PKU_0w_8w.HAMD_0w(index_35_65));std(cov_PKU_0w_8w.HAMD_0w(index_65))];
HAMD_8w_std=[std(cov_PKU_0w_8w.HAMD_8w(index_35));std(cov_PKU_0w_8w.HAMD_8w(index_35_65));std(cov_PKU_0w_8w.HAMD_8w(index_65))];

figure();
parallelcoords([HAMD_0w_mean,HAMD_8w_mean], 'LineStyle', '-',...
'Marker', '.', 'MarkerSize', 10,'group',[1;2;3]);
%title('androgynous','Interpreter', 'none');
axis([0.5,2.5,0,30]);
hold on;
err=[HAMD_0w_std,HAMD_8w_std]./sqrt(size(HAMD_0w_mean,1));
x=[1,1,1;2,2,2]';
y=[HAMD_0w_mean';HAMD_8w_mean']';
b=errorbar(x(1,:),y(1,:),err(1,:),'linestyle','None','Color',[0,0.45,0.74])
errorbar(x(2,:),y(2,:),err(2,:),'linestyle','None','Color',[0.85,0.33,0.10])
errorbar(x(3,:),y(3,:),err(3,:),'linestyle','None','Color',[0.93,0.69,0.13])

set(gca,'Xticklabel',{'baseline','8 week'});
ylabel('HAMD','Interpreter', 'none');  
xlabel('Group')
legend({'female-like','androgynous','male-like'});

saveas(gcf,[outdir,'/HAMD_0w_8w_mean'],'png');
saveas(gcf,[outdir,'/HAMD_0w_8w_mean'],'fig');

end