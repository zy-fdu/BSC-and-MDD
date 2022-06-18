function display_0w_8w_mean(data_0w,data_8w,cluster,outdir,ylabel_str)

axis_min=floor(min([data_0w;data_8w]));
axis_max=ceil(max([data_0w;data_8w]));

cluster_num=unique(cluster);
data_0w_mean=zeros(size(cluster_num));
data_8w_mean=zeros(size(cluster_num));
data_0w_std=zeros(size(cluster_num));
data_8w_std=zeros(size(cluster_num));
for i=1:length(cluster_num)
    data_0w_mean(i,1)=mean(data_0w(find(cluster==cluster_num(i))));
    data_8w_mean(i,1)=mean(data_8w(find(cluster==cluster_num(i))));
    data_0w_std(i,1)=std(data_0w(find(cluster==cluster_num(i))));
    data_8w_std(i,1)=std(data_8w(find(cluster==cluster_num(i))));
end

h=figure();
set(h,'DefaultAxesColorOrder',[ 0 0 1;1 0 0; 0 1 0]); 
parallelcoords([data_0w_mean,data_8w_mean], 'LineStyle', '-','LineWidth',2, ...
'Marker', '.', 'MarkerSize', 15,'group',[1;2;3]);
axis([0.5,2.5,axis_min,axis_max]);
hold on;

%gscatter(delta_gc_regressed,delta_HAMD_regressed,group);
if strcmp(ylabel_str,'HAMD')
    x1=ones(size(data_0w))+randi(size(data_0w,1),size(data_0w))/1000;
    x2=2*ones(size(data_8w))+randi(size(data_0w,1),size(data_0w))/1000;
%     axis_min=0;
%     axis_max=35;
%     axis([0.5,2.5,axis_min,axis_max]);
else
    x1=ones(size(data_0w));
    x2=2*ones(size(data_8w));
    
end

gscatter(x1,data_0w,cluster,[ 0 0 1;1 0 0; 0 1 0],'.',10);
gscatter(x2,data_8w,cluster,[ 0 0 1;1 0 0; 0 1 0],'.',10);

err=[data_0w_std,data_8w_std]./sqrt(size(data_0w_mean,1));
x=[1,1,1;2,2,2]';
y=[data_0w_mean';data_8w_mean']';
b=errorbar(x(1,:),y(1,:),err(1,:),'linestyle','None','Marker','*','LineWidth',1);
errorbar(x(2,:),y(2,:),err(2,:),'linestyle','None','Marker','*','LineWidth',1);
errorbar(x(3,:),y(3,:),err(3,:),'linestyle','None','Marker','*','LineWidth',1);

set(gca,'Xticklabel',{'baseline','8 week'});
ylabel(ylabel_str,'Interpreter', 'none');  
xlabel('Time')
legend({['gc<0.35, nSubj=',num2str(numel(find(cluster==1)))],['0.35<=gc<=0.65, nSubj=',num2str(numel(find(cluster==2)))],['gc>0.65, nSubj=',num2str(numel(find(cluster==3)))]});
saveas(gcf,[outdir,ylabel_str,'_0w_8w_',date],'png')
saveas(gcf,[outdir,ylabel_str,'_0w_8w_',date],'fig')

% h=figure()
% set(h,'DefaultAxesColorOrder',[1 0 0; 0 0 1; 0 1 0]); 
% coordLineStyle = 'k.';
% group=cell(size(cluster));
% group(cluster==1)={'cluster1'};
% group(cluster==1)={'cluster2'};
% group(cluster==1)={'cluster3'};
% %boxplot([data_0w,data_8w], 'Symbol', coordLineStyle,'colors','k'); hold on;
% parallelcoords([data_0w,data_8w], 'LineStyle', '-','LineWidth',2, ...
% 'Marker', '.', 'MarkerSize', 15,'group',cluster);
% axis([0.5,2.5,axis_min,axis_max]);
% hold on;
% 
% title(ylabel_str,'Interpreter', 'none');
% set(gca,'Xticklabel',{'baseline','8 week'});
% set(gca,'FontSize',16);
% ylabel(ylabel_str,'Interpreter', 'none');  
% xlabel('Time')

end



