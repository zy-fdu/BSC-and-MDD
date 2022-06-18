function display_deltaHAMD3groups(delta,BSCgroup,outdir)

axis_min=floor(min([delta]));
axis_max=ceil(max([delta]));

group_num=unique(BSCgroup);

    delta_mean(1)=mean(delta(find(BSCgroup==1)));
    delta_std(1)=std(delta(find(BSCgroup==1)));
    delta_mean(2)=mean(delta(find(BSCgroup==2)));
    delta_std(2)=std(delta(find(BSCgroup==2)));
    delta_mean(3)=mean(delta(find(BSCgroup==3)));
    delta_std(3)=std(delta(find(BSCgroup==3)));
   
    
    h=figure();
    set(h,'DefaultAxesColorOrder',[ 0 0 1;1 0 0; 0 1 0]);
    
    b=bar([1 2 3],delta_mean);
    b.FaceColor = 'w';
    b.LineWidth=1;
    b.EdgeColor = 'flat';
    
    b.CData(:,:)=[ 0 0 1;1 0 0; 0 1 0];
    
    ylabel('Reduction of HAMD','Interpreter', 'none');
    xlabel('BSC group')
    set(gca,'Xticklabel',{'female-like','androgynous','male-like'});
    set(gca,'LineWidth',0.8);
    hold on;
  %  BSCgroup_rand=zeros(size(BSCgroup));
    
    for i=1:numel(BSCgroup)
        num=(numel(find(BSCgroup==BSCgroup(i)&delta==delta(i)))-1);
        if num>0
            if mod((numel(find(BSCgroup(1:i)==BSCgroup(i)&delta(1:i)==delta(i)))),2)==1
                num_curr=-floor((numel(find(BSCgroup(1:i)==BSCgroup(i)&delta(1:i)==delta(i))))/2);
            else
                num_curr=numel(find(BSCgroup(1:i)==BSCgroup(i)&delta(1:i)==delta(i)))/2;
            end
            BSCgroup_rand(i)=BSCgroup(i)+num_curr/10;
        else
            BSCgroup_rand(i)=BSCgroup(i);
        end
    end
    
    gscatter(BSCgroup_rand,delta,BSCgroup,[ 0 0 1;1 0 0; 0 1 0],'.',10);
    axis([0,4,-1,50]);
    set(gca,'yTick',[0:10:50]);
    legend('off');
    saveas(gcf,[outdir,'HAMD_reduction_threegroups_',date],'png')
    saveas(gcf,[outdir,'HAMD_reduction_threegroups_',date],'fig')

end

