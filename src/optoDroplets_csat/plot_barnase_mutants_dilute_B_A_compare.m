clear all;

% Plot Figure S2B

dir1={'Barnase_additional_mutants'};
myloc={'cytoplasm'};

da=load(['../../data/optoDroplets_csat/' dir1{1} '_' myloc{1} '_mean_csat.mat']);
mInt=da.mInt;
sInt=da.sInt;
diro=da.dir2;
dG=-da.dG;

bins=min(dG):0.5:max(dG)+0.5;
mycolor=cool(length(bins));
figure;
x=[0:1e3:4e4];
y=[0:1e3:4e4];
count=0;
for d=1:length(diro)
    da=load(['../../data/optoDroplets_csat/' dir1{1} '_' myloc{1} '_' diro{d} '_dilute_ints_outliers_removed.mat']);
    fullDiluteIntB=da.fullDiluteIntB;
    fullDiluteIntA=da.fullDiluteIntA;
    %return
    count=count+1;
    subplot(2,7,count)
    plot(x,y,'-k','linewidth',3); hold on; 
    
    
    subplot(2,7,count)
    [~,idx]=histc(dG(d),bins);
    s=scatter(fullDiluteIntB,fullDiluteIntA,60,'o','markerfacecolor',mycolor(idx,:),'markeredgecolor','k'); hold on;
    plot([0 mInt(d)],[mInt(d) mInt(d)],'--','color',mycolor(idx,:),'linewidth',1);
    title(diro{d});
    
    if isnan(mInt(d))==0
        pos=find(fullDiluteIntB>=mInt(d));
        [f,gof]=fit(fullDiluteIntB(pos)'-mInt(d),fullDiluteIntA(pos)','poly1','Lower',[0.1 mInt(d)],'Upper',[1 mInt(d)]);
        clear pos; 
        pos=find(x>=mInt(d));
        x2=x(pos);
        y2=f.p1.*(x2-f.p2)+f.p2;
        plot(x2,y2,'--','color',mycolor(idx,:),'linewidth',1);
        clear pos; clear x2; clear y2; clear f; 
    else
        pos=find((fullDiluteIntB-fullDiluteIntA)./fullDiluteIntB<0.10);
        m=max(fullDiluteIntA(pos));
        mincsat(d)=m;
        plot([0 m],[m m],':','color',mycolor(idx,:),'linewidth',1);
        if strcmp(diro{d},'WT')==1 
            for dl=1:length(diro)
                subplot(2,7,dl)
                plot([0 m],[m m],':','color',mycolor(idx,:),'linewidth',1); hold on; 
                m
            end
        end
        %return
    end
    
    xlim([0 4e4]);
    ylim([0 4e4]);
    clear fullDiluteIntA; clear fullDiluteIntB;
    %return
end

figure; 
colormap(mycolor)
cb=colorbar;
set(cb,'YTick',bins(1):2:bins(end))
caxis([bins(1) bins(end)])