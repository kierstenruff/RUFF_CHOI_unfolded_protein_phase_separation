clear all;

% Plot Figure S3

dir1={'Barnase_sticker_mutants'};
myloc={'cytoplasm'};

da=load(['../../data/optoDroplets_csat/' dir1{1} '_' myloc{1} '_mean_csat.mat']);
mInt=da.mInt;
sInt=da.sInt;
diro=da.dir2;


%mycolor=cool(length(mInt));
mycolor=[0.8980 0.5176 0.5216; 0.8980 0.5176 0.5216; 0.8980 0.5176 0.5216; 0.8980 0.5176 0.5216; 0.8 0.8 0.8; 0.7059 0.8510 0.6235; 0.7059 0.8510 0.6235; 0.7059 0.8510 0.6235; 0.7059 0.8510 0.6235; 0.7059 0.8510 0.6235; 0.2980 0.4431 0.8000];
figure;
x=[0:1e3:5e4];
y=[0:1e3:5e4];
count=0;
for d=1:length(diro)
    da=load(['../../data/optoDroplets_csat/' dir1{1} '_' myloc{1} '_' diro{d} '_dilute_ints_outliers_removed.mat']);
    fullDiluteIntB=da.fullDiluteIntB;
    fullDiluteIntA=da.fullDiluteIntA;
    
    count=count+1;
    subplot(2,6,count)
    plot(x,y,'-k','linewidth',3); hold on; 
 
    
    subplot(2,6,count)
    s=scatter(fullDiluteIntB,fullDiluteIntA,60,'o','markerfacecolor',mycolor(d,:),'markeredgecolor','k'); hold on;
    plot([0 mInt(d)],[mInt(d) mInt(d)],'--','color',mycolor(d,:),'linewidth',1);
    title(diro{d});
    
    if isnan(mInt(d))==0
        pos=find(fullDiluteIntB>=mInt(d));
        [f,gof]=fit(fullDiluteIntB(pos)'-mInt(d),fullDiluteIntA(pos)','poly1','Lower',[0.1 mInt(d)],'Upper',[1 mInt(d)]);
        clear pos; 
        pos=find(x>=mInt(d));
        x2=x(pos);
        y2=f.p1.*(x2-f.p2)+f.p2;
        plot(x2,y2,'--','color',mycolor(d,:),'linewidth',1);
        clear pos; clear x2; clear y2; clear f; 
    else
        pos=find((fullDiluteIntB-fullDiluteIntA)./fullDiluteIntB<0.10);
        m=max(fullDiluteIntA(pos));
        mincsat(d)=m;
        plot([0 m],[m m],':','color',mycolor(d,:),'linewidth',1);
        %return
    end
    
    
    xlim([0 5e4]);
    ylim([0 5e4]);
    
    
    if d==5
        for i=[1 2 3 4 6 7 8 9 10 11]
            subplot(2,6,i)
            plot([0 mInt(d)],[mInt(d) mInt(d)],'--','color',mycolor(d,:),'linewidth',1); hold on; 
        end
    end
    
    clear fullDiluteIntA; clear fullDiluteIntB;
    %return
end

%figure; 
%colormap(mycolor)
%cb=colorbar;
%set(cb,'YTick',1:2:length(bins))
%caxis([bins(1) bins(end)])