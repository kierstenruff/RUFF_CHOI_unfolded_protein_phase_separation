clear all; 

% Plot Figure 5D

da=load(['../../data/optoDroplets_csat/Barnase_hsp_coexpression_cytoplasm_mean_csat.mat']);
mInt=da.mInt;
sInt=da.sInt;
diro=da.dir2;

da2=load(['../../data/optoDroplets_csat/Barnase_hsp_coexpression_nucleus_mean_csat.mat']);
mInt2=da2.mInt;
sInt2=da2.sInt;
diro2=da2.dir2;

csat(:,1)=mInt2([5 7]);
csat(:,2)=mInt([5 7]);
csaterr(:,1)=sInt2([5 7]);
csaterr(:,2)=sInt([5 7]);

figure;
bar(csat); hold on; 
set(gca,'XTickLabel',{'L14A','L14A + Hsp70'});
legend('nucleus','cytoplasm');
ylabel('Barnase Concentration (AU)');

errorbar([1 1],csat(1,:),csaterr(1,:),'o'); hold on; 
errorbar([2 2],csat(2,:),csaterr(2,:),'o'); hold on; 

clear da; 
dir1={'Barnase_hsp_coexpression'};
myloc={'nucleus'};

mycolor=cool(4);
figure;
x=[0:1e3:4e4];
y=[0:1e3:4e4];
for d=5:8%1:length(dir2)
    da=load(['../../data/optoDroplets_csat/' dir1{1} '_' myloc{1} '_' diro2{d} '_dilute_ints_outliers_removed.mat']);
    fullDiluteIntB=da.fullDiluteIntB;
    fullDiluteIntA=da.fullDiluteIntA;
    
    subplot(2,2,d-4)
    plot(x,y,'-k','linewidth',3); hold on; 
    
    if d==5
        for i=2:4
            subplot(2,2,i)
            s=scatter(fullDiluteIntB,fullDiluteIntA,60,'o','markerfacecolor',mycolor(d-4,:),'markeredgecolor','k'); hold on;
            alpha(s,0.5)
            plot([0 mInt2(d)],[mInt2(d) mInt2(d)],'--','color',mycolor(d-4,:),'linewidth',1);
            
            pos=find(fullDiluteIntB>=mInt2(d));
            [f,gof]=fit(fullDiluteIntB(pos)'-mInt2(d),fullDiluteIntA(pos)','poly1','Lower',[0.1 mInt2(d)],'Upper',[1 mInt2(d)]);
            clear pos; 
            pos=find(x>=mInt2(d));
            x2=x(pos);
            y2=f.p1.*(x2-f.p2)+f.p2;
            plot(x2,y2,'--','color',mycolor(d-4,:),'linewidth',1);
            clear pos; clear x2; clear y2; clear f; 
        end
    end
    
    subplot(2,2,d-4)
    s=scatter(fullDiluteIntB,fullDiluteIntA,60,'o','markerfacecolor',mycolor(d-4,:),'markeredgecolor','k'); hold on;
    plot([0 mInt2(d)],[mInt2(d) mInt2(d)],'--','color',mycolor(d-4,:),'linewidth',1);
    
    if isnan(mInt2(d))==0
        pos=find(fullDiluteIntB>=mInt2(d));
        [f,gof]=fit(fullDiluteIntB(pos)'-mInt2(d),fullDiluteIntA(pos)','poly1','Lower',[0.1 mInt2(d)],'Upper',[1 mInt2(d)]);
        clear pos; 
        pos=find(x>=mInt2(d));
        x2=x(pos);
        y2=f.p1.*(x2-f.p2)+f.p2;
        plot(x2,y2,'--','color',mycolor(d-4,:),'linewidth',1);
        clear pos; clear x2; clear y2; clear f; 
    else
        pos=find((fullDiluteIntB-fullDiluteIntA)./fullDiluteIntB<0.10);
        m=max(fullDiluteIntA(pos));
        plot([0 m],[m m],'-.','color',mycolor(d-4,:),'linewidth',1);
        if d==6
            csat(3,1)=m;
            csaterr(3,1)=0;
        elseif d==8
            csat(4,1)=m;
            csaterr(4,1)=0;
        end
        %return
        clear pos; clear m; 
    end
    
    xlim([0 4e4]);
    ylim([0 4e4]);
    clear fullDiluteIntA; clear fullDiluteIntB; clear da; 
    %return
end

dir1={'Barnase_hsp_coexpression'};
myloc={'cytoplasm'};

figure;
x=[0:1e3:4e4];
y=[0:1e3:4e4];
for d=5:8%1:length(dir2)
    da=load(['../../data/optoDroplets_csat/' dir1{1} '_' myloc{1} '_' diro{d} '_dilute_ints_outliers_removed.mat']);
    fullDiluteIntB=da.fullDiluteIntB;
    fullDiluteIntA=da.fullDiluteIntA;
    
    subplot(2,2,d-4)
    plot(x,y,'-k','linewidth',3); hold on; 
    
    if d==5
        for i=2:4
            subplot(2,2,i)
            s=scatter(fullDiluteIntB,fullDiluteIntA,60,'o','markerfacecolor',mycolor(d-4,:),'markeredgecolor','k'); hold on;
            alpha(s,0.5)
            plot([0 mInt(d)],[mInt(d) mInt(d)],'--','color',mycolor(d-4,:),'linewidth',1);
            
            pos=find(fullDiluteIntB>=mInt(d));
            [f,gof]=fit(fullDiluteIntB(pos)'-mInt(d),fullDiluteIntA(pos)','poly1','Lower',[0.1 mInt(d)],'Upper',[1 mInt(d)]);
            clear pos; 
            pos=find(x>=mInt(d));
            x2=x(pos);
            y2=f.p1.*(x2-f.p2)+f.p2;
            plot(x2,y2,'--','color',mycolor(d-4,:),'linewidth',1);
            clear pos; clear x2; clear y2; clear f; 
        end
    end
    
    subplot(2,2,d-4)
    s=scatter(fullDiluteIntB,fullDiluteIntA,60,'o','markerfacecolor',mycolor(d-4,:),'markeredgecolor','k'); hold on;
    plot([0 mInt(d)],[mInt(d) mInt(d)],'--','color',mycolor(d-4,:),'linewidth',1);
    
    if isnan(mInt(d))==0
        pos=find(fullDiluteIntB>=mInt(d));
        [f,gof]=fit(fullDiluteIntB(pos)'-mInt(d),fullDiluteIntA(pos)','poly1','Lower',[0.1 mInt(d)],'Upper',[1 mInt(d)]);
        clear pos; 
        pos=find(x>=mInt(d));
        x2=x(pos);
        y2=f.p1.*(x2-f.p2)+f.p2;
        plot(x2,y2,'--','color',mycolor(d-4,:),'linewidth',1);
        clear pos; clear x2; clear y2; clear f; 
    else

        pos=find((fullDiluteIntB-fullDiluteIntA)./fullDiluteIntB<0.10);
        m=max(fullDiluteIntA(pos));
        plot([0 m],[m m],'-.','color',mycolor(d-4,:),'linewidth',1);
        if d==6
            csat(3,2)=m;
            csaterr(3,2)=0;
        elseif d==8
            csat(4,2)=m;
            csaterr(4,2)=0;
        end
        %return
    end
    
    xlim([0 4e4]);
    ylim([0 4e4]);
    clear fullDiluteIntA; clear fullDiluteIntB; clear da; 
    %return
end

figure; 
barh(csat); hold on; 
set(gca,'YDir','Reverse')
legend('nucleus','cytoplasm');
set(gca,'yticklabel',{'L14A','L14A + Hsp70','L14A + Hsp40','L14A + Hsp40 + Hsp70'});

errorbar(csat(1,1:2),[1 1],[],[],csaterr(1,1:2),csaterr(1,1:2),'o'); hold on; 
errorbar(csat(2,1:2),[2 2],[],[],csaterr(2,1:2),csaterr(2,1:2),'o'); hold on;


f=figure;
f.Position=[100 100 600 200];
bar(csat/1000); hold on; 
legend('nucleus','cytoplasm');
set(gca,'xticklabel',{'L14A','L14A + Hsp70','L14A + Hsp40','L14A + Hsp40 + Hsp70'});
ylim([0 25])
set(gca,'YTick',[5:10:25]);
ylabel('csat (a.u.)');

errorbar([1 1],csat(1,1:2)/1000,csaterr(1,1:2)/1000,'o'); hold on; 
errorbar([2 2],csat(2,1:2)/1000,csaterr(2,1:2)/1000,'o'); hold on;

%X=[csat(:,1)' csat(:,2)'; csaterr(:,1)' csaterr(:,2)']
%xlswrite('csat_hsp_data.xls',X)
%return