clear all; 

% Plot Figures 4D and 4E

da=load(['../../data/optoDroplets_csat/Barnase_sticker_mutants_cytoplasm_mean_csat.mat']);
mInt=da.mInt;
sInt=da.sInt;
diro=da.dir2;

mInt(1:4)=[4.8581e4 4.4742e4 3.6323e4 4.3997e4];

%X=[mInt; sInt]
%xlswrite('csat_sticker_mutants_data.xls',X)
%return

figure;
for i=1:length(mInt)
    if i<5
        mycolor=[0.8980 0.5176 0.5216];
        h=barh(i, mInt(i)-mInt(5)); hold on; 
        %h=barh(i, mInt(i)/mInt(5)); hold on; 
        set(h, 'FaceColor', mycolor) 
    elseif i>5 & i<=length(mInt)-1
        mycolor=[0.7059 0.8510 0.6235];
        h=barh(i, mInt(i)-mInt(5)); hold on;
        %h=barh(i, mInt(i)/mInt(5)); hold on; 
        set(h, 'FaceColor', mycolor) 
    elseif i==length(mInt)
        mycolor=[0.2980 0.4431 0.8000];
        h=barh(i+1, mInt(i)-mInt(5));
        %h=barh(i, mInt(i)/mInt(5)); hold on; 
        set(h, 'FaceColor', mycolor) 
    end
end
%xlim([-12000 12000])
%plot([sInt(5) sInt(5)],[0 length(mInt)+3],'-k'); hold on; 
%plot([2000 2000],[0 length(mInt)+3],'-k'); hold on;
%plot([-2000 -2000],[0 length(mInt)+3],'-k'); hold on;
set(gca,'YTick',[1 2 3 4 6 7 8 9 10 12]);
set(gca,'YTickLabel',{'8xA:9xS','8xA:7xS','8xA:5xS','8xA:3xS','8xA:F7Y','8xA:F106Y','8xA:F7Y, F106Y','8xA:F7Y,F82Y,F106Y','8xA:F7Y,F56Y,F82Y,F106Y','8xA:Q2Y, G40Y, G65Y, T100Y'});  
set(gca,'YDir','Reverse');


A50(1,:)=[3.489 3.498 3.486]; %L14A
A50(2,:)=[3.44 3.37 3.375]; %I25A, I96G
A50(3,:)=[3.531 3.535 3.505]; %8xAla
A50(4,:)=[3.647 3.566 3.528]; %8xAla + F7Y, F56Y, F82Y, F106Y
A50(5,:)=[3.213 3.192 3.187]; %8xAla + Q2Y, G40Y, G65Y, T100Y

mA50(1,:)=mean(A50,2);
sA50(1,:)=std(A50,0,2);

da2=load(['../../data/optoDroplets_csat/Barnase_additional_mutants_cytoplasm_mean_csat.mat']);
mInt2=da2.mInt;
sInt2=da2.sInt;
diro2=da2.dir2;

mcsat(1:2)=[mInt2(6) mInt2(13)];
mcsat(3:5)=[mInt(5) mInt(10) mInt(11)];
scsat(1:2)=[sInt2(6) sInt2(13)];
scsat(3:5)=[sInt(5) sInt(10) sInt(11)];

%mycolor=[0.2980 0.4431 0.8000; 0.7059 0.8510 0.6235; 0.8980 0.5176 0.5216; 0.82 0.82 0.82; 0 0 0]
mycolor=[0.8000 0.2980 0.3216; 0.4078 0.2667 0.5765; 0 0 0; 0.2235 0.4863 0.2353; 0.1922 0.2392 0.5961]
figure; 

mcsat=mcsat/1000;
scsat=scsat/1000;
[fitvals,s]=polyfit(mcsat',mA50',1);
%x=min(mcsat)-0.5:0.1:max(mcsat)+0.5;
x=6.5:0.1:19.5;
y=fitvals(1)*x+fitvals(2);
plot(x,y,'-k','linewidth',4); hold on; 
[yfit, dy]=polyconf(fitvals,x,s,'predopt','curve');
plot(x,yfit - dy,'-k');
plot(x,yfit + dy,'-k');
patch([x, fliplr(x)], [yfit - dy fliplr(yfit + dy)], [211 211 211]/255, 'EdgeColor','none', 'FaceAlpha',0.5)

%plot(mcsat,mA50,'o','color',mycolor,'markersize',10); hold on; 
%errorbar(mcsat,mA50,sA50,sA50,scsat,scsat,'o','color',mycolor,'markersize',10); hold on; 
for i=1:length(mcsat)
    errorbar(mcsat(i),mA50(i),sA50(i),sA50(i),scsat(i),scsat(i),'o','markerfacecolor',mycolor(i,:),'color',mycolor(i,:),'markersize',10); hold on;
end
legend('L14A','I25A, I96G','8xA','8xA:F7Y, F56Y, F82Y, F106Y','8xA:Q2Y, G40Y, G65Y, T100Y');
    
[f gof]=fit(mcsat',mA50','poly1')
