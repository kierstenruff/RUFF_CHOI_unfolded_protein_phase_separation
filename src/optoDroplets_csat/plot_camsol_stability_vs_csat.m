clear all;

% Plot Figure 3E

constructs={'WT','I55G','V45T','I25A','I88A','L14A','I51A','L89G','I96G','I55G, L89G','I88G','V45T, I88G','I25A, I96G','L14A,I51A,I88A','L14G,I51G,I88G','L14S,I51S,I88S','L14D,I51D,I88D','8xAla','8xSer','8xAsp','8A,7S','8A,Y24S,Y90S,Y97S','8A,F7Y','8A,F106Y','8A,F7Y,F106Y','8A,F7Y,F82Y,F106Y','8A,F7Y,F56Y,F82Y,F106Y','8A,Q2Y,G40Y,G65Y,T100Y'};
mytitle={'WT','I55G','V45T','I25A','I88A','L14A','I51A','L89G','I96G','I55G,L89G','I88G','V45T,I88G','I25A,I96G','3xA','3xG','3xS','3xD','8xA','8xS','8xD','7S','3S','F7Y','F106Y','2FY','3FY','4FY','4Y'};

da=importdata('../../data/optoDroplets_csat/barnase_variants_dG_mcsat_scsat_camsol.txt');
names=da.textdata(1:end,1);
dG=da.data(:,1);
csat=da.data(:,2);
csol=da.data(:,4);

R=8.131;
T=293;
dGoffset=1.2875e+04; 
cstarU=1.0830e+04;


csatlist=[];
dGlist=[];
csollist=[];
pUlist=[];
pUlistnooff=[];
predcsat=[];
for c=1:length(constructs)

    pos=find(strcmp(names,constructs{c})==1);
    csatlist=[csatlist csat(pos)];
    dGlist=[dGlist dG(pos)];
    csollist=[csollist csol(pos)];

    tmpdG=dG(pos)*1000;
    tmp=(cstarU/(exp(-(-tmpdG-dGoffset)/(R*T))/(1+exp(-(-tmpdG-dGoffset)/(R*T)))))/1000;
    predcsat=[predcsat tmp];

    pU=(exp(-(-tmpdG-dGoffset)/(R*T))/(1+exp(-(-tmpdG-dGoffset)/(R*T))));
    pUlist=[pUlist pU];
    pUnooff=(exp(-(-tmpdG)/(R*T))/(1+exp(-(-tmpdG)/(R*T))));
    pUlistnooff=[pUlistnooff pUnooff];

    clear tmp; clear tmpdG; clear pF; 
end

ncsollist=csollist/csollist(1);


% Which constructs do we want to fit over?
%pos=find(isnan(csatlist)==0); % All mutants
%pos=4:13; % Initial mutants
%pos=4:13; % Initial mutants
%pos=14:19; % Extreme mutants
%pos=14:18; % Extreme mutants without 8xS
%pos=4:19; % Non-sticker mutants
pos=4:18; % Non-sticker mutants & 8xS

% Make color scale for pU
bins=0.45:0.05:1.00;
mycolor=jet(length(bins));


%% Stability
f=figure;
%f.Position=[100 100 500 500];
%f.Position=[100 100 1200 300];
f.Position=[100 100 250 600];
%subplot(1,3,1)
subplot(3,1,1)
pos=4:18;
[fitvals,s]=polyfit(csatlist(pos),1./pUlist(pos),1);
x=min(csatlist(pos))-0.5:0.1:max(csatlist(pos))+0.5;
y=fitvals(1)*x+fitvals(2);
plot(x,y,'-k','linewidth',4); hold on; 
[yfit, dy]=polyconf(fitvals,x,s,'predopt','curve');
plot(x,yfit - dy,'-k');
plot(x,yfit + dy,'-k');
patch([x, fliplr(x)], [yfit - dy fliplr(yfit + dy)], [211 211 211]/255, 'EdgeColor','none', 'FaceAlpha',0.5)
pos=4:19;
for p=1:length(pos)
    [~,idx]=histc(pUlist(pos(p)),bins);
    %idx
    %pUlist(pos(p))
    plot(csatlist(pos(p)),1/pUlist(pos(p)),'o','markersize',8,'color',mycolor(idx,:),'markerfacecolor',mycolor(idx,:)); hold on;
    text(csatlist(pos(p)),1/pUlist(pos(p)),mytitle(pos(p)),'color','k');
end
xlim([min(x)-0.5 max(x)+0.5])
ylabel('1/pU');
xlabel('csat');

pos=4:18;
x=csatlist(pos);
y=1./pUlist(pos);
tmp=fitlm(x,y)
r2=tmp.Rsquared.Ordinary
text(min(x),max(y),['R^2=' num2str(round(r2,2))]);
%tmp2=corrcoef(x,y);
%text(min(x),max(y)-0.25,['r=' num2str(round(tmp2(1,2),2))]);

%% Solubility
%f=figure;
%f.Position=[100 100 500 500];
%subplot(1,3,2)
subplot(3,1,2)
pos=4:18;
[fitvals,s]=polyfit(csatlist(pos),ncsollist(pos),1);
x=min(csatlist(pos))-0.5:0.1:max(csatlist(pos))+0.5;
y=fitvals(1)*x+fitvals(2);
plot(x,y,'-k','linewidth',4); hold on; 
[yfit, dy]=polyconf(fitvals,x,s,'predopt','curve');
plot(x,yfit - dy,'-k');
plot(x,yfit + dy,'-k');
patch([x, fliplr(x)], [yfit - dy fliplr(yfit + dy)], [211 211 211]/255, 'EdgeColor','none', 'FaceAlpha',0.5)
pos=4:19;
for p=1:length(pos)
    [~,idx]=histc(pUlist(pos(p)),bins);
    plot(csatlist(pos(p)),ncsollist(pos(p)),'o','markersize',8,'color',mycolor(idx,:),'markerfacecolor',mycolor(idx,:)); hold on;
    text(csatlist(pos(p)),ncsollist(pos(p)),mytitle(pos(p)),'color','k');
end
xlim([min(x)-0.5 max(x)+0.5])
ylabel('CamSol/CamSol_{WT}');
xlabel('csat');

pos=4:18;
x=csatlist(pos);
y=ncsollist(pos);
tmp=fitlm(x,y)
r2=tmp.Rsquared.Ordinary
text(min(x),max(y),['R^2=' num2str(round(r2,2))]);
%tmp2=corrcoef(x,y);
%text(min(x),max(y)-0.05,['r=' num2str(round(tmp2(1,2),2))]);

%% Stability and Solubility
%f=figure;
%f.Position=[100 100 500 500];
%subplot(1,3,3);
subplot(3,1,3)
pos=4:18;
[fitvals,s]=polyfit(csatlist(pos),1./pUlist(pos)+(pUlist(pos).*ncsollist(pos)),1);
x=min(csatlist(pos))-0.5:0.1:max(csatlist(pos))+0.5;
y=fitvals(1)*x+fitvals(2);
plot(x,y,'-k','linewidth',4); hold on; 
[yfit, dy]=polyconf(fitvals,x,s,'predopt','curve');
plot(x,yfit - dy,'-k');
plot(x,yfit + dy,'-k');
patch([x, fliplr(x)], [yfit - dy fliplr(yfit + dy)], [211 211 211]/255, 'EdgeColor','none', 'FaceAlpha',0.5)
pos=4:19;
for p=1:length(pos)
    [~,idx]=histc(pUlist(pos(p)),bins);
    plot(csatlist(pos(p)),1/pUlist(pos(p))+(pUlist(pos(p))*ncsollist(pos(p))),'o','markersize',8,'color',mycolor(idx,:),'markerfacecolor',mycolor(idx,:)); hold on;
    text(csatlist(pos(p)),1/pUlist(pos(p))+(pUlist(pos(p))*ncsollist(pos(p))),mytitle(pos(p)),'color','k');
end
xlim([min(x)-0.5 max(x)+0.5])
ylabel('1/pU+pU(CamSol/CamSol_{WT})');
xlabel('csat');

pos=4:18;
x=csatlist(pos);
y=1./pUlist(pos)+(pUlist(pos).*ncsollist(pos));
tmp=fitlm(x,y)
r2=tmp.Rsquared.Ordinary
text(min(x),max(y),['R^2=' num2str(round(r2,2))]);
%tmp2=corrcoef(x,y);
%text(min(x),max(y)-0.25,['r=' num2str(round(tmp2(1,2),2))]);

