clear all;

% Plot Figures S2C and 3C

da=load(['../../data/optoDroplets_csat/Barnase_additional_mutants_cytoplasm_mean_csat.mat']);
dGo=da.dG;
mInt=da.mInt;
sInt=da.sInt;
diro=da.dir2;
pos=find(strcmp(diro,'Ex2')==1);
diro{pos}='V45T, I88G';
pos=find(strcmp(diro,'Ex3')==1);
diro{pos}='I55G, L89G';
pos=find(strcmp(diro,'Ex4')==1);
diro{pos}='I25A, I96G';
%return

%X=[dGo; mInt; sInt]
%xlswrite('csat_initial_mutants_data.xls',X)
%return

%% Plot dGU vs csat and fit based on unfolded concentration

dG=-dGo*1000;
f=figure(142);
f.Position = [100 100 540 700];
errorbar(mInt,dG/1000,sInt,'horizontal','o','markersize',8,'color',[0.1882 0.2275 0.60],'markerfacecolor',[0.1882 0.2275 0.60]); hold on;
for i=1:length(mInt)
    text(mInt(i)+1000,(dG(i)/1000)+0.5,diro{(i)},'color',[0.1882 0.2275 0.60]);
end


R=8.131;
T=293;
xdata=dG(4:end);
ydata=mInt(4:end); 
fun=@(x,xdata)(x(1)./(exp(-(xdata+x(2))./(R*T))./(1+exp(-(xdata+x(2))./(R*T)))));
x0=[18000 -15000];
x=lsqcurvefit(fun,x0,xdata,ydata);

xfull=[-80000:100:20000];
y=(x(1)./(exp(-(xfull+x(2))./(R*T))./(1+exp(-(xfull+x(2))./(R*T)))));
figure(142)
plot(y,xfull./1000,'--','color',[0.1882 0.2275 0.60],'linewidth',2);
ylabel('\DeltaG_U');
xlabel('c_{sat} (AU)');

text(6000,22000/1000,['Offset: ' num2str(round(x(2)/1e3,2)) 'kJ/mol']);
text(6000,18000/1000,['c*_U: ' num2str(round(x(1)/1e3,2)) 'e3 (AU)']);
ylim([-80 25])
xlim([0.5e4 3e4])

dGoffset=x(2);
cstarU=x(1);
%return

%% Bootstrap analysis with probabilities based on if hydrophobic blobs are changing 
ntrials=1000; % set to 5 if need to by pass this
%ntrials=5;

probmat=load('../../data/optoDroplets_csat/Barnase_blob_change_pick_probability.mat');
seqnames=probmat.seqsnames;
problist=probmat.qblobpickprob;

count=0; 
for i=4:length(mInt)
    count=count+1; 
    currname=diro{i};
    pos=find(strcmp(seqnames,currname)==1);
    blobpickprob(count)=problist(pos);
    %return
end

%return
for n=1:ntrials
    csumprobs=cumsum(blobpickprob);
    myval=csumprobs(end)*rand(1,length(blobpickprob));
    for i=1:length(myval)
        pos(i)=find(myval(i)<csumprobs,1,'first');
    end
    %return
    xdataloo=xdata(pos);
    ydataloo=ydata(pos);
    fun=@(x,xdataloo)(x(1)./(exp(-(xdataloo+x(2))./(R*T))./(1+exp(-(xdataloo+x(2))./(R*T)))));
    x0=[18000 -15000];
    x=lsqcurvefit(fun,x0,xdataloo,ydataloo);
    allfit(n,:)=x;
    
    clear x; clear xdataloo; clear ydataloo; 
end   

% Actually get mean and std of bootstrap fit and use this to define range
% of potential values
x=mean(allfit);
tmp=std(allfit);
%y=(x(1)./(exp(-(xfull+x(2))./(R*T))./(1+exp(-(xfull+x(2))./(R*T)))));
%plot(y,xfull./1000,'--','color',[0 0.8 0.8],'linewidth',2);

x(1)=mean(allfit(:,1))+tmp(1);
x(2)=mean(allfit(:,2))-tmp(2);
maxcstarU=x(:,1);
maxdGoffset=x(2);
y=(x(1)./(exp(-(xfull+x(2))./(R*T))./(1+exp(-(xfull+x(2))./(R*T)))));
figure(142)
plot(y,xfull./1000,'-','color',[0.8 0.8 0.8],'linewidth',2);

%return

x(1)=mean(allfit(:,1))-tmp(1);
x(2)=mean(allfit(:,2))+tmp(2);
mincstarU=x(:,1);
mindGoffset=x(2);
y=(x(1)./(exp(-(xfull+x(2))./(R*T))./(1+exp(-(xfull+x(2))./(R*T)))));
figure(142)
plot(y,xfull./1000,'-','color',[0.8 0.8 0.8],'linewidth',2);

%return


%% Get stacked bars
% All of this was written as if dG was of folding not unfolding so must
% change to negative values
dG=-dG; 
dGoffset=-dGoffset;

[dG, idx]=sort(dG);
pF=exp(-(dG+dGoffset)./(R*T))./(1+exp(-(dG+dGoffset)./(R*T)));
pU=1-pF;

% Load in min possible csat for the constructs with no fitted mInt
% these will need to be a single bar with an
% arrow instead of a stacked bar
mInt(1:3)=[29777 24449 37026]; % values extracted from
%plot_barnase_mutants_dilute_B_A_compare.m where these are the maximum
%after activation values that correspond to a less than 10% change from
%dilute phase intensity before

% Split bars based on unfolded and folded populations
csatbar(:,1)=mInt(idx).*pU;
csatbar(:,2)=mInt(idx).*pF;

for i=1:length(dG)
    cnames{i}=diro{idx(i)};
end

%return
da2=load(['../../data/optoDroplets_csat/Barnase_unfolding_mutants_cytoplasm_mean_csat.mat']);
dGo2=da2.dG;
mInt2=da2.mInt;
sInt2=da2.sInt;
diro2=da2.dir2;
%diro2={'L14D, L42D, I51D, L63D, I76D, I88D, L89D, I96D','L14S, L42S, I51S, L63S, I76S, I88S, L89S, I96S','L14A, L42A, I51A, L63A, I76A, I88A, L89A, I96A','L14G, L42G, I51G, L63G, I76G, I88G, L89G, I96G','L14D, I51D, I88D','L14S, I51S, I88S','L14A, I51A, I88A','L14G, I51G, I88G'};

%D=[dGo2; mInt2; sInt2]
%xlswrite('csat_unfolding_mutants_data.xls',D)
%return

dG2=dGo2;
figure(142);
errorbar(mInt2,-dG2,sInt2,'horizontal','ok','markersize',8,'markerfacecolor','k'); hold on;
for i=1:length(mInt2)
    text(mInt2(i)+1000,-dG2(i)+0.5,diro2{(i)},'color','k');
end

%return

% Stacked bar based on unfolded population
[dG2, idx]=sort(dG2);
pF2=exp(-(dG2+dGoffset)./(R*T))./(1+exp(-(dG2+dGoffset)./(R*T)));
pU2=1-pF2;

mInt2(end)=32208; % Minimum possible csat for 7D from plot_barnase_unfolding_mutants_dilute_B_A_compare.m

csatbar(length(pU)+3:length(pU)+length(pU2)+2,1)=mInt2(idx).*pU2;
csatbar(length(pF)+3:length(pF)+length(pF2)+2,2)=mInt2(idx).*pF2;
for i=1:length(dG2)
    cnames{length(pU)+2+i}=diro2{idx(i)};
end


f2=figure;
f2.Position = [100 100 540 700];
mycolor=[0.8275 0.8275 0.8275; 0.7059 0.8510 0.6235];
H=barh(csatbar,'stacked'); hold on;
for i = 1:2
    H(i).FaceColor = 'flat';
    H(i).CData = mycolor(i,:);
end
set(gca,'YDir','reverse');
plot([cstarU cstarU],[0 length(csatbar)+2],'--k','linewidth',2); 
plot([mincstarU mincstarU],[0 length(csatbar)+2],'-','color',[0.8 0.8 0.8]); 
plot([maxcstarU maxcstarU],[0 length(csatbar)+2],'color',[0.8 0.8 0.8]);

xlim([0 max(sum(csatbar,2))+5000]);
ylim([-2 size(csatbar,1)+1])
xlabel('c_{sat}');
set(gca,'ytick',[1:1:length(csatbar)]);
set(gca,'yticklabel',cnames);

%% Get predicted csat values
bins=[0.4:0.05:1];
mycolor=cool(length(bins));

ypred=(cstarU./(exp(-(-dG-dGoffset)./(R*T))./(1+exp(-(-dG-dGoffset)./(R*T)))));
ypred2=(cstarU./(exp(-(-dG2-dGoffset)./(R*T))./(1+exp(-(-dG2-dGoffset)./(R*T)))));

% When saved figure:
% maxcstarU = 1.1353e04, maxdGoffset = -1.3170e04
% mincstarU = 9.1033e04, mindGoffset = -1.1997e04
ypredmax=[(maxcstarU./(exp(-(-dG+maxdGoffset)./(R*T))./(1+exp(-(-dG+maxdGoffset)./(R*T))))) (maxcstarU./(exp(-(-dG2+maxdGoffset)./(R*T))./(1+exp(-(-dG2+maxdGoffset)./(R*T)))))]/1000;
ypredmin=[(mincstarU./(exp(-(-dG+mindGoffset)./(R*T))./(1+exp(-(-dG+mindGoffset)./(R*T))))) (maxcstarU./(exp(-(-dG2+maxdGoffset)./(R*T))./(1+exp(-(-dG2+maxdGoffset)./(R*T)))))]/1000;

ypredall=[ypred(4:end) ypred2]/1000;
pUall=[pU(4:end) pU2];
mIntall=[mInt(4:end) mInt2]/1000;
sIntall=[sInt(4:end) sInt2]/1000;


%% Plot in reference to predicted csat - color by camsol
camsoldata=importdata('../../data/optoDroplets_csat/barnase_variants_dG_mcsat_scsat_camsol.txt');
normcamsol=camsoldata.data(:,4)./camsoldata.data(1,4);
varlist=camsoldata.rowheaders;

f=figure;
f.Position=[100 100 550 200];
ymin=0.6;
ymax=1.4;
ylim([ymin ymax])
xlim([0 1])
rectangle('Position',[mean([pU(3),pU(4)]) ymin mean([pU(6),pU(7)])-mean([pU(3),pU(4)]) ymax-ymin],'facecolor',[211 211 211]/255,'edgecolor',[211 211 211]/255); hold on;
rectangle('Position',[mean([pU(6),pU(7)]) ymin 1-mean([pU(6),pU(7)]) ymax-ymin],'facecolor',[100 139 216]/255,'edgecolor',[100 139 216]/255); hold on;
rectangle('Position',[0 ymin mean([pU(3),pU(4)]) ymax-ymin],'facecolor',[229 132 133]/255,'edgecolor',[229 132 133]/255); hold on;
plot([0 1],[1 1],'-k','linewidth',4);
plot(pU(4:end),ypredmax(4:13)./ypredall(1:10),'--k');
plot(pU(4:end),ypredmin(4:13)./ypredall(1:10),'--k');


bins=1.01:0.02:1.09;
%mycolor3=jet(length(bins));
mycolor3=[255 255 255; 180 217 159; 79 124 50; 51 94 23; 29 61 9]./[255 255 255];
%bins=nanmin(blobchange):0.015:nanmax(blobchange);
%mycolor3=[255 255 255; 180 217 159; 79 124 50; 51 94 23; 29 61 9]./[255 255 255];
%mycolor3=bluewhitered(length(bins));

for i=4:length(mInt)
    pos=find(strcmp(varlist,diro(i))==1);
    currcamsol=normcamsol(pos)
    if pU(i)<mean([pU(6),pU(7)])
        plot(pU(i),mInt(i)/ypred(i),'o','markersize',10,'markerfacecolor','k','markeredgecolor','k'); hold on;
        text(pU(i),mInt(i)/ypred(i),diro{i});
    else
        if currcamsol<=bins(1)
            idx=1;
        elseif currcamsol>=bins(end)
            idx=length(bins);
        else
            [~,idx]=histc(currcamsol,bins);
        end
        plot(pU(i),mInt(i)/ypred(i),'o','markersize',10,'markerfacecolor',mycolor3(idx,:),'markeredgecolor','k'); hold on;
        text(pU(i),mInt(i)/ypred(i),diro{i});
    end
    clear currcamsol; clear pos; 
end
ylabel('Measured csat/Predicted csat');
xlabel('pU');

colormap(mycolor3);
cbh = colorbar ; %Create Colorbar
caxis([min(bins) max(bins)])
cbh.Ticks = min(bins):0.02:max(bins); %Create 8 ticks from zero to 1
cbh.TickLabels = 1.01:0.02:1.09;

