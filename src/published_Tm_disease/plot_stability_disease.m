clear all;

% Plots Figure 1A and 1B

%% Load in data

[da,txt]=xlsread('../../data/published_Tm_disease/Table_S1.xlsx','ProThermDB_Homo_Sapiens');
iD1=txt(2:end,6);
Tm=da(:,26);

[da2,txt2]=xlsread('../../data/published_Tm_disease/Table_S1.xlsx','UNIPROT_to_KEGG');
uiDmap=txt2(:,1);
kiDmap=txt2(:,2);

[da3,txt3]=xlsread('../../data/published_Tm_disease/Table_S1.xlsx','KEGG_disease_genes');
kgenes=txt3(:,1);

%% Get unstable vs medium stable vs stable proteins
[sTm,idx]=sort(Tm);
siD1=iD1(idx);

% Tm data come from Leuenberger et al., DOI: 10.1126/science.aai7825
% Leuenberger et al., classified the 10% of proteins with the lowest Tms as unstable, the 10% of proteins with the highest Tms as stable, and the remaining proteins as having medium stability

num10=round(length(Tm)*0.1);
for i=1:length(sTm)
    if i<=num10
        type{i}='U';
    elseif i>num10 & i<=length(Tm)-num10
        type{i}='M';
    else
        type{i}='S';
    end
end

%% Plot Tm probability density estimates - Figure 1A
mycolor=[211 211 211; 180 217 159; 57 125 61]/255;

f=figure;
%f.Position=[100 100 150 60];
pts=35:0.1:80;

% Unstable proteins
pos=find(strcmp(type,'U')==1);
length(pos)
[f,xi] = ksdensity(Tm(pos),pts); 
subplot(3,1,3)
yp = (f + abs(f))/2;
area(xi,yp,'FaceColor',mycolor(1,:),'edgecolor',[1 1 1])
xlim([35 80]);
ylim([0 0.3]);

% Medium stable proteins
pos=find(strcmp(type,'M')==1);
length(pos)
[f,xi] = ksdensity(Tm(pos),pts);
subplot(3,1,2)
yp = (f + abs(f))/2;
area(xi,yp,'FaceColor',mycolor(2,:),'edgecolor',[1 1 1])
xlim([35 80]);
ylim([0 0.3]);

% Stable proteins
pos=find(strcmp(type,'S')==1);
length(pos)
[f,xi] = ksdensity(Tm(pos),pts); 
subplot(3,1,1)
yp = (f + abs(f))/2;
area(xi,yp,'FaceColor',mycolor(3,:),'edgecolor',[1 1 1])
xlim([35 80]);
ylim([0 0.3]);
%return


%% Find genes that are associated with disease and have Tm data for them
count=0;
for i=1:length(kiDmap)
    pos=find(strcmp(kgenes,kiDmap{i})==1);
    if isempty(pos)==0 & length(pos)==1
        count=count+1;
        diseaseacc{count}=uiDmap{i};
    end
    clear pos; 
end

%% Get Tm of accs associated with disease and those not
count=0;
count2=0;
for i=1:length(sTm)
    pos=find(strcmp(diseaseacc,siD1{i})==1);
    if isempty(pos)==0 & length(pos)==1
        count=count+1;
        disTm(count)=sTm(i);
        poslist(count)=i;
        disacc{count}=siD1{i};
        distype{count}=type{i};
    else
        count2=count2+1;
        nodisTm(count2)=sTm(i);
        nodistype{count2}=type{i};
        nodisacc{count2}=siD1{i};
    end
end

%% Find number of unstable, medium stable, and stable in the diease and non disease class - Figure 1B
clear fracdata; 

numdU=length([find(strcmp(distype,'U')==1)]);
numndU=length([find(strcmp(nodistype,'U')==1)]);
fracdata(1,1)=numdU/(numdU+numndU);
fracdata(1,2)=numndU/(numdU+numndU);

numdM=length([find(strcmp(distype,'M')==1)]);
numndM=length([find(strcmp(nodistype,'M')==1)]);
fracdata(2,1)=numdM/(numdM+numndM);
fracdata(2,2)=numndM/(numdM+numndM);

numdS=length(find(strcmp(distype,'S')==1));
numndS=length(find(strcmp(nodistype,'S')==1));
fracdata(3,1)=numdS/(numdS+numndS);
fracdata(3,2)=numndS/(numdS+numndS);

ctable=[numdU numndU; numdS numndS];
[h, p, stats]=fishertest(ctable)

f=figure;
%f.Position=[100 100 150 200];
%f.Position=[100 100 150 120];
mycolor=[211 211 211; 180 217 159; 57 125 61]/255;
for i=1:3
    h = bar(i, fracdata(i,1)); hold on; 
    set(h, 'FaceColor', mycolor(i,:));
end
ylabel('Fraction Disease Proteins');
set(gca,'XTick',1:1:3);
set(gca,'XTickLabel',{'U','M','S'});
%print -painters -depsc 'Fraction_Disease_Proteins_by_stability.eps'

%return
%% Get unstable disease protein ids
pos=[find(strcmp(distype,'U')==1)];
currlist=disacc(pos);
currTm=disTm(pos);

