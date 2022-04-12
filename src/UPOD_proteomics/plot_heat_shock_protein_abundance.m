clear all;

%% Plot Figure 6C
% Data files needed: 
% Table_S5.xlsx - abundance of top 100 proteins
% Shemesh_NatCom_2021_Supp_1.xlsx - list of known chaperones from https://doi.org/10.1038/s41467-021-22369-9
% Gene_mapping_list.txt - UniProt ID to gene name mapping

% Load in proteomics data
[da,txt]=xlsread('../../data/UPOD_proteomics/Table_S5.xlsx','top_100');
muts=txt(1,2:end);
proteins=txt(2:end,1);

for i=1:length(muts)
    tmp=strsplit(muts{i},'_');
    muts{i}=tmp{1};
end

% Load in known chaperones list
[da2,txt2]=xlsread('../../data/UPOD_proteomics/Shemesh_NatCom_2021_Supp_1.xlsx','Sheet1');
gnames=txt2(:,2);

% Load in gene names of top 100 proteins
da3=importdata('../../data/UPOD_proteomics/Gene_mapping_list.txt');
for i=1:length(da3)
    tmp=strsplit(da3{i},'\t');
    gnamesp{i}=tmp{2};
    acc{i}=tmp{1};
end

intersect(gnames,gnamesp)
%return

% If just focus on high abundance chaperones
%hspid={'mTFP1','P0DMV9;P0DMV8','P04792','Q99832','Q8NBS9'};
%hspname={'mTFP','HSPA1A/B','HSPB1','CCT7','TXNDC5'};

% If just focus on chaperones with significant changes in enrichment
hspid={'mTFP1','P0DMV9;P0DMV8','P04792','Q99832','P48723'};
hspname={'mTFP','HSPA1A/B','HSPB1','CCT7','HSPA13'};
mycolor2=[183 34 37; 193 118 41; 186 183 51; 58 127 62; 0 0 0; 205 76 83; 141 92 166; 46 68 150]./[250 250 250]; 

bvar={'WT','9S','3S','FY','8ala','L14A','Ex4','4Y'};
f=figure;
f.Position=[100 100 800 100];
for i=1:length(hspid)
    pos=find(strcmp(proteins,hspid(i))==1);
    subplot(1,5,i);
    for m=1:length(bvar)
        pos3=find(strcmp(muts,bvar(m))==1);
        tmp=da(pos,pos3);
        mabund(i,m)=mean(log10(tmp));
        sabund(i,m)=std(log10(tmp))/sqrt(4);
        if m==1 | m==6 | m==7
            errorbar(m,mabund(i,m),sabund(i,m),'s','markeredgecolor','k','markerfacecolor',mycolor2(m,:),'color',mycolor2(m,:),'markersize',9); hold on; 
        else
            errorbar(m,mabund(i,m),sabund(i,m),'o','markeredgecolor','k','markerfacecolor',mycolor2(m,:),'color',mycolor2(m,:),'markersize',9); hold on; 
        end
        clear tmp; clear pos3; 
    end
    ylabel('log_{10}(Abundance)');
    set(gca,'xtick',1:1:length(bvar));
    set(gca,'xticklabel',bvar);
    xlim([0 length(bvar)+1]);
    ylim([min(mabund(i,:))-0.15 max(mabund(i,:))+0.15]);
end