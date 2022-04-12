clear all;

% Plots Figure 1C

[da,txt]=xlsread('../../data/published_Tm_disease/Table_S1.xlsx','Unstable_proteins_disease');

iD=txt(2:end,1);
distype=txt(2:end,5);
gene=txt(2:end,4);
Tm=da;

udistype=unique(distype);

mycolor=[0 0 0; 211 211 211; 57 125 61; 180 217 159; 204 76 82; 232 163 166; 48 64 150; 121 149 204; 104 69 147]./255;

count=0;
f=figure;
f.Position=[100 100 800 200];
for d=1:length(udistype)
    pos=find(strcmp(distype,udistype{d}));
    currTm=Tm(pos);
    [scurrTm,idx]=sort(currTm);
    scurrIDs=iD(pos(idx));
    scurrgenes=gene(pos(idx));
    for g=1:length(pos)
        count=count+1;
        h = bar(count, currTm(g)); hold on; 
        set(h, 'FaceColor', mycolor(d,:));
        genelist(count)=scurrgenes(g);
    end
    count=count+1;
end
ylim([37 46])
set(gca,'Xtick',1:1:count)
set(gca,'XtickLabel',genelist)
xtickangle(90)

%print -painters -depsc 'unstable_proteins_associated_w_disease_bar.eps'
