clear all;

da=load(['../../data/optoDroplets_csat/Barnase_hsp_inhibition_mean_csat.mat']);
mInt=da.mInt;
sInt=da.sInt;
diro=da.dir2;

csat=mInt/1000;
csaterr=sInt/1000;

mycolor=[0.8 0.8 0.8];

f=figure;
f.Position=[100 100 600 200];
h=bar(1, csat(1)); hold on; 
set(h, 'FaceColor', [211 211 211]/255)
errorbar(1,csat(1),csaterr(1),'ok'); hold on; 

h=bar(3, csat(2)); hold on; 
set(h, 'FaceColor', [225 242 214]/255)
errorbar(3,csat(2),csaterr(2),'ok'); hold on; 
h=bar(4, csat(3)); hold on; 
set(h, 'FaceColor', [181 218 160]/255) 
errorbar(4,csat(3),csaterr(3),'ok'); hold on; 
h=bar(5, csat(4)); hold on; 
set(h, 'FaceColor', [58 126 62]/255) 
errorbar(5,csat(4),csaterr(4),'ok'); hold on; 

h=bar(7, csat(5)); hold on; 
set(h, 'FaceColor', [200 183 226]/255) 
errorbar(7,csat(5),csaterr(5),'ok'); hold on; 
h=bar(8, csat(6)); hold on; 
set(h, 'FaceColor', [154 130 191]/255) 
errorbar(8,csat(6),csaterr(6),'ok'); hold on; 
h=bar(9, csat(7)); hold on; 
set(h, 'FaceColor', [104 69 148]/255)
errorbar(9,csat(7),csaterr(7),'ok'); hold on; 

xlim([0 10])
set(gca,'XTick',[1 3:5 7:9]);
set(gca,'XTickLabel',{'L14A','1','5','10','1','5','10'});
ylabel('csat (a.u.)');
ylim([0 35])
set(gca,'YTick',[5:10:35]);

plot([0 10],[csat(1) csat(1)],'-k');



