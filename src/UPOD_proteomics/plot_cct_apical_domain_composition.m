clear all;

%% Plot Figure S5C

% Data files needed:
% cct.fasta
% cct_apical_domain_start_and_end.txt

%mycolors={'#0f8041', '#118c42', '#159e4c', '#1cb75a','#20d368','#000000','#262626','#3f3f3f','#565656','#707070','#3b54a4','#4568bf','#ea2729','#f94851','#ffce09','#ffd939','#ffe164','#b8539f','#d663be','#e86dce'}

aas={'F','Y','W','L','I','V','M','A','Q','S','N','T','G','K','R','E','D','P','H','C'};
mycolors=[255, 206, 9; 255, 217, 57; 255, 225, 100; 0 0 0; 38, 38, 38; 63, 63, 63; 86, 86, 86; 112, 112, 112; 15, 128, 65; 17, 140, 66; 21, 158, 76; 28, 183, 90; 32, 211, 104; 59, 84, 164; 69, 104, 191; 234, 39, 41; 249, 72, 81; 184, 83, 159; 214, 99, 190; 232, 109, 206]./[255 255 255];

[hs, seqs]=fastaread('../../data/UPOD_proteomics/cct.fasta'); % sequence in fasta format

da=importdata('../../data/UPOD_proteomics/cct_apical_domain_start_and_end.txt');
uniID=da.textdata(2:end,3);
st=da.data(:,1);
en=da.data(:,2);

names={'CCT1','CCT2','CCT3','CCT4','CCT5','CCT6a','CCT7','CCT8'};

for s=1:length(seqs)
    currseq=seqs{s};
    fullhead=hs{s};
    tmp=strsplit(fullhead,'|');
    unihead=tmp{2};
    pos=find(strcmp(uniID,unihead)==1);
    subseq=currseq(st(pos):en(pos));
    for a=1:length(aas)
        tmp2=strfind(subseq,aas(a));
        fracaa(s,a)=length(tmp2)/length(subseq);
        clear tmp2;
    end
    %return
    clear currseq; clear fullhead; clear tmp; clear unihead; clear pos; 
    clear subseq;
end

[~,idx]=sort(sum(fracaa(:,1:3),2));
for i=1:length(idx)
    sfracaa=fracaa(idx,:);
end

figure;
ba=bar(sfracaa,'stacked', 'FaceColor', 'Flat');
for a=1:length(aas)
    ba(a).FaceColor=mycolors(a,:);
end
set(gca,'xticklabel',names(idx));
legend(aas)
ylabel('Amino Acid Fraction')


