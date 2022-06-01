clear all; 
aas='ARNDCQEGHILKMFPSTWYV';
KD=[0.700 0.000 0.111 0.111 0.778 0.111 0.111 0.456 0.144 1.000 0.922 0.067 0.711 0.811 0.322 0.411 0.422 0.400 0.356 0.967];

binsz=3;
minblsz=4;
seqspace=20;

[h seqs]=fastaread('../../data/optoDroplets_csat/Barnase_seqs_conservative.txt');

f=figure(40); 
f.Position = [100 100 540 700];
for s=1:length(seqs)
    if s==1
        wtseq=seqs{s};
    end
    currseq=seqs{s}
    
    nummutations(s)=length(find((wtseq==currseq)==0));
    
    for i=1:length(currseq)
        pos=strfind(aas,currseq(i));
        hydrovec(i)=KD(pos);
        clear pos; 
    end

    for i=1:length(currseq)-binsz+1
        meanHydro(i)=mean(hydrovec(i:i+2));
    end
    
    pos=find(meanHydro>0.37);
    mean01=zeros(1,length(meanHydro));
    mean01(pos)=1;
    clear pos;
    
    %figure(39);
    %plot(meanHydro); hold on;
    if s==1
        wtvec=meanHydro;
    end
    tmpvec=meanHydro;
    mmtmp(s)=mean(wtvec-tmpvec);
    
    figure(40);
    plot([1 length(meanHydro)],[seqspace*(s-1) seqspace*(s-1)],'-k','linewidth',2);
    %% Get hydrophobic blobs
    count=0;
    for i=1:length(mean01)
        if i==1
            lenhydro=0;
        end
        if mean01(i)==1
            if lenhydro==0
                sthydro=i;
            end
            lenhydro=lenhydro+1;
        else 
            % Check to see if length hydro is above minblsz
            if lenhydro>=minblsz
                circpos=[sthydro seqspace*(s-1)-(lenhydro-1)/2 (lenhydro-1) (lenhydro-1)];
                rectangle('Position',circpos,'Curvature',[1 1],'EdgeColor',[0.2235 0.4824 0.2314],'FaceColor',[0.2235 0.4824 0.2314]); hold on;
                count=count+1;
                szhblob(s,count)=lenhydro;
                sthblob(s,count)=sthydro; 
                seqblob01(s,sthydro:sthydro+lenhydro-1)=1;
            end
            lenhydro=0;
            sthydro=0;
        end
    end
    
    %% Get polar blobs
    count=0;
    for i=1:length(mean01)
        if i==1
            lenpolar=0;
        end
        if mean01(i)==0
            if lenpolar==0
                stpolar=i;
            end
            lenpolar=lenpolar+1;
        else 
            % Check to see if length hydro is above minblsz
            if lenpolar>=minblsz
                circpos=[stpolar seqspace*(s-1)-(lenpolar-1)/2 (lenpolar-1) (lenpolar-1)];
                rectangle('Position',circpos,'Curvature',[1 1],'EdgeColor',[0.1882 0.2275 0.60],'FaceColor',[0.1882 0.2275 0.60]); hold on;
                count=count+1;
                szpblob(s,count)=lenpolar;
                stpblob(s,count)=stpolar; 
                seqblob01(s,stpolar:stpolar+lenpolar-1)=-1;
            end
            lenpolar=0;
            stpolar=0;
        end
    end
    xlim([0 length(meanHydro)+1])
    
    
    clear hydrovec; clear meanHydro; clear mean01; 
end

set(gca,'YTick',[0:seqspace:seqspace*(s-1)]);
set(gca,'YTickLabel',h);
set(gca,'YDir','reverse');
ylim([0-seqspace seqspace*(s)])

maxchange=10;
for s=2:length(seqs)
    tmp=sum(seqblob01(1,:)-seqblob01(s,:));
    qblobpickprob(s-1)=(maxchange-tmp)/maxchange; 
    valchange(s-1)=tmp;
end

for s=2:length(h)
    seqsnames{s-1}=h{s};
end

% valchange for dG vs sticker disruption schematic
vcall=[0 valchange]+2*nummutations;

%save('../../data/optoDroplets_csat/Barnase_blob_change_pick_probability.mat','seqsnames','qblobpickprob');