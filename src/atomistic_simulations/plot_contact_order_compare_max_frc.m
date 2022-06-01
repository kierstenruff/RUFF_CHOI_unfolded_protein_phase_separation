clear all;

% Plot Figure 4A

constructs={'WT'};
mytitle={'WT'};

dthresh=5;

Ts=[335];
reps=[1 2 3 4 5];
repsfrc=[1 2 3];

types{1}={'GSTNQ'};
types{2}={'RKH'};
types{3}={'DE'};
types{4}={'CP'}; 
types{5}={'FWY'};
types{6}={'MVILA'};
mycolor3=[0 0.5 0; 0 0 1; 1 0 0; 1 0 1; 1 0.8 0; 0 0 0];


[h myseqs]=fastaread(['../../data/atomistic_simulations/barnase_seqs.fasta']);

for c=1:length(constructs)
    figure(802+c);
    pos=find(strcmp(h,constructs{c})==1);
    myseq=myseqs{pos};
    
    % Full Hamiltonian
    for r=1:length(reps)
        conord(r,:)=load(['../../data/atomistic_simulations/' constructs{c} '/' num2str(Ts) '/' num2str(reps(r)) '/ana/contact_order_' num2str(dthresh) '.000.csv']);
        mcor(r)=mean(conord(r,:));
    end

    mconord(c,:)=mean(conord);
    sconord(c,:)=std(conord);

    subplot(4,1,1)
    bar(mconord(c,:)); hold on;

    clear mcor; 

    % For FRC
    for r=1:length(repsfrc)
        conord(r,:)=load(['../../data/atomistic_simulations/' constructs{c} '/FRC/' num2str(repsfrc(r)) '/ana/contact_order_' num2str(dthresh) '.000.csv']);
        mcofrcr(r)=mean(conord(r,:));
    end

    mconordfrc(c,:)=mean(conord);
    sconordfrc(c,:)=std(conord);

    subplot(4,1,2)
    bar(mconordfrc(c,:)); hold on;
    clear mcofracr; 

    subplot(4,1,3)
    bar(mconord(c,:)-mconordfrc(c,:)); hold on;

    subplot(4,1,4)
    bar(mconord(c,:)./mconordfrc(c,:)); hold on;

    for i=1:length(types)
        tmp=types{i};
        tmp=tmp{1};
        for j=1:length(tmp)
            pos=strfind(myseq,tmp(j));
            subplot(4,1,1)
            plot(pos,mconord(c,pos),'o','markerfacecolor',mycolor3(i,:),'markeredgecolor','k','markersize',6); hold on;

            subplot(4,1,2)
            plot(pos,mconordfrc(c,pos),'o','markerfacecolor',mycolor3(i,:),'markeredgecolor','k','markersize',6); hold on;

            subplot(4,1,3)
            plot(pos,(mconord(c,pos)-mconordfrc(c,pos)),'o','markerfacecolor',mycolor3(i,:),'markeredgecolor','k','markersize',6); hold on;
            subplot(4,1,4)
            plot(pos,(mconord(c,pos)./mconordfrc(c,pos)),'o','markerfacecolor',mycolor3(i,:),'markeredgecolor','k','markersize',6); hold on;
           clear pos; 
        end
    end



    % Find the contact order that corresponds to a given normalized cut off
    minval=max(mconordfrc(c,:));
    bins=0:0.01:minval;
    mycolor=gray(length(bins));

    % Sets contact order equivalent to max frc to 10
    maxco=max(mconordfrc(c,:));
    subplot(4,1,1)
    plot([0.5 length(myseq)+0.5],[maxco maxco],'--');
    myscal=10/maxco;

    f=figure(903+c);
    f.Position=[100 100 1500 150];

    for i=1:length(types)
        tmp=types{i};
        tmp=tmp{1};
        for j=1:length(tmp)
            pos=strfind(myseq,tmp(j));
            for p=1:length(pos)
                plot(pos(p),0.5,'o','markerfacecolor',mycolor3(i,:),'markeredgecolor','k','markersize',myscal*(mconord(c,pos(p)))); hold on;
                if mconord(c,pos(p))>max(mconordfrc(c,:))
                    idx=1;
                    text(pos(p)-0.5,1,myseq(pos(p)),'color',mycolor(idx,:),'fontweight','bold');
                else
                    tmp2=minval-(mconord(c,pos(p)));
                    if tmp2>bins(end)
                        idx=length(bins);
                        text(pos(p)-0.5,1,myseq(pos(p)),'color',mycolor(idx,:));
                    else
                        [~,idx]=histc(tmp2,bins);
                        text(pos(p)-0.5,1,myseq(pos(p)),'color',mycolor(idx,:));
                    end
                    clear tmp2;  
	        end
            end
            clear pos; 
        end
        clear tmp; 
    end
    xlim([0.5 length(myseq)+0.5])
    %ylim([0.25 0.75])

    pos=find(mconord(c,:)>maxco)
    myseq(pos)

    % Create legend
    figure;
    plot(0.5,1,'o','markeredgecolor','k','markersize',0.03*myscal); hold on; 
    plot(0.75,1,'o','markeredgecolor','k','markersize',0.05*myscal);
    plot(1,1,'o','markeredgecolor','k','markersize',0.07*myscal);
    title(mytitle(c));
    set(gca,'xtick',[0.5 0.75 1]);
    set(gca,'xticklabel',{'0.03','0.05','0.07'})
end




