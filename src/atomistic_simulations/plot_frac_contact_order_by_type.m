clear all;

% Plot Figure 7A

dthresh=5;

types{1}={'GSTNQ'};
types{2}={'RKH'};
types{3}={'DE'};
types{4}={'CP'}; 
types{5}={'FWY'};
types{6}={'MVILA'};
mycolor3=[0 0.5 0; 0 0 1; 1 0 0; 1 0 1; 1 0.8 0; 0.0 0.0 0.0];

seqs={'AQVINTFDGVADYLQTYHKLPDNYITKSEAQALGWVASKGNLADVAPGKSIGGDIFSNREGKLPGKSGRTWREADINYTSGFRNSDRILYSSDWLIYKTTDAYQTFTKIR',...
'ATKAVCVLKGDGPVQGIINFEQKESNGPVKVWGSIKGLTEGLHGFHVHEFGDNTAGCTSAGPHFNPLSRKHGGPKDEERHVGDLGNVTADKDGVADVSIEDSVISLSGDHCIIGRTLVVHEKADDLGKGGNEESTKTGNAGSRLACGVIGIAQ',...
'MTERRVPFSLLRGPSWDPFRDWYPHSRLFDQAFGLPRLPEEWSQWLGGSSWPGYVRPLPPAAIESPAVAAPAYSRALSRQLSSGVSEIRHTADRWRVSLDVNHFAPDELTVKTKDGVVEITGKHEERQDEHGYISRCFTRKYTLPPGVDPTQVSSSLSPEGTLTVEAPMPKLATQSNEITIPVTFESRAQLGGPEAAKSDETAAK',...
'MPPYTVVYFPVRGRCAALRMLLADQGQSWKEEVVTVETWQEGSLKASCLYGQLPKFQDGDLTLYQSNTILRHLGRTLGLYGKDQQEAALVDMVNDGVEDLRCKYISLIYTNYEAGKDDYVKALPGQLKPFETLLSQNQGGKTFIVGDQISFADYNLLDLLLIHEVLAPGCLDAFPLLSAYVGRLSARPKLKAFLASPEYVNLPINGNGKQ',...
'MSSGNAKIGHPAPNFKATAVMPDGQFKDISLSDYKGKYVVFFFYPLDFTFVCPTEIIAFSDRAEEFKKLNCQVIGASVDSHFCHLAWVNTPKKQGGLGPMNIPLVSDPKRTIAQDYGVLKADEGISFRGLFIIDDKGILRQITVNDLPVGRSVDETLRLVQAFQFTDKHGEVCPAGWKPGSDTIKPDVQKSKEYFSKQK',...
'MDTSRVQPIKLARVTKVLGRTGSQGQCTQVRVEFMDDTSRSIIRNVKGPVREGDVLTLLESEREARRLR',...
'MARTKQTARKSTGGKAPRKQLATKAARKSAPSTGGVKKPHRYRPGTVALREIRRYQKSTELLIRKLPFQRLVREIAQDFKTDLRFQSAAIGALQEASEAYLVGLFEDTNLCAIHAKRVTIMPKDIQLARRIRGERA',...
'ATLEKLMKAFESLKSFQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQPPPPPPPPPPPQLPQPPPQAQPLLPQPQPPPPPPPPPPGPAVAEEPLHRP'};

mydirs={'../../data/atomistic_simulations/WT/335',...
'../../data/atomistic_simulations/SOD1/335',...
'../../data/atomistic_simulations/HSPB1/335',...
'../../data/atomistic_simulations/GSTP1/335',...
'../../data/atomistic_simulations/PRDX1/335',...
'../../data/atomistic_simulations/RPS28/335',...
'../../data/atomistic_simulations/H33A/335',...
'../../data/atomistic_simulations/Q49'};
numreps=[5 5 5 5 5 5 5 5 5];

figure;
for d=1:length(mydirs) 
        if d<length(mydirs)
	    for r=1:numreps(d)
	        conord(r,:)=load([mydirs{d} '/' num2str(r) '/ana/contact_order_' num2str(dthresh) '.000.csv']);
	    end
        else
            for r=1:numreps(d)
                conord(r,:)=load([mydirs{d} '/' num2str(r) '/335/ana/contact_order_' num2str(dthresh) '.000.csv']);
            end
        end

	mconord=mean(conord);
	sconord=std(conord);

        myseq=seqs{d};
	for i=1:length(types)
	    tmp=types{i};
	    tmp=tmp{1};
	    pos=[];
	    for j=1:length(tmp)
		pos=[pos strfind(myseq,tmp(j))];
	    end
	    mcoaas(d,i)=sum(mconord(pos))/sum(mconord);
            tmp2=length(pos)*(sum(mconord)/length(myseq));
            expmco(d,i)=tmp2/sum(mconord);
            fractype(d,i)=length(pos)/length(myseq);

	    if mcoaas(d,i)>0
		plot(d,i,'o','markersize',100*mcoaas(d,i),'markerfacecolor',mycolor3(i,:),'markeredgecolor','k'); hold on; 
                plot(d,i,'o','markersize',100*expmco(d,i),'markeredgecolor',[0.5 0.5 0.5]); hold on; 
	    end
	%return
	clear tmp; clear pos; clear tmp2;  
	end
	clear conord; clear mconord; clear sconord; 
end


set(gca,'ytick',1:1:length(types));
set(gca,'yticklabel',{'Polar','Positive','Negative','Cys / Pro','Aromatic','Aliphatic'});
ylim([0 7])

xlim([0 length(mydirs)+1]);
set(gca,'xtick',1:1:length(mydirs));
set(gca,'xticklabel',{'Barnase','SOD1','HSPB1','GSTP1','PRDX1','RPS28','H3-3A','Httex1-Q49'});
xtickangle(45);

