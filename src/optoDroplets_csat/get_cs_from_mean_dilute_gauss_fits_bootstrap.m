clear all;

% Uses data generated from get_dilute_phase_intensity_PS_by_gauss_fit.m

myloc={'cytoplasm'};
minco=1000; % used to determine if fullDiluteIntB-fullDiluteIntA is around 1to1 line
nTrials=5%50;

%dir1={'Barnase'};
%dir2={'Ex2','Ex4','I25A','I88G','L14A','WT'};
%myorder=[1 2 3 4 5 6];

dir1={'Barnase_additional_mutants'};
dir2={'WT','I55G','V45T','I25A','I88A','L14A','I51A','L89G','I96G','Ex3','I88G','Ex2','Ex4'};
%dG=[-24.9 -18.7 -17.0 -13.0 -12.1 -11.7 -10.3 -9.6 -8.3 -3.7 -3.2 -0.7 0.8];  
myorder=[1:1:length(dir2)];
%dGo=dG(myorder);
%diro=dir2(myorder);

%dir1={'Barnase_sticker_mutants'};
%dir2={'9S','7S','5S','3S','8Ala','F7Y','F106Y','F7Y_F106Y','F7Y_F82Y_F106Y','FY','4Y'};
%myorder=[1:1:length(dir2)];

%dir1={'Barnase_unfolding_mutants'}; 
%dir2={'3D','3S','3A','3G','7A','7S','7D'}; 
%dG=[8.3 10.4 10.5 25.2 63.0 64.4 74.4];
%[~,myorder]=sort(dG,'ascend');
%myorder=[1:1:length(dir2)];
%dGo=dG(myorder);
%diro=dir2(myorder);

%dir1={'Barnase_hsp_coexpression'};
%dir2={'WT_invEm','WT_hsp40','WT_hsp70','WT_hsp70_hsp40','L14A_invEm','L14A_hsp40','L14A_hsp70','L14A_hsp70_hsp40'};
%myorder=[1:1:length(dir2)];


figure(200);
for d=1:length(dir2)
    da=load(['../../data/optoDroplets_csat/' dir1{1} '/Mean_Dilute_Phase_Intensities_For_Phase_Separated_Cells_by_Gauss_Fit_' dir2{d} '_' myloc{1} '.mat']);
    fullDiluteIntB=da.fullDiluteIntB;
    fullDiluteIntA=da.fullDiluteIntA;
    %return
    if strcmp(dir2{d},'WT')==1 | strcmp(dir2{d},'I25A')==1 | strcmp(dir2{d},'L14A')==1 | strcmp(dir2{d},'Ex2')==1 | strcmp(dir2{d},'Ex4')==1 | strcmp(dir2{d},'I88G')==1
        da2=load(['../../data/optoDroplets_csat/Barnase/Mean_Dilute_Phase_Intensities_For_Phase_Separated_Cells_by_Gauss_Fit_' dir2{d} '_' myloc{1} '.mat']);
        fullDiluteIntA=[fullDiluteIntA da2.fullDiluteIntA];
        fullDiluteIntB=[fullDiluteIntB da2.fullDiluteIntB];
    end
    
    % Make vectors for fitting
    myslope=nan(1,nTrials);
    mycsat=nan(1,nTrials);

    % Make sure no dilute phase concentration before is less than zero
    pos=find(fullDiluteIntB<0);
    fullDiluteIntB(pos)=[];
    fullDiluteIntA(pos)=[];
    clear pos;
    
    pt=find(myorder==d);
    subplot(3,5,pt)
    plot(fullDiluteIntB,fullDiluteIntA,'ok','color','k'); hold on;
    title(dir2{d})
    xlim([0 5e4]);
    ylim([0 5e4]);
    %diro{myorder(d)}=dir2{d};
    %dGo(myorder(d))=dG(d);
    
    % For all points off 1 to 1 line find the 90% of max difference in
    % dilute concentration before and after
    pos=find((fullDiluteIntB-fullDiluteIntA)>minco);
    maxco=0.9*mean(fullDiluteIntB(pos)-fullDiluteIntA(pos));
    clear pos; 
    %return
    
    % Remove outliers
    % First remove points far off diagonal (>90% of max difference) but still in 1by1 regime
    pos2=find(fullDiluteIntB-fullDiluteIntA<minco);
    if isempty(pos2)==0
        pos=find(fullDiluteIntB<=max(fullDiluteIntB(pos2)) & fullDiluteIntB-fullDiluteIntA>maxco);
        plot(fullDiluteIntB(pos),fullDiluteIntA(pos),'om','markerfacecolor','m'); hold on;
        %return
        fullDiluteIntB(pos)=[];
        fullDiluteIntA(pos)=[];
        clear pos; 
    end
    clear pos2; 
    %return
    
    % Next remove remaining on diagonal or intermediate off diagonal points that don't
    % fit
    pos2=find(fullDiluteIntB-fullDiluteIntA<maxco);
    if isempty(pos2)==0
        pos=find(fullDiluteIntB<=max(fullDiluteIntB(pos2)));
        mdl=fitlm(fullDiluteIntB(pos),fullDiluteIntA(pos))
        tmp=find(mdl.Diagnostics.CooksDistance>5*mean(mdl.Diagnostics.CooksDistance))
        plot(fullDiluteIntB(pos(tmp)),fullDiluteIntA(pos(tmp)),'og','markerfacecolor','g'); hold on;
        fullDiluteIntB(pos(tmp))=[];
        fullDiluteIntA(pos(tmp))=[];
        clear pos; clear tmp; clear mdl; 
    end
    clear pos2; 
    
    % Lastly find points way off the diagonal but don't agree with a linear fit
    % of the other points off the diagonal 
    pos2=find(fullDiluteIntB-fullDiluteIntA>=maxco);
    if isempty(pos2)==0
        mdl=fitlm(fullDiluteIntB(pos2),fullDiluteIntA(pos2))
        tmp=find(mdl.Diagnostics.CooksDistance>5*mean(mdl.Diagnostics.CooksDistance))
        plot(fullDiluteIntB(pos2(tmp)),fullDiluteIntA(pos2(tmp)),'ob','markerfacecolor','b'); hold on;
        %return
        fullDiluteIntB(pos2(tmp))=[];
        fullDiluteIntA(pos2(tmp))=[];
        clear mdl; clear tmp; 
    end
    clear pos2; 
    
    %return
    x=0:1e4:5e4;
    y=0:1e4:5e4;
    plot(x,y,'-k'); hold on;
    
    pos=find(isnan(fullDiluteIntB)==1);
    fullDiluteIntB(pos)=[];
    fullDiluteIntA(pos)=[];
    clear pos; 
    %return
    
    % Must have at least 3 off diagonal points to get threshold
    % concentration 
    if length(find(fullDiluteIntB-fullDiluteIntA>=maxco))<=2 
        mInt(myorder(d))=nan;
        sInt(myorder(d))=nan;
    else
        % Identify where to start splitting the data for on diagonal and
        % off diagonal fits, takes into account both on diagonal and off
        % diagonal points 
        p2=find(fullDiluteIntA./fullDiluteIntB>0.9);
        p3=find(fullDiluteIntB-fullDiluteIntA>maxco);
        minsplit=min(0.9*max(fullDiluteIntB(p2)),0.9*min(fullDiluteIntB(p3)));
        
        % Sample 90% of points off the diagonal
        poff=find((fullDiluteIntB-fullDiluteIntA)>=maxco);
        numSamples=round(0.9*length(poff));
        for n=1:nTrials

            myc=0; % counter for unique points
            splitVals=[];
            while myc<=2 | isempty(splitVals)==1 % Need at least 3 unique points sampled to be a good trial
                if length(poff)==3
                    poff2=randsample(length(poff),numSamples,'false');
                else
                    poff2=randsample(length(poff),numSamples,'true');
                end
                % only fitting points really on the diagonal (<minco) and
                % really of the diagonal (>=maxco), intermediate points
                % can make fits go to lower csat values since makes fits
                % more vertical
                pos=[poff(poff2) find((fullDiluteIntB-fullDiluteIntA)<minco)];
                diluteB=fullDiluteIntB(pos);
                diluteA=fullDiluteIntA(pos);           
                myc=length(find(unique(diluteB-diluteA)>maxco));
                
                sdiluteB=sort(diluteB);
                p1=find(diluteB-diluteA<minco);
                splitVals=minsplit:100:sdiluteB(end-2);
                %figure; 
                %plot(diluteB,diluteA,'o'); 
                %return
                clear poff2; clear pos; clear p1; 
            end

            
            %sdiluteB=sort(diluteB);
            %splitVals=sdiluteB(11):100:sdiluteB(end-10);
            %splitVals=sdiluteB(3):100:sdiluteB(end-5);
            
            % Next try to fit sampled data for each potential split in on
            % and off diagonal points to identify where csat is
            for i=1:length(splitVals)
               pos1=find(diluteB<splitVals(i));
               pos2=find(diluteB>=splitVals(i));
               if length(unique(diluteB(pos2)))>=3 & length(unique(diluteB(pos1)))>=1 % must have at least 3 points to fit off diagonal and 1 to fit on diagonal, only need one on diagonal because we are fitting to the 1 to 1 line directly so is more of a error from this rather than a fit
                   % Fit off diagonal points starting at splitVals - fit is
                   % shifted in the x direction by splitVals such that y
                   % intercept is splitVals and only fitting for the slope
                   % off the 1 to 1 line
                   [f2,gof2]=fit(diluteB(pos2)'-splitVals(i),diluteA(pos2)','poly1','Lower',[0.1 splitVals(i)],'Upper',[1 splitVals(i)]);
                   fitErr(i)=gof2.sse;
                   myfits(i,:)=[f2.p1 f2.p2];
                   r2(i)=gof2.rsquare;
                   tmpme(i)=max(diluteB(pos1));
               else
                   fitErr(i)=nan;
                   myfits(i,:)=[nan nan];
                   r2(i)=nan;
               end
               %return
               clear pos1; clear pos2;  clear f2; clear gof2; 
            end

            if isnan(min(fitErr))==0
                % minimize both fitErr and 1-r2 since want the best fit it
                % terms of sum of squares and r2 
                [~,i1]=sort(fitErr);
                [~,i2]=sort(1-r2);
                for f=1:length(fitErr)
                    totp(f)=find(f==i1)+find(f==i2);
                end
                %return
                clear i1; clear i2; 
                
                idxB=find(totp==min(totp));
                pos=find(diluteB>=splitVals(idxB(end)));
                y2=myfits(idxB(end),1).*(x-myfits(idxB(end),2))+myfits(idxB(end),2);
                %figure(200+d);
                %subplot(2,5,n)
                %plot(diluteB,diluteA,'ob'); hold on; 
                %plot(diluteB(pos),diluteA(pos),'or');
                %plot(x,y2,'-r'); 
                %plot(x,x,'-b');
                %return

                figure(200);
                plot(x,y2,'-r'); 

                %find the csat
                x0 = myfits(idxB(end),2);
                myslope(n) = myfits(idxB(end),1);
                mycsat(n) = myfits(idxB(end),2);
                myInt(n,myorder(d))=x0;
            else
                myInt(n,myorder(d))=nan;
            end


            
            %return
            clear fitErr; clear myfits; clear idxB; clear pos; clear x0; clear y2; 
            clear diluteB; clear diluteA; clear splitVals; clear sdiluteB; clear r2; 
            clear totp;   
        end
        mInt(myorder(d))=nanmean(myInt(:,myorder(d)));
        sInt(myorder(d))=nanstd(myInt(:,myorder(d)));%/sqrt(nTrials);
        clear p2; clear p3; clear minsplit; clear poff; 
    end    
    %return
    %save(['../../data/optoDroplets_csat/' dir1{1} '_' myloc{1} '_' dir2{d} '_dilute_ints_outliers_removed.mat'],'fullDiluteIntA','fullDiluteIntB','myslope','mycsat');
    clear da; clear fullDiluteIntA; clear fullDiluteIntB; clear da2; clear maxco; clear myslope; clear mycsat; 
end

figure; bar(mInt); hold on; 
errorbar([1:1:length(mInt)],mInt,sInt,'o');
ylabel('Phase Separation Intensity (AU)')

%return
%save(['../../data/optoDroplets_csat/' dir1{1} '_' myloc{1} '_mean_csat.mat'],'mInt','sInt','dir2','dG','nTrials');
%save(['../../data/optoDroplets_csat/' dir1{1} '_' myloc{1} '_mean_csat.mat'],'mInt','sInt','dir2','nTrials');

