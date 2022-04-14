clear all;

% Get dilute phase intensities before and after activation
type={'Barnase_additional_mutants'};
dir={'Ex2','Ex3','Ex4','I25A','I51A','I55G','I88A','I96G','L14A','L89G','V45T','WT'};
cloc={'cytoplasm'};
ibins=5:0.1:15; 

for d=1:length(dir)
    fullDiluteIntBfromMean=[];
    fullDiluteIntBfromStd=[];
    fullDiluteIntB=[];
    fullDiluteIntA=[];
    fullDenseIntA=[];
    fullgof=[];
    fullHighIntensityCells=[];

    %% Load in data
    da=readtable(['../../data/optoDroplets_csat/' type{1} '/' dir{d} '_' cloc{1} '.csv']);
    uimagenames=unique(da.image_name);
    numImgs=length(uimagenames);
    mycolor=jet(numImgs);

    for i=1:numImgs  
        posc=find(strcmp(da.image_name,uimagenames{i})==1);
        niBeforeAct=(da.before_highPMT(posc)); % Intensity before activation 
        niAfterAct=(da.after_highPMT(posc)); % Intensity after activation 
        iBeforeAct=log(da.before_highPMT(posc)); % Natural log intensity before activaiton 
        iAfterAct=log(da.after_highPMT(posc)); % Natural log intensity after activation 
        cels=da.cell_number(posc); % Get cell names 
        ucels=unique(cels); % Get unique cells
        xpos=da.x(posc);
        ypos=da.y(posc);
    
        max(iBeforeAct)
    
        % Plot cells based on natural log intensities
        for p=1:length(iBeforeAct) 
            bmat(xpos(p),ypos(p))=iBeforeAct(p);
            amat(xpos(p),ypos(p))=iAfterAct(p);
        end
        
        % Plot check
        figure(2051+d); imagesc(bmat);
        caxis([5 12]);
        figure(3051+d); imagesc(amat);
        caxis([5 12]);
        
        % Create a dilute intensities vector
        diluteIntB=nan*ones(1,length(ucels));
        diluteIntA=nan*ones(1,length(ucels));
        denseIntA=nan*ones(1,length(ucels));
        mygof=nan*ones(1,length(ucels));
        highIntCells=nan*ones(1,length(ucels));
        
        for c=1:length(ucels)
            pos=find(strcmp(cels,ucels{c})==1);
            diluteIntBfromMean(c)=mean(exp(iBeforeAct(pos)));
            diluteIntBfromStd(c)=std(exp(iBeforeAct(pos)));
        
            % Plot check
            figure(2051+d);
            text(ypos(pos(1)),xpos(pos(1)),num2str(c),'color','w');
            figure(3051+d);
            text(ypos(pos(1)),xpos(pos(1)),num2str(c),'color','w'); 
        
            % Get histograms of natural log intensities 
            hB=hist(iBeforeAct(pos),ibins);
            hA=hist(iAfterAct(pos),ibins);
            nhB=hB/sum(hB);
            nhA=hA/sum(hA);  
            %return
            % Remove cells that have a large percentage of the intensity at max of microscope
            if sum(nhB(61:end))>0.25
                highIntCells(c)=1;
            else
        
                % Since before activation cells should be more or less the same intensity remove pixels with intensities not around the maximum peak intensity - use gaussian fit to find the mean 
                fun = @(p,ibins)p(1) * exp(-((ibins - p(2))/p(3)).^2);
                %options = optimoptions('lsqcurvefit','OptimalityTolerance',1e-20);
                x0=[1 5 1];
                lb = [];
                ub = [];
                p = lsqcurvefit(fun,x0,ibins,nhB,lb,ub);
                %p = lsqcurvefit(fun,x0,ibins,nhB,lb,ub,options);
                mi=p(1); % max height
                lp=find(ibins>=p(2)); % find bin position corresponding to mean position
                ib=find(nhB(1:lp(1))<0.2*mi,1,'last');
                ie=find(nhB(lp(1):end)<0.2*mi,1,'first');
                tmp=ibins(ib:1:lp(1)+ie-1);
                
                %% Plot check
                %figure; plot(ibins,nhB,'.r'); hold on; title(c);
                %plot(ibins,fun(p,ibins),'-b');
                %return

                if isempty(tmp)==0
                    tmp2=find(iBeforeAct(pos)>=tmp(1) & iBeforeAct(pos)<=tmp(end) & iBeforeAct(pos)<=ibins(61)); % Also remove pixels that are past the max intensity of the microscope
                
                    %% Calculate Dilute phase mean values with filtered data using raw intensity histograms
                    rbins=500:500:max(niAfterAct);
                    hBR=hist(niBeforeAct(pos(tmp2)),rbins);
                    hAR=hist(niAfterAct(pos(tmp2)),rbins);
                    nhBR=hBR/sum(hBR);
                    nhAR=hAR/sum(hAR);

                    % Make sure before and after activation there are at least 5 bins, and at least 100 pixels to work with
                    if length(find(hBR~=0))>5 & length(find(hAR~=0))>5 & length(tmp2)>=100 
                        % Find dilute mean before activation
                        fun = @(p,rbins)p(1) * exp(-((rbins - p(2))/p(3)).^2);
                        x0=[0.6 exp(ibins(lp(1))) 300];
                        lb=[];
                        ub=[];
                        p = lsqcurvefit(fun,x0,rbins,nhBR,lb,ub);
                        %p = lsqcurvefit(fun,x0,rbins,nhBR,lb,ub,options);
                        diluteIntB(c)=p(2);

                        % Find dilute mean after activation
                        % Make sure dilute mean after activation can't be larger than before activation and that second peak occurs at high intensity
                        options2 = fitoptions('gauss2','Lower',[-Inf -Inf 0 -Inf 60000 0],'Upper',[Inf diluteIntB(c) Inf Inf 65535 Inf]);
                        [f gof]=fit(rbins',nhAR','gauss2',options2);
                        diluteIntA(c)=f.b1;
                        denseIntA(c)=f.b2;
                        mygof(c)=gof.rsquare;
                        %return
                        
                        % Check to make sure fits are good - cannot have a negative value and must have r2 gte 0.85
                        if diluteIntA(c)<0 | diluteIntB(c)<0 | gof.rsquare<0.85
                            diluteIntA(c)=nan;
                            diluteIntB(c)=nan;
                            denseIntA(c)=nan;
                        end
                        
                        % Plot check                        
                        %figure; plot(f,rbins,nhAR,'o'); hold on; title(c);
                        %plot(rbins,nhBR,'.r'); 
                        %plot(rbins,fun(p,rbins),'-b');
                        %return
				
                        clear x0; clear p; clear options2; clear f; 
                    end
                    
                clear rbins; clear hBR; clear hAR; clear nhBR; clear nhAR; clear tmp2; 
                end
                
                clear p; clear lp;  clear mi; clear ib; clear ie; clear tmp; 
            end
                
            clear pos; clear hB; clear hA; clear nhB; clear nhA; 
        end
        
        fullDiluteIntBfromMean=[fullDiluteIntBfromMean diluteIntBfromMean];
        fullDiluteIntBfromStd=[fullDiluteIntBfromStd diluteIntBfromStd];
        fullDiluteIntB=[fullDiluteIntB diluteIntB];
        fullDiluteIntA=[fullDiluteIntA diluteIntA];
        fullDenseIntA=[fullDenseIntA denseIntA];
        fullgof=[fullgof mygof];
        fullHighIntensityCells=[fullHighIntensityCells highIntCells];
        numCells(i)=length(ucels); 
        
        figure(274+d);
        x=0:1e4:4e4;
        y=0:1e4:4e4;
        plot(diluteIntB,diluteIntA,'o','color',mycolor(i,:)); hold on;
        for c=1:length(diluteIntB)
            text(diluteIntB(c),diluteIntA(c),num2str(c),'color',mycolor(i,:));
        end
        plot(x,y,'-k'); hold on;
        
        %return
        clear posc; clear xpos; clear ypos; 
        clear niBeforeAct; clear niAfterAct; clear iBeforeAct; clear iAfterAct;
        clear cels; clear ucels; clear bmat; clear amat;  
        clear diluteIntB; clear diluteIntA; clear denseIntA; clear mygof; 
        clear diluteIntBfromMean; clear diluteIntBfromStd; clear highIntCells; 
    end
    
    
    %return
    save(['../../data/optoDroplets_csat/' type{1} '/Mean_Dilute_Phase_Intensities_For_Phase_Separated_Cells_by_Gauss_Fit_' dir{d} '_' cloc{1} '.mat'],'fullDiluteIntB','fullDiluteIntA','fullDenseIntA','fullgof','numCells','fullDiluteIntBfromMean','fullDiluteIntBfromStd','fullHighIntensityCells');
    
    fclose('all');
    clear da; clear uimagenames; clear numImgs; 
    clear fullDiluteIntB; clear fullDiluteIntA; clear fullDenseIntA; clear numCells; clear fullgof; clear numCells; 
    clear fullDiluteIntBfromMean; clear fullDiluteIntBfromStd; clear fullHighIntensityCells; 
end
 
