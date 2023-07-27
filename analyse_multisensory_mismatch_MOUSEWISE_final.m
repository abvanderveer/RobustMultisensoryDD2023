clcl
%%%% this is the one for the paper!. subplot 2F or something
for file_load=1;

   files= {};

   mouseids1=[];
   sexs1=[];
end

limtrials1=15; %limit the number of trials analysed for each stimulus type(control/redundandt/deviant).
baseBEGIN=10; baseEND=14; %what you call the baseline
BEGINo=15; ENDo=28; %these are the analysis windows. note: stimulus starts at 14 and finishes at 28
onestim=0; %if 0, then plot responses to both stimuli (e.g. 45 and 135 deg). if 1, then plot only to first orientation. if 2, then plot to second

for file=1:size(files,1);
    load(strcat('CA_workspace_SORTED_',files{file}));
    %%%%%%%this compiles all the files into three matrices
    for buildingmatrices=1;
        singtrialsall=[];
        count=0; counttest=0;
        effectsall=[];  %this is the main varable
        actcont=[];
        countall=0;
    end
   
    
   load(strcat('CA_workspace_SORTEDdifferent_',files{file}),'vectrials')
            
    framerate=28;

    for limit_some_trials=1;
        if limit_some_trials==1; 
            if limtrials1>0;
                limtrials=limtrials1;
                for ss1=1:2;
                    vecstim=vectrials{ss1,2};
                    tmps1=[];
                    for tr1=1:size(vecstim,2);
                        
                        if and(sum(vecstim(tr1)==tmps1)<2,vecstim(tr1)~=999)
                            tmps1=horzcat(tmps1,vecstim(tr1));
                            vecstim(tr1)=99;
                        end
                        
                    end
                    for additionally_remove_first_few_controls=1;
%                         tmps1=[]; gh=size(vecstim,2)+1;
%                         for tr1=1:size(vecstim,2);
%                             if and(vecstim(gh-tr1)>15,vecstim(gh-tr1)<22)
%                                 if sum(tmps1==1)<(limtrials)
%                                     tmps1=horzcat(tmps1,1);
%                                 else
%                                     vecstim(gh-tr1)=99;
%                                 end
%                             end
%                         end
                    end
                    vectrials{ss1,2}=vecstim;
                end
            else
                for ss1=1:2;
                    vecstim=vectrials{ss1,2};
                    tmps1=[];
                    for tr1=1:size(vecstim,2);
                        
                        if and(sum(vecstim(tr1)==tmps1)==0,vecstim(tr1)~=15)
                            tmps1=horzcat(tmps1,vecstim(tr1));
                            vecstim(tr1)=99;
                        end
                        
                    end
                    vectrials{ss1,2}=vecstim;
                end
            end
        end
    end

    for analysis=1
        
        samp=1/framerate;
        timeaxis=-(samp*(framerate/2)):samp:(samp*framerate); 
        
        for stim=1:2;
            limtrials=limtrials1;
            vecdat=vectrials{stim,1};
            vecstim=vectrials{stim,2};
            if limtrials1>0;
                for checkingfortoo_few_trials=1;
                    ss=[1 2 3 4 5 15];
                    for gss=1:6;
                        sss(gss)=sum(vecstim==ss(gss));
                    end
                    if min(sss)<limtrials1;
                        limtrials=min(sss);
                    end
                    disp(limtrials);
                end
            
            effects=[];
            for z1=1:8
                tmp1=vecdat(:,:,vecstim==z1);
                if size(tmp1,3)>(limtrials-1)
                    effects(:,:,z1+1)=mean(tmp1(:,:,1:limtrials),3);
                else
                    effects(:,:,z1+1)=mean(tmp1(:,:,1:end),3);
                end
            end
            tmp1=vecdat(:,:,vecstim==15);
            if size(tmp1,3)>(limtrials-1)
                effects(:,:,10)=mean(tmp1(:,:,1:limtrials),3);
            else
                effects(:,:,10)=mean(tmp1(:,:,1:end),3);
            end
            tmp1=vecdat(:,:,vecstim==phzstd(stim));
            
                if size(tmp1,3)>(limtrials-1)
                    effects(:,:,1)=mean(tmp1(:,:,1:limtrials),3);
                     else
                    effects(:,:,1)=mean(tmp1(:,:,1:end),3);
                end
                
            else
                 vecdat=vectrials{stim,1};
                vecstim=vectrials{stim,2};
                effects=[];
                stimax=[phzstd(stim) 1 2 3 4 5 6 7 8 15];
                for z1=1:10
                    effects(:,:,z1)=mean(vecdat(:,:,vecstim==stimax(z1)),3);
                end
            end
                
            
            for c1=1:size(effects,1) %baseline correction
                basos=[1 2 2 2 2 2 2 2 2 2; 1 2 2 2 2 2 2 2 2 2];
                at=squeeze(effects(c1,baseBEGIN:baseEND,:));
                
                for z1=1:size(effects,3)
                    mnn=mean(effects(c1,baseBEGIN:baseEND,z1),2);
                     stt=std(at(:));
                    
                    if stt<.001; stt=std(effects(c1,ENDo:end,z1)); end
                    if stt<.001; stt=10000; end
                    for t1=1:size(effects,2)
                        effects(c1,t1,z1)=(effects(c1,t1,z1)-mnn)./stt;
                    end
                    mnns(file,stim,c1,z1)=mnn;
                    stts(file,stim,c1,z1)=stt;
                end
                
            end
            
            
            for getspontaneousactivity=1
                if stim==1
                    for c1=1:size(effects,1)
                        tmp=[];
                        for ph=1:3
                            
                            tmp2=dfofanalyse(c1,phases==ph);
                            tmp3=stimvisMAIN(1,phases==ph);
                            t1=0; t2=0; %t1 becomes teh beginning of stims. t2 becomes the end.
                            for t=100:size(tmp3,2)-100
                                if and(t1==0,tmp3(1,t)>0)
                                    t1=t;
                                end
                                if and(t2==0,t1>0)
                                    if max(tmp3(1,t-99:t+99))==0
                                        t2=t-99;
                                    end
                                end
                            end
                            if t2==0; t2=size(tmp3,2)-99; end
                            
                          tmp22=tmp2(1,tmp2(1,:)>0); pr=prctile(tmp22,95); tmp2=tmp2/std(tmp22<pr);
                            
                            for alltimepoints=0 %if you want to just analyse interstimulus activity, make this =0
                                if alltimepoints==1
                                    bin1=floor((t2-t1)/10);
                                    for b=1:10;
                                        tmp=horzcat(tmp,mean(tmp2(1,((b-1)*bin1)+t1+1:(b*bin1)+t1)));
                                    end
                                else
                                    tmp2=tmp2(1,t1:t2); tmp3=tmp3(1,t1:t2);
                                    tmp2=tmp2(1,tmp3==0);
                                    bin1=floor(size(tmp2,2)/10);
                                    for b=1:10
                                        tmp=horzcat(tmp,mean(tmp2(1,((b-1)*bin1)+1:(b*bin1))));
                                    end
                                end
                            end
                            
                        end
                        actcont=vertcat(actcont,tmp);
                    end
                end
            end
            
            %for getting single trial data. make sure to comment out the
            %"exclude first few trials" section above if you want to really see the truth!
            tmpcells=[];
            for cc=1:size(vecdat,1);
                
                tmp00=dfofanalyse(cc,30:449)./stts(file,stim,cc,1); 
                tmp1=squeeze(vecdat(cc,:,vecstim==phzstd(stim)))./stts(file,stim,cc,1);
                tmp0=[]; stim0=[];
                for g=1:limtrials;
                    tmp0=vertcat(tmp0,tmp1(:,g));
                    stim0=horzcat(stim0,zeros(1,13),ones(1,14),zeros(1,16));
                end
                stim0=horzcat(zeros(1,420),stim0); tmp0=vertcat(tmp00',tmp0); 
                tmp1=squeeze(vecdat(cc,:,vecstim==1))./stts(file,stim,cc,2); ; tmp2=squeeze(vecdat(cc,:,vecstim==2))./stts(file,stim,cc,3);
                tmp3=squeeze(vecdat(cc,:,vecstim==3))./stts(file,stim,cc,4);  tmp4=squeeze(vecdat(cc,:,vecstim==4))./stts(file,stim,cc,5); 
                tmp5=squeeze(vecdat(cc,:,vecstim==5))./stts(file,stim,cc,6);  tmp6=squeeze(vecdat(cc,:,vecstim==6))./stts(file,stim,cc,7); 
                tmp7=squeeze(vecdat(cc,:,vecstim==7))./stts(file,stim,cc,8);  tmp8=squeeze(vecdat(cc,:,vecstim==8))./stts(file,stim,cc,9); 
                tmp15=squeeze(vecdat(cc,:,vecstim==15))./stts(file,stim,cc,10);  ggg=0;
                for g=1:5;
                    try
                        tmp0=vertcat(tmp0,tmp1(:,g),tmp2(:,g),tmp3(:,g),tmp4(:,g),tmp5(:,g),tmp6(:,g),tmp7(:,g),tmp8(:,g),tmp15(:,g));
                        
                        stim0=horzcat(stim0,zeros(1,13),ones(1,14)-.5,zeros(1,16),zeros(1,13),ones(1,14)-.5,zeros(1,16),zeros(1,13),ones(1,14)-.5,zeros(1,16),zeros(1,13),ones(1,14)-.5,zeros(1,16),zeros(1,13),ones(1,14)-.5,zeros(1,16),zeros(1,13),ones(1,14)-.5,zeros(1,16),zeros(1,13),ones(1,14)-.5,zeros(1,16),zeros(1,13),ones(1,14)-.5,zeros(1,16),zeros(1,13),ones(1,14)+.5,zeros(1,16));
                        ggg=g;
                    end
                end
                if ggg==5;
                tmpcells=vertcat(tmpcells,tmp0');
                stim10=stim0;
                end
            end
            
            
            
            if onestim==0;
                effectsall((1+countall):(countall+size(effects,1)),1:size(effects,2),2:10,stim)=effects(:,1:end,2:10);
                effectsall((1+countall):(countall+size(effects,1)),1:size(effects,2),1,stim)=effects(:,1:end,1);
                
                singtrialsall((1+countall):(countall+size(effects,1)),1:size(tmpcells,2),stim)=tmpcells;
            elseif onestim==1;
                if stim==1;
                    effectsall((1+countall):(countall+size(effects,1)),1:size(effects,2),2:10,stim)=effects(:,1:end,2:10);
                    effectsall((1+countall):(countall+size(effects,1)),1:size(effects,2),1,stim)=effects(:,1:end,1);
                    
                    singtrialsall((1+countall):(countall+size(effects,1)),1:size(tmpcells,2),stim)=tmpcells;
                else
                    effectsall((1+countall):(countall+size(effects,1)),1:size(effects,2),2:10,stim)=effects(:,1:end,2:10)-effects(:,1:end,2:10);
                    effectsall((1+countall):(countall+size(effects,1)),1:size(effects,2),1,stim)=effects(:,1:end,1)-effects(:,1:end,1);
                    singtrialsall((1+countall):(countall+size(effects,1)),1:size(tmpcells,2),stim)=tmpcells-tmpcells;
                end
            elseif onestim==2;
                if stim==2;
                    effectsall((1+countall):(countall+size(effects,1)),1:size(effects,2),2:10,stim)=effects(:,1:end,2:10);
                    effectsall((1+countall):(countall+size(effects,1)),1:size(effects,2),1,stim)=effects(:,1:end,1);
                    singtrialsall((1+countall):(countall+size(effects,1)),1:size(tmpcells,2),stim)=tmpcells;
                else
                    effectsall((1+countall):(countall+size(effects,1)),1:size(effects,2),2:10,stim)=effects(:,1:end,2:10)-effects(:,1:end,2:10);
                    effectsall((1+countall):(countall+size(effects,1)),1:size(effects,2),1,stim)=effects(:,1:end,1)-effects(:,1:end,1);
                    singtrialsall((1+countall):(countall+size(effects,1)),1:size(tmpcells,2),stim)=tmpcells-tmpcells;
                end
            end
            
        end
        
        countall=countall+size(effects,1);
    end
    
       disp(size(vecdat,2));


    
titles={'Training';' ';' ';' ';'controls';' ';' ';' ';' ';'DEVIANT'};

if size(effectsall,1)<15;
    include_cutoff=0; %number of standard deviations above baseline for inclusion
else
    include_cutoff=1;
end
redundantnumberA=0;
exclude_outlier=100; %get rid of the cell if any responses are this big

for stringentTEST=0;
        %%%%%%%%%figures out relative responses to different stimuli
        simpleeffects=[]; c1=1; c2=1; cols=[]; newindices=[];notindices=[];
        for stim=1:2
            for c=1:size(effectsall,1)
                if redundantnumberA==0; 
                	r=max(mean(effectsall(c,BEGINo:ENDo,2:6,stim),2));
                else
                    r=mean(effectsall(c,BEGINo:ENDo,redundantnumberA+1,stim));
                end
                cn=mean(effectsall(c,BEGINo:ENDo,1,stim));
                d=mean(effectsall(c,BEGINo:ENDo,10,stim));cn=d; 
                %d=cn; r=cn;
                if and(or(or(r>include_cutoff,cn>include_cutoff),d>include_cutoff),and(and(r<exclude_outlier,cn<exclude_outlier),d<exclude_outlier))
                    simpleeffects(c1,1)=(cn);
                    simpleeffects(c1,2)=(r);
                    simpleeffects(c1,3)=(d);
                    newindices(c1,1)=c;
                    newindices(c1,2)=stim;
                    c1=1+c1;
                else
                    notindices(c2,1)=c;
                    notindices(c2,2)=stim;
                    c2=1+c2;
                    
                end
                
            end
        end
        for makingcolors=1;
            for c=1:(c1-1)
                s=simpleeffects(c,1)./max(simpleeffects(c,:));
                cols(c,2)=s;
                s=simpleeffects(c,3)./max(simpleeffects(c,:));
                cols(c,1)=s;
                s=simpleeffects(c,2)./max(simpleeffects(c,:));
                cols(c,3)=s;
            end
            for c=1:(c1-1);
                a=cols(c,:)-min(cols(c,:));
                cols(c,:)=cols(c,:)./max(cols(c,:));
            end
        end
   
end

for countcellstotal=1; % for only focusing on one response per cell
    numcells=size(unique(newindices(:,1)),1);
    disp(strcat('all cells=',num2str(size(effectsall,1)),'; responsive=',num2str(numcells),'; percent=',num2str(numcells/size(effectsall,1))));
    ex1=zeros(size(simpleeffects,1),1);
    for c=1:size(effectsall,1);
        if sum(newindices(:,1)==c)>1;
            clear tmp1 tmp2; cnt1=1;
            for c1=1:size(newindices,1);
                if newindices(c1,1)==c;
                    tmp1(cnt1,1)=mean(simpleeffects(c1,[1 2 3]));
                    tmp2(cnt1,1)=c1;
                    cnt1=1+cnt1;
                end
            end
            if tmp1(1,1)>tmp1(2,1);
                ex1(tmp2(2,1),1)=1;
            else
                ex1(tmp2(1,1),1)=1;
            end
        end
    end
    for onlylookatonePERcell=1;
        if onlylookatonePERcell==1;
            simpleeffects=simpleeffects(ex1==0,:);
            newindices=newindices(ex1==0,:);
        end
    end
    for makingcolors=1;
    c1=size(simpleeffects,1)+1; clear cols
    for c=1:(c1-1)
        s=simpleeffects(c,1)./max(simpleeffects(c,:));
        cols(c,2)=s;
        s=simpleeffects(c,3)./max(simpleeffects(c,:));
        cols(c,1)=s;
        s=simpleeffects(c,2)./max(simpleeffects(c,:));
        cols(c,3)=s;
        
    end
    for c=1:(c1-1);
        a=cols(c,:)-min(cols(c,:));
        cols(c,:)=cols(c,:)./max(cols(c,:));
        cols(c,2)=0;cols(c,3)=0;%cols(c,1)=0;
    end
end

end


redundantnumber=(redundantnumberA+1); %if you want to plot (but not "select by") a redudnant later in the stream, use this line
for plottingeverything=1; 
%%%%for selecting all or subset of cells for subesequent plots
%%%for selecting all or subset of cells for subesequent plots
plotanti=0; %if =1, then plot everything outside of the selected box
selections=1; %if the cells you want to select are scattered more than a single box can capture, make this 2 or 3 and select multiple regions.
pind=1:size(simpleeffects,1); pind=pind';
close all

    for plotbOLTH=1;
        effectsallBOTH=effectsall; 
       tmpeff=[];
    for c=1:numcells;
        tmpeff(c,:,:)=effectsallBOTH(newindices(pind(c),1),:,:,(newindices(pind(c),2)));
        for t=1:size(effectsallBOTH,2)
            mn=mean(mean(effectsallBOTH(newindices(pind(c),1),t,[1 redundantnumber 10],newindices(pind(c),2)),3),4);
            a1(c,t)=mean(mean(effectsallBOTH(newindices(pind(c),1),t,redundantnumber,newindices(pind(c),2)),3),4);%-mn;
            a2(c,t)=mean(effectsallBOTH(newindices(pind(c),1),t,1,newindices(pind(c),2)),4);%-mn;
            a3(c,t)=mean(effectsallBOTH(newindices(pind(c),1),t,10,newindices(pind(c),2)),4);%-mn;
        end
    end


        x=timeaxis;
        y1=mean(mean(mean(tmpeff(:,:,redundantnumber),3),4),1); y1a=y1;%-mean(y1(1,BEGINo-5:BEGINo-1));
        err1=std(a1,0,1)./sqrt(numcells);

        y2=mean(mean(mean(tmpeff(:,:,1),3),4),1);y2a=y2;%-mean(y2(1,BEGINo-5:BEGINo-1));
        err2=std(a2,0,1)./sqrt(numcells);

        y3=mean(mean(mean(tmpeff(:,:,10),3),4),1); y3a=y3;%-mean(y3(1,BEGINo-5:BEGINo-1));
        err3=std(a3,0,1)./sqrt(numcells);
        if numcells>0;
          % figure; shadedErrorBar(x,y1a,err1,'lineProps','k'); hold on;shadedErrorBar(x,y2a,err1,'lineprops','k','patchSaturation',0);hold on;shadedErrorBar(x,y3a,err1,'lineProps','r');
       %  xlim([-.2 1]); title ('Average of cell responses', 'FontSize', 16, 'FontWeight', 'bold');ylabel('spike rate','FontSize',16,'FontWeight','bold'); xlabel('time (sec)','FontSize',16,'FontWeight','bold'); set(gca,'FontSize',16,'FontWeight','bold');
         %   set(gcf,'Color','w');

            
            for statstat=1;
                r1=mean(mean(tmpeff(:,BEGINo:ENDo,redundantnumber),3),2);%-mean(mean(tmpeff(:,BEGINo-5:BEGINo-1,redundantnumber),3),2);
                c1=mean(mean(tmpeff(:,BEGINo:ENDo,1),3),2);%-mean(mean(tmpeff(:,BEGINo-5:BEGINo-1,1),3),2);
                d1=mean(mean(tmpeff(:,BEGINo:ENDo,10),3),2);%-mean(mean(tmpeff(:,BEGINo-5:BEGINo-1,10),3),2);

                %[h,p1,ci,stats1]=ttest(r1,c1);
                [h,p2,ci,stats2]=ttest(d1,r1);
            end
            yy=ylim;figure;
            errorbar_groups([mean(c1) mean(r1) mean(d1)]',[std(c1)/sqrt(size(c1,1)-1) std(r1)/sqrt(size(r1,1)-1) std(d1)/sqrt(size(d1,1)-1)]', 'bar_colors',[1 1 1; .5 .5 .5; 1 0 0]);
            ylabel('spike rate','FontSize',16,'FontWeight','bold'); xlabel('time (sec)','FontSize',16,'FontWeight','bold'); set(gca,'FontSize',16,'FontWeight','bold');
            set(gcf,'Color','w');
            title(strcat('T_D_D','=',num2str(round(stats2.tstat*100)/100),', p=',num2str(round(p2*1000)/1000)));
        else
            figure; plot(x,y1,'b','LineWidth',3); hold on;plot(x,y2,'k','LineWidth',3); hold on;plot(x,y3,'r','LineWidth',3);yy=ylim;
            xlim([-.2 1]); title ('Average of cell responses', 'FontSize', 16, 'FontWeight', 'bold');
        end
        close all force
            mouseavs1(file,1)=mean(r1); disp(size(r1,1));
            mouseavs1(file,2)=mean(d1); disp(size(d1,1));
        for condol=1:10;
            for c=1:size(tmpeff,1);
                tmpeff(c,:,condol)=tmpeff(c,:,condol)-mean(tmpeff(c,BEGINo-10:BEGINo-3,condol));
            end
            mouseavs_time1(file,:,condol)=mean(tmpeff(:,:,condol),1);%-mean(min(tmpeff(:,BEGINo-10:BEGINo-1,condol),[],2),1);
        end
        responsivecells(file,1)=numcells;
       

       
    end

   
clear r1 d1 effects tmpeff
end
end

inc=responsivecells>0;
mouseavs=mouseavs1(inc,:);
mouseavs_time=mouseavs_time1(inc,:,:);
mouseids=mouseids1(inc);
suba1=0;
for suba=1:max(mouseids);
    if sum(mouseids==suba)>0
        suba1=suba1+1; 
        newdat(suba1,:)=mean(mouseavs(mouseids==suba,:),1);
        newdat_time(suba1,:,:)=mean(mouseavs_time(mouseids==suba,:,:),1);
    end
end

conds=[2 2 2; 3 3 3; 4 4 4; 5 5 5; 6 6 6; 7 8 9; 10 10 10];
        for manyplotthing=1;
            x=timeaxis;
            figure;
            for cnd=1:7;
                y1=mean(mean(newdat_time(:,:,conds(cnd,:)),1),3);
                y1=y1-min(y1(1,baseBEGIN:baseEND));
              
                    err=std(newdat_time(:,:,conds(cnd,1)),0,1)./sqrt(size(newdat_time(:,:,conds(cnd,1)),1));

                subplot(1,7,cnd); shadedErrorBar(x,y1,err,'lineProps','k'); set(gcf,'Color','w'); xlim([-.2 .8]);
                if cnd==1;
                    ylabel('spike rate','FontSize',16,'FontWeight','bold'); xlabel('time (sec)','FontSize',16,'FontWeight','bold'); set(gca,'FontSize',16,'FontWeight','bold');
                    title(titles{cnd}, 'FontSize', 18, 'FontWeight', 'bold');
                else
                    title(titles{cnd}, 'FontSize', 18, 'FontWeight', 'bold'); yticks([]);set(gca,'FontSize',16,'FontWeight','bold');
                end
            end
            for fixaxes=1

                ymin1=1; ymax1=-1;
                for cnd=1:7
                    subplot(1,7,cnd); yy=ylim;
                    if yy(1)<ymin1; ymin1=yy(1); end
                    if yy(2)>ymax1; ymax1=yy(2); end

                end
                for cnd=1:7
                    subplot(1,7,cnd); ylim([ymin1 ymax1]); 
                end

            end
        end
  
make_eps_savable
        for manyplotthing=1;
            clear tvals pvals
            for cnd=1:7;
                for cnd2=1:7;
                    if cnd~=cnd2;
                        y1=mean(mean(newdat_time(:,BEGINo:ENDo,conds(cnd,:)),2),3)-mean(mean(newdat_time(:,baseBEGIN:baseEND,conds(cnd,:)),2),3);
                        y2=mean(mean(newdat_time(:,BEGINo:ENDo,conds(cnd2,:)),2),3)-mean(mean(newdat_time(:,baseBEGIN:baseEND,conds(cnd2,:)),2),3);;
                        [h,p,ci,stats]=ttest(y1,y2);
                        tvals(cnd,cnd2)=stats.tstat;
                        pvals(cnd,cnd2)=1-p;
                        if cnd2==7;
                            cohD(cnd)=(mean(y2)-mean(y1))./std(vertcat(y1,y2));
                        end
                        if cnd==4;
                            mouseyav(:,1)=y1';
                        end
                        if cnd==7;
                            mouseyav(:,2)=y1';
                        end
                    end
                end
            end
        end
        figure; 
        subplot(1,2,1); imagesc(tvals); subplot(1,2,2); imagesc(pvals);
        figure; plot(cohD); make_eps_savable

               suba1
       
figure; scatter((sexs1*0)+1,mouseyav(:,2)',100,'r','fill'); xlim([-.5 1.5]);hold on; scatter((sexs1*0)+0,mouseyav(:,1)',100,'k','fill');
for x=1:suba1; line([0 1],[mouseyav(:,1) mouseyav(:,2)],'Color','k'); hold on; end;
make_eps_savable
