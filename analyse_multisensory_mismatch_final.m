clcl
%%%% this is the one for the paper!
%%used first file 2021_0222green for demo
for file_load=1;
files= {;};
end

limtrials1=15; %limit the number of trials analysed for each stimulus type(control/redundandt/deviant). 15 for paper
baseBEGIN=10; baseEND=14; %what you call the baseline. 10 to 14 for paper
BEGINo=15; ENDo=28; %these are the analysis windows. note: stimulus starts at 15 and finishes at 28. used that for paper.
onestim=0; 

%%%%%%%this compiles all the files into three matrices
for buildingmatrices=1;
      singtrialsall=[];
    count=0; counttest=0;
    effectsall=[];  %this is the main varable
    effectsallSTD=[];  
    actcont=[];
    countall=0;
end
for file=1:size(files,1);
    load(strcat('CA_workspace_SORTED_',files{file}));


    

    for different_df=1;
        oasis_setup;
        try load(strcat('CA_workspace_SORTEDdifferent_',files{file}),'vectrials')
        catch
            for calculateDFOF=1;
                framerate=28;
                dfofall=data;
                deltaf=data-data;
                clear data
                for c=1:size(dfofall,1);

                    % Y{1}=dfofall(c,:);
                    % [ca_foopsi,cb,c1,~,~,spikes_foopsi] = constrained_foopsi(y,[],[],g2);
                    % SAMPLES = cont_ca_sampler(dfofall(:,:)')

                    [cx, sx, options] = deconvolveCa(dfofall(c,:),'constrained_foopsi', 'ar1');% 'optimize_pars', true, 'optimize_b', true)
                    deltaf(c,:)=smoothdata(sx,'gaussian',5);
                end
                dfofanalyse2=deltaf;
                dfofanalyse=deltaf;
                % for c=1:size(deltaf,1); for t=1:size(deltaf,2); if deltaf(c,t)<0; deltaf(c,t)=0; end; end; end



            end

            for phz=1:2;
                stimvistypeall=stimvistypeMAIN;
                for vectorize_populations=1;
                    for stimulus=1:2;
                        vecdat=[];
                        vecstim=[];
                        tr=0;
                        %controls
                        tmpvis=stimvistypeall;
                        for t=1:size(tmpvis,2); if and(tmpvis(1,t)~=0,tmpvis(1,t)~=stimulus); tmpvis(1,t)=0; end; end;
                        tmpvis=tmpvis(phases==1); dfoftmp=dfofanalyse(:,phases==1);
                        for t=framerate:size(dfoftmp,2)-(framerate*2);
                            if (tmpvis(1,t)-tmpvis(1,t-1))>.5;
                                tr=1+tr;
                                vecdat(:,:,tr)=dfoftmp(:,t-(framerate*1):t+framerate*2);
                                vecstim(tr)=(tmpvis(1,t)+19);
                            end
                        end

                        %reduds and devios
                        tmpvis=stimvistypeall;
                        for t=1:size(tmpvis,2); if and(tmpvis(1,t)~=0,and(tmpvis(1,t)~=stimulus,tmpvis(1,t)~=(stimulus+2))); tmpvis(1,t)=0; end; end;
                        tmpvis=tmpvis(phases==2); dfoftmp=dfofanalyse(:,phases==2);
                        cnts=0;
                        for t=framerate:size(dfoftmp,2)-(framerate*2);
                            if (tmpvis(1,t)-tmpvis(1,t-1))>.5
                                tr=1+tr;
                                if (tmpvis(1,t)-tmpvis(1,t-1))<2.5;
                                    cnts=1+cnts;
                                    vecdat(:,:,tr)=dfoftmp(:,t-(framerate*1):t+framerate*2);
                                    if cnts>14; cnts=14; end;
                                    vecstim(tr)=(cnts);
                                else
                                    cnts=0;
                                    vecdat(:,:,tr)=dfoftmp(:,t-(framerate*1):t+framerate*2);
                                    vecstim(tr)=(15);
                                end

                            end
                        end

                        vectrials{stimulus,1}=vecdat;
                        vectrials{stimulus,2}=vecstim;
                    end



                end
            end





            for fixframerate=1;
                if framerate==56;
                    vectrials1=vectrials;
                    disp(strcat('file ',num2str(file),'has wonky FR'));
                    for stimx=1:2;
                        vecdat1=vectrials1{stimx,1}; vecdat2=[];
                        for cx=1:size(vecdat1,1);
                            for trx=1:size(vecdat1,3);
                                vecdat2(cx,:,trx)=vecdat1(cx,28:98,trx);
                            end
                        end
                        vectrials{stimx,1}=vecdat2;
                    end
                    clear vectrials1 vecdat1 vecdat2
                    framerate=28;
                end

            end

            for rejitter=1;
                dohaloscore=1;
                if rejitter==1;
                    try
                        load(strcat('CA_rejitter2_',files{file}),'x');
                        disp(strcat('CA_rejitter2_',files{file},' jit',num2str(x),'sample'))

                    catch
                        vec1=vectrials{1,1}; st1=vectrials{1,2};
                        celltmp=[];
                        for c=1:size(vec1,1);
                            celltmp(c,:)=mean(vec1(c,:,st1<16),3);
                            celltmp(c,:)=celltmp(c,:)-min(celltmp(c,16:27));
                        end
                        tmp=mean(celltmp(:,28:42),2);
                        [ix,b]=sort(tmp);
                        cut=round(size(tmp,1)./3);
                        celltmp1=celltmp(b(cut:end),:);
                        celltmp2=celltmp1(1,:)-celltmp1(1,:);
                        for t=(framerate*2):(size(stimvistypeMAIN,2)-(framerate*2));
                            if and(stimvistypeMAIN(1,t)>2.1,stimvistypeMAIN(1,t-1)==0);
                                celltmp2=celltmp2+mean(halos(:,t-(framerate*1):t+framerate*2),1);
                            end
                        end
                        if dohaloscore==1; celltmp1=celltmp2; end;

                        figure; plot(mean(celltmp1,1));title('select stim begin time.. first frame of rampup'); [x,y]=ginput(1);

                        x=round(x);
                        if or(x<0,x>141);
                            x=28;
                        end
                        save(strcat('CA_rejitter2_',files{file}),'x');

                    end
                    for stimx=1:2;
                        vecdata=vectrials{stimx,1};vecdata1=[];
                        for c=1:size(vecdata,1);
                            for tr=1:size(vecdata,3);
                                vecdata1(c,:,tr)=vecdata(c,x-14:x+28,tr);
                            end
                        end
                        vectrials{stimx,1}=vecdata1;
                    end
                else
                    x=28;
                    for stimx=1:2;
                        vecdata=vectrials{stimx,1};vecdata1=[];
                        for c=1:size(vecdata,1);
                            for tr=1:size(vecdata,3);
                                vecdata1(c,:,tr)=vecdata(c,x-14:x+28,tr);
                            end
                        end
                        vectrials{stimx,1}=vecdata1;
                    end
                end
                close all force
            end



            save(strcat('CA_workspace_SORTEDdifferent_',files{file}),'vectrials');

        end


    end

    
    

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
    
    for flippymixymatch=1;
        if flippymixymatch==1;
            vectrialsalt=vectrials;
            clear vectrials
               
                vectmp1=vectrialsalt{1,1};
                stimtmp1=vectrialsalt{1,2};
                vectmp2=vectrialsalt{2,1};
                stimtmp2=vectrialsalt{2,2};
                aa1=0; bb1=0; clear tmpdevs1 tmpalls1 tmpalls1stim tmpdevs2 tmpalls2 tmpalls2stim 
                for t=1:size(stimtmp1,2);
                    if stimtmp1(1,t)==15;
                        aa1=1+aa1;
                        tmpdevs1(:,:,aa1)=vectmp1(:,:,t); 
                    else
                        bb1=1+bb1;
                        tmpalls1(:,:,bb1)=vectmp1(:,:,t);
                        tmpalls1stim(1,bb1)=stimtmp1(1,t);
                    end
                end

                aa2=0; bb2=0;
                for t=1:size(stimtmp2,2);
                    if stimtmp2(1,t)==15;
                        aa2=1+aa2;
                        tmpdevs2(:,:,aa2)=vectmp2(:,:,t); 
                    else
                        bb2=1+bb2;
                        tmpalls2(:,:,bb2)=vectmp2(:,:,t);
                        tmpalls2stim(1,bb2)=stimtmp2(1,t);
                    end
                end

                tmpalls1stim(1,(bb1+1):(bb1+aa2))=zeros(1,aa2)+15;
                tmpalls1(:,:,(bb1+1):(bb1+aa2))=tmpdevs2;
                tmpalls2stim(1,(bb2+1):(bb2+aa1))=zeros(1,aa1)+15;
                tmpalls2(:,:,(bb2+1):(bb2+aa1))=tmpdevs1;

                vectrials{1,1}=tmpalls1; 
                vectrials{1,2}=tmpalls1stim;
                  vectrials{2,1}=tmpalls2; 
                vectrials{2,2}=tmpalls2stim;
        end
    end

    for analysis=1
        
        samp=1/framerate;
        timeaxis=-(samp*(framerate/2)):samp:(samp*framerate); timeaxis=timeaxis+.035; 
        
        for stim=1:2;
            limtrials=limtrials1;
            vecdat=vectrials{stim,1};
            vecstim=vectrials{stim,2};
            if limtrials1>0; 
            for checkingfortoo_few_trials=1;
                ss=[1 2 3 4 5 15 phzstd(stim)];
                for gss=1:7; 
                    sss(gss)=sum(vecstim==ss(gss));
                end
                if min(sss)<limtrials1;
                   limtrials=min(sss); 
                end
                disp(limtrials); 
            end
            
            effects=[]; effectsSTD=[];
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
                at=squeeze(effects(c1,baseBEGIN:baseEND,:)); %std from avg baseline
                
                for z1=1:size(effects,3)
                    mnn=mean(effects(c1,baseBEGIN:baseEND,z1),2); 
                    stt=std(at(:)); 
                    
                    if stt<.001; stt=std(effects(c1,ENDo:end,z1)); end
                    if stt<.001; stt=10000; end
                    for t1=1:size(effects,2)
                        effects(c1,t1,z1)=(effects(c1,t1,z1)-mnn);
                        effectsSTD(c1,t1,z1)=(effects(c1,t1,z1))./stt;
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
                effectsallSTD((1+countall):(countall+size(effectsSTD,1)),1:size(effectsSTD,2),2:10,stim)=effectsSTD(:,1:end,2:10);
                effectsallSTD((1+countall):(countall+size(effectsSTD,1)),1:size(effectsSTD,2),1,stim)=effectsSTD(:,1:end,1);

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
end

    
titles={'Training';' ';' ';' ';'controls';' ';' ';' ';' ';'DEVIANT'};
include_cutoff=1; %number of standard deviations above baseline for inclusion; 1 for paper
exclude_outlier=100; %get rid of the cell if any responses are this big
viewSTD=0; %zero for lineplots and stats, 1 for raster

simpleeffects=[]; c1=1; c2=1; cols=[]; newindices=[];notindices=[]; effectsallBOTH=[];  effectsallBOTHstd=[];
for c=1:size(effectsallSTD,1)
    r1=mean(mean(effectsallSTD(c,BEGINo:ENDo,2:4,1)));
    cn1=mean(effectsallSTD(c,BEGINo:ENDo,1,1));
    d1=mean(effectsallSTD(c,BEGINo:ENDo,10,1));cn1=d1;

    r2=mean(mean(effectsallSTD(c,BEGINo:ENDo,2:4,2)));
    cn2=mean(effectsallSTD(c,BEGINo:ENDo,1,2));
    d2=mean(effectsallSTD(c,BEGINo:ENDo,10,2));cn2=d2;

    if1=and(or(or(r1>include_cutoff,cn1>include_cutoff),d1>include_cutoff),and(and(r1<exclude_outlier,cn1<exclude_outlier),d1<exclude_outlier));
    if2=and(or(or(r2>include_cutoff,cn2>include_cutoff),d2>include_cutoff),and(and(r2<exclude_outlier,cn2<exclude_outlier),d2<exclude_outlier));
    if and(if1==1,if2==1);
        simpleeffects(c1,1)=(cn1+cn2)/2;
        simpleeffects(c1,2)=(r1+r2)/2;
        simpleeffects(c1,3)=(d1+d2)/2;
        effectsallBOTH(c1,:,:)=mean(effectsall(c,:,:,:),4);
        effectsallBOTHstd(c1,:,:)=mean(effectsallSTD(c,:,:,:),4);
        c1=1+c1;
    elseif and(if1==1,if2==0);
        simpleeffects(c1,1)=(cn1+cn1)/2;
        simpleeffects(c1,2)=(r1+r1)/2;
        simpleeffects(c1,3)=(d1+d1)/2;
        effectsallBOTH(c1,:,:)=mean(effectsall(c,:,:,1),4);
        effectsallBOTHstd(c1,:,:)=mean(effectsallSTD(c,:,:,1),4);
        c1=1+c1;
    elseif and(if1==0,if2==1);
        simpleeffects(c1,1)=(cn2+cn2)/2;
        simpleeffects(c1,2)=(r2+r2)/2;
        simpleeffects(c1,3)=(d2+d2)/2;
        effectsallBOTH(c1,:,:)=mean(effectsall(c,:,:,2),4);
        effectsallBOTHstd(c1,:,:)=mean(effectsallSTD(c,:,:,2),4);
        c1=1+c1;
    end
end
numcells=c1-1;
disp(strcat('all cells=',num2str(size(effectsall,1)),'; responsive=',num2str(numcells),'; percent=',num2str(numcells/size(effectsall,1))));

redundantnumber=4;% %plotting by 3rd redundant (i.e. 4) for the paper
for plottingeverything=1;

    close all

    for plotbOLTH=1;
        if viewSTD==1; effectsallBOTH=effectsallBOTHstd; end;
        tmpeff=[];
        for c=1:numcells;
            tmpeff(c,:,:)=effectsallBOTH(c,:,:);
            for t=1:size(effectsallBOTH,2)
                mn=mean(mean(effectsallBOTH(c,t,[1 redundantnumber 10]),3),4);
                a1(c,t)=mean(effectsallBOTH(c,t,redundantnumber),3);%-mn;
                a2(c,t)=(effectsallBOTH(c,t,1));%-mn;
                a3(c,t)=(effectsallBOTH(c,t,10));%-mn;
            end
        end


        x=timeaxis;
        y1=mean(mean(mean(tmpeff(:,:,redundantnumber),3),4),1); y1b=min(y1(1,BEGINo-5:BEGINo-1));
        err1=std(a1,0,1)./sqrt(numcells);

        y2=mean(mean(mean(tmpeff(:,:,1),3),4),1);y2b=min(y2(1,BEGINo-5:BEGINo-1));
        err2=std(a2,0,1)./sqrt(numcells);

        y3=mean(mean(mean(tmpeff(:,:,10),3),4),1); y3b=min(y3(1,BEGINo-5:BEGINo-1));
        err3=std(a3,0,1)./sqrt(numcells);

        y1a=y1-((y1b+y3b)/2);
        y2a=y2-((y1b+y3b)/2);
        y3a=y3-((y1b+y3b)/2);

        if numcells>1;
            figure; shadedErrorBar(x,y1a,err1,'lineProps','k'); hold on;shadedErrorBar(x,y3a,err1,'lineProps','r');%hold on;shadedErrorBar(x,y2a,err1,'lineprops','k','patchSaturation',0);
            xlim([-.21 1]); title ('Average of cell responses', 'FontSize', 16, 'FontWeight', 'bold');ylabel('z-scores','FontSize',16,'FontWeight','bold'); xlabel('time (sec)','FontSize',16,'FontWeight','bold'); set(gca,'FontSize',16,'FontWeight','bold');
            set(gcf,'Color','w');
            clear pvals_tmp aa1 bb1
            for t=1:size(tmpeff,2);
                y1=mean(mean(mean(tmpeff(:,:,redundantnumber),3),4),1); %y1a=y1-mean(y1(1,BEGINo-10:BEGINo-1));
                aa1=tmpeff(:,t,redundantnumber)-((y1b+y1b)/2);%-mean(y1(1,BEGINo-10:BEGINo-1));
                y3=mean(mean(mean(tmpeff(:,:,10),3),4),1); %y3a=y3-mean(y3(1,BEGINo-10:BEGINo-1));
                bb1=tmpeff(:,t,10)-((y2b+y3b)/2);%-mean(y3(1,BEGINo-10:BEGINo-1));;
                [h,p,ci,stats]=ttest(aa1,bb1);
                pvals_tmp(t)=p;
            end
            figure; plot(x,pvals_tmp);

            for statstat=1;
                r1=mean(mean(tmpeff(:,BEGINo:ENDo,redundantnumber),3),2);%-mean(mean(tmpeff(:,BEGINo-5:BEGINo-1,redundantnumber),3),2);
                c1=mean(mean(tmpeff(:,BEGINo:ENDo,1),3),2);%-mean(mean(tmpeff(:,BEGINo-5:BEGINo-1,1),3),2);
                d1=mean(mean(tmpeff(:,BEGINo:ENDo,10),3),2);%-mean(mean(tmpeff(:,BEGINo-5:BEGINo-1,10),3),2);

                %[h,p1,ci,stats1]=ttest(r1,c1);
                [h,p2,ci,stats2]=ttest2(d1,r1);
            end
            yy=ylim;figure;
            errorbar_groups([mean(c1) mean(r1) mean(d1)]',[std(c1)/sqrt(size(c1,1)-1) std(r1)/sqrt(size(r1,1)-1) std(d1)/sqrt(size(d1,1)-1)]', 'bar_colors',[1 1 1; .5 .5 .5; 1 0 0]);
            ylabel('z-scores','FontSize',16,'FontWeight','bold'); xlabel('time (sec)','FontSize',16,'FontWeight','bold'); set(gca,'FontSize',16,'FontWeight','bold');
            set(gcf,'Color','w');
            title(strcat('T_D_D','=',num2str(round(stats2.tstat*100)/100),', p=',num2str(round(p2*1000)/1000)));
        else
            figure; plot(x,y1,'b','LineWidth',3); hold on;plot(x,y2,'k','LineWidth',3); hold on;plot(x,y3,'r','LineWidth',3);yy=ylim;
            xlim([-.2 1]); title ('Average of cell responses', 'FontSize', 16, 'FontWeight', 'bold');
        end

        for manyplotthing=1;
            figure;
            for cnd=1:10;
                y1=mean(tmpeff(:,:,cnd),1);
                y1=y1-min(y1(1,1:baseEND));
                if cnd==1; err1=std(a1,0,1)./sqrt(numcells);
                    err=std(a1,0,1)./sqrt(numcells);

                elseif cnd==10;
                    err=std(a3,0,1)./sqrt(numcells);
                else
                    err=std(a2,0,1)./sqrt(numcells);
                end
                subplot(1,10,cnd); shadedErrorBar(x,y1,err,'lineProps','k'); set(gcf,'Color','w'); xlim([-.2 .8]);
                if cnd==1;
                    ylabel('z-scores','FontSize',16,'FontWeight','bold'); xlabel('time (sec)','FontSize',16,'FontWeight','bold'); set(gca,'FontSize',16,'FontWeight','bold');
                    title(titles{cnd}, 'FontSize', 18, 'FontWeight', 'bold');
                else
                    title(titles{cnd}, 'FontSize', 18, 'FontWeight', 'bold'); yticks([]);set(gca,'FontSize',16,'FontWeight','bold');
                end
            end
            for fixaxes=1

                ymin1=1; ymax1=-1;
                for cnd=1:10
                    subplot(1,10,cnd); yy=ylim;
                    if yy(1)<ymin1; ymin1=yy(1); end
                    if yy(2)>ymax1; ymax1=yy(2); end

                end
                for cnd=1:10
                    subplot(1,10,cnd); ylim([ymin1 ymax1]);
                end

            end
        end

        for manyplotthing=2;
            figure;
            ccnds={1 redundantnumber 10};
            ttls=[1 5 10];
            colars={'k' 'b' 'r'};
            for cnd=1:3;
                y1=mean(mean(tmpeff(:,:,ccnds{cnd}),3),1);
                y1=y1;%-min(y1(1,1:baseEND));
                if cnd==1; err1=std(a1,0,1)./sqrt(numcells/2);
                    err=std(a1,0,1)./sqrt(numcells/2);

                elseif cnd==3;
                    err=std(a3,0,1)./sqrt(numcells/2);
                else
                    err=std(a2,0,1)./sqrt(numcells/2);
                end
                if cnd==1;
                    subplot(1,3,cnd); shadedErrorBar(x,y1,err,'lineprops','k','patchSaturation',0); set(gcf,'Color','w'); xlim([-.2 .8]);
                end
                if cnd==2;
                    subplot(1,3,cnd); shadedErrorBar(x,y1,err,'lineProps','k'); xlim([-.2 .8]);
                end
                if cnd==3;
                    subplot(1,3,cnd); shadedErrorBar(x,y1,err,'lineProps','r'); xlim([-.2 .8]);
                end
                if cnd==1;
                    ylabel('z-scores','FontSize',16,'FontWeight','bold'); xlabel('time (sec)','FontSize',16,'FontWeight','bold'); set(gca,'FontSize',16,'FontWeight','bold');
                    title(titles{ttls(cnd)}, 'FontSize', 18, 'FontWeight', 'bold');
                else
                    if cnd==2;
                        title(titles{5}, 'FontSize', 18, 'FontWeight', 'bold'); yticks([]);set(gca,'FontSize',16,'FontWeight','bold');
                    else
                        title(titles{ttls(cnd)}, 'FontSize', 18, 'FontWeight', 'bold'); yticks([]);set(gca,'FontSize',16,'FontWeight','bold');

                    end
                end
            end
            for fixaxes=1

                ymin1=1; ymax1=-1;
                for cnd=1:3
                    subplot(1,3,cnd); yy=ylim;
                    if yy(1)<ymin1; ymin1=yy(1); end
                    if yy(2)>ymax1; ymax1=yy(2); end

                end
                for cnd=1:3
                    subplot(1,3,cnd); ylim([ymin1 ymax1]);
                end

            end
        end


    end



    for rasters=1
        sortby=1:10;
        BOOSTscalefactor=50; %between 0 and 100. higher is brighter
        figure;
        effects1=tmpeff;
        ak1=mean(effects1(:,:,sortby),3); [ix,b]=sort(mean(ak1(:,BEGINo:ENDo),2),'descend');
        ak1=mean(effects1(:,:,sortby(:)),3); [mmm,idx]=max(ak1(:,BEGINo:ENDo),[],2); [ix,b]=sort(idx);

        effects2=effects1(b,:,:);

        for z1=1:size(effects2,3)
            subplot(1,size(effects2,3),z1); imagesc(timeaxis,1:size(effects2,1),effects2(:,:,z1)); colormap('Parula'); set(gca,'FontSize',14,'FontWeight','bold');
            xlim([-.2 .8]);
            title(titles{z1});

            if z1==1; xlabel('sec','FontSize',14,'FontWeight','bold');
                ylabel('neurons','FontSize',14,'FontWeight','bold');
            end
        end
        for settingcolorscale=1;
            scalefactor=(BOOSTscalefactor/100); if scalefactor==1; scalefactor=.99; end
            tmpmax=0; tmpmin=0;
            for z1=1:size(effects2,3)
                subplot(1,size(effects2,3),z1);
                cb=caxis; if cb(1)<tmpmin; tmpmin=cb(1); end
                if cb(2)>tmpmax; tmpmax=cb(2); end
            end
            for z1=1:size(effects2,3)
                subplot(1,size(effects2,3),z1);
                cb=[tmpmin tmpmax];
                caxis(cb*(1-scalefactor))
            end
        end
    end

    for limited_rasters=1;

        sortby=[4 10];
        BOOSTscalefactor=77; %between 0 and 100. higher is brighter
        figure;
        titles2={'standard pairing';'oddball pairing'};
        for rasters=1

            effects1=tmpeff;
           ak1=max(effects1(:,:,sortby),[],3); [ix,b]=sort(mean(ak1(:,BEGINo:ENDo),2),'descend');
           %  ak1=mean(effects1(:,:,sortby(:)),3); [mmm,idx]=max(ak1(:,BEGINo:ENDo),[],2); [ix,b]=sort(idx);
            effects2=effects1(b,:,:);

            for z1=1:size(sortby,2)
                subplot(1,size(sortby,2),z1); imagesc(timeaxis,1:size(effects2,1),effects2(:,:,sortby(z1))); colormap('Parula'); set(gca,'FontSize',14,'FontWeight','bold');
                xlim([-.1 .8]);
                title(titles2{z1});

                if z1==1; xlabel('sec','FontSize',14,'FontWeight','bold');
                    ylabel('neurons','FontSize',14,'FontWeight','bold');
                end
            end
            for settingcolorscale=1;
                scalefactor=(BOOSTscalefactor/100); if scalefactor==1; scalefactor=.99; end
                tmpmax=0; tmpmin=0;
                for z1=1:size(sortby,2)
                    subplot(1,size(sortby,2),z1);
                    cb=caxis; if cb(1)<tmpmin; tmpmin=cb(1); end
                    if cb(2)>tmpmax; tmpmax=cb(2); end
                end
                for z1=1:size(sortby,2)
                    subplot(1,size(sortby,2),z1);
                    cb=[tmpmin tmpmax];
                    caxis(cb*(1-scalefactor))
                    caxis([-1 12]);
                end
            end
        end
    end
    make_eps_saveable


    for sortseparately=1;
        if sortseparately==1;

            for limited_rasters=1;

                sortby=[4 10];
                BOOSTscalefactor=77; %between 0 and 100. higher is brighter
                figure;
                titles2={'standard pairing';'oddball pairing'};
                for rasters=1

                   

                    for z1=1:size(sortby,2)
                         effects1=tmpeff;
                        ak1=mean(effects1(:,:,sortby(z1)),3); [ix,b]=sort(max(ak1(:,BEGINo:ENDo+5),[],2),'descend');
                    %  ak1=mean(effects1(:,:,sortby(:)),3); [mmm,idx]=max(ak1(:,BEGINo:ENDo),[],2); [ix,b]=sort(idx);
                        effects2=effects1(b,:,:);
                        subplot(1,size(sortby,2),z1); imagesc(timeaxis,1:size(effects2,1),effects2(:,:,sortby(z1))); colormap('Parula'); set(gca,'FontSize',14,'FontWeight','bold');
                        xlim([-.1 .9]);
                        title(titles2{z1});

                        if z1==1; xlabel('sec','FontSize',14,'FontWeight','bold');
                            ylabel('neurons','FontSize',14,'FontWeight','bold');
                        end
                    end
                    for settingcolorscale=1;
                        scalefactor=(BOOSTscalefactor/100); if scalefactor==1; scalefactor=.99; end
                        tmpmax=0; tmpmin=0;
                        for z1=1:size(sortby,2)
                            subplot(1,size(sortby,2),z1);
                            cb=caxis; if cb(1)<tmpmin; tmpmin=cb(1); end
                            if cb(2)>tmpmax; tmpmax=cb(2); end
                        end
                        for z1=1:size(sortby,2)
                            subplot(1,size(sortby,2),z1);
                            cb=[tmpmin tmpmax];
                            caxis(cb*(1-scalefactor))
                            caxis([0 14]);
                        end
                    end
                end
            end
            make_eps_saveable

        end
    end

    
end


% [p,table] = anova_rm({datV1multi1 datPTLpmulti1});