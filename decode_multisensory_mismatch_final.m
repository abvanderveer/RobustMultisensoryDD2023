
clcl
%%% decoding analysis from the paper
for file_load=1;
    %%% ptlp files
    files= {};


end

for file=1:size(files,1);
    clearvars -except file files mouseavs nummacells mousevec mouseavsTIME mouseavs_ran
    limtrials1=10; %limit the number of trials analysed for each stimulus type(control/redundandt/deviant).
    baseBEGIN=10; baseEND=14; %what you call the baseline
    BEGINo=15; ENDo=28; %these are the analysis windows. note: stimulus starts at 15 and finishes at 28
    onestim=0; %if 0, then plot responses to both stimuli (e.g. 45 and 135 deg). if 1, then plot only to first orientation. if 2, then plot to second

    %%%%%%%this compiles all the files into three matrices
    for buildingmatrices=1;

        count=0; counttest=0;
        effectsall=[];  %this is the main varable
        effectsallTR=[];  %this is the main varable
        actcont=[];
        countall=0;
    end

    load(strcat('CA_workspace_SORTED_',files{file}));

    load(strcat('CA_workspace_SORTEDdifferent_',files{file}),'vectrials')


    framerate=28;

    disp(files{file});
    for remove_first_trial_of_each_type=1;
        if remove_first_trial_of_each_type==1;
            if limtrials1>0;
                limtrials=limtrials1;
                for ss1=1:2;
                    vecstim=vectrials{ss1,2};
                    tmps1=[];
                    for tr1=1:size(vecstim,2);

                        if and(sum(vecstim(tr1)==tmps1)==0,vecstim(tr1)~=15)
                            tmps1=horzcat(tmps1,vecstim(tr1));
                            %vecstim(tr1)=99;
                        end

                    end
                    for additionally_remove_first_few_controls=1;
                        tmps1=[]; gh=size(vecstim,2)+1;
                        disp(sum(vecstim==(20+(ss1-1))));
                        for tr1=1:size(vecstim,2);
                            if and(vecstim(gh-tr1)>15,vecstim(gh-tr1)<22)
                                if sum(tmps1==1)<(limtrials)
                                    tmps1=horzcat(tmps1,1);
                                else
                                    vecstim(gh-tr1)=99;
                                end
                            end
                        end
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
        timeaxis=-(samp*(framerate/2)):samp:(samp*framerate); timeaxis=timeaxis+.088;

        for stim=1:2; %gotta combiny stimmos
            limtrials=limtrials1;
            %vecdat1=vectrials{1,1};vecdat2=vectrials{1,1}; vecdat=vecdat1; vecdat(:,:,(size(vecdat1,3)+1):(size(vecdat1,3)+(size(vecdat2,3))))=vecdat2;
            vecdat=vectrials{stim,1};
            disp(size(vecdat,1));
            mousetmp=zeros(size(vecdat,1),1)+file;
            %vecstim=horzcat(vectrials{1,2},vectrials{1,2}); for t=1:size(vecstim,2); if vecstim(1,t)==20; vecstim(1,t)=21; end; end;
            vecstim=vectrials{stim,2};
            for renamevecstim=0;
                if renamevecstim==1;
                    vecstim1=vecstim;
                    for t=1:size(vecstim,2);
                        if or(or(vecstim(1,t)==1,vecstim(1,t)==1),vecstim(1,t)==1)
                            vecstim1(1,t)=1;
                        elseif or(or(vecstim(1,t)==2,vecstim(1,t)==2),vecstim(1,t)==2)
                            vecstim1(1,t)=2;
                        elseif or(or(vecstim(1,t)==3,vecstim(1,t)==3),vecstim(1,t)==4)
                            vecstim1(1,t)=3;
                        elseif or(or(vecstim(1,t)==5,vecstim(1,t)==5),vecstim(1,t)==6)
                            vecstim1(1,t)=4;
                        end
                    end
                    vecstim=vecstim1;
                end
            end

            if limtrials1>0;
                for checkingfortoo_few_trials=1;
                    ss=[1 2 3 15 phzstd(stim)];
                    for gss=1:5;
                        sss(gss)=sum(vecstim==ss(gss));
                    end
                    if min(sss)<limtrials1;
                        limtrials=min(sss);
                    end
                    disp(limtrials);
                end

                effects=[]; effectsTR=[];
                for z1=1:8
                    tmp1=vecdat(:,:,vecstim==z1);
                    if size(tmp1,3)>(limtrials-1)
                        effects(:,:,z1+1)=mean(tmp1(:,:,1:limtrials),3);
                    else
                        effects(:,:,z1+1)=mean(tmp1(:,:,1:end),3);
                    end
                end
                for z1=1:4
                    tmp1=vecdat(:,:,vecstim==z1);
                    try
                        effectsTR(:,:,:,z1+1)=(tmp1(:,:,1:limtrials));
                    catch
                        aab=size(tmp1,3);
                        effectsTR(:,:,1:aab,z1+1)=(tmp1(:,:,1:aab));
                        aac=limtrials-aab;
                        effectsTR(:,:,(aab+1):limtrials,z1+1)=(tmp1(:,:,1:aac));
                    end
                end

                tmp1=vecdat(:,:,vecstim==15);
                if size(tmp1,3)>(limtrials-1)
                    effects(:,:,10)=mean(tmp1(:,:,1:limtrials),3);
                else
                    effects(:,:,10)=mean(tmp1(:,:,1:end),3);
                end
                tmp1=vecdat(:,:,vecstim==15);
                effectsTR(:,:,:,6)=(tmp1(:,:,1:limtrials));

                tmp1=vecdat(:,:,vecstim==phzstd(stim));

                if size(tmp1,3)>(limtrials-1)
                    effects(:,:,1)=mean(tmp1(:,:,1:limtrials),3);
                else
                    effects(:,:,1)=mean(tmp1(:,:,1:end),3);
                end
                tmp1=vecdat(:,:,vecstim==phzstd(stim));
                effectsTR(:,:,:,1)=(tmp1(:,:,1:limtrials));



            end


            for c1=1:size(effects,1) %baseline correction
                basos=[1 2 2 2 2 2 2 2 2 2; 1 2 2 2 2 2 2 2 2 2];
                at=squeeze(effects(c1,baseBEGIN:baseEND,:));
                for t=1:size(phases,2); if phases(1,t)>2; phases(1,t)=2; end; end;
                for z1=1:size(effects,3)
                    mnn=mean(effects(c1,baseBEGIN:baseEND,z1),2); %mnn=0;
                    stt=std(at(:));

                    if stt<.001; stt=std(effects(c1,ENDo:end,z1)); end
                    if stt<.001; stt=10000; end
                    for t1=1:size(effects,2)
                        effects(c1,t1,z1)=(effects(c1,t1,z1)-mnn)./stt;
                    end
                    if z1<7
                        for tri=1:size(effectsTR,3);
                            % effectsTR(c1,:,tri,z1)=effectsTR(c1,:,tri,z1)-min(effectsTR(c1,baseBEGIN:baseEND,tri,z1))./stt;
                            effectsTR(c1,:,tri,z1)=(effectsTR(c1,:,tri,z1)-mnn)./stt;
                        end
                    end
                    mnns(file,stim,c1,z1)=mnn;
                    stts(file,stim,c1,z1)=stt;
                end

            end



            effectsall((1+countall):(countall+size(effects,1)),1:size(effects,2),2:10,stim)=effects(:,1:end,2:10);
            effectsall((1+countall):(countall+size(effects,1)),1:size(effects,2),1,stim)=effects(:,1:end,1);
            effectsallTR((1+countall):(countall+size(effects,1)),1:size(effects,2),:,2:6,stim)=effectsTR(:,1:end,:,2:6);
            effectsallTR((1+countall):(countall+size(effects,1)),1:size(effects,2),:,1,stim)=effectsTR(:,1:end,:,1);


        end

        countall=countall+size(effects,1);
    end



    close all force

    %stimflip=[2;1];
    clear devvec typvec ranvec
    for stim=1:2;
        clear effectsallTR1
        %%%analyse
        effectsallTR1(:,:,:,1)=mean(effectsallTR(:,:,:,4,stim),4);effectsallTR1(:,:,:,2)=mean(effectsallTR(:,:,:,6,stim),4);

        %
        time1=BEGINo; time2=ENDo;%time2=BEGINo+3;%
        trs=size(effectsallTR1,3);
        sampa=5; %how many you use to create the average

        for ttt=1; %1:size(effectsallTR1,2);%BEGINo:ENDo;%1:size(effectsallTR1,2);
            typacc=[]; devacc=[]; ranacc=[];
            for z=1:1000;
                avset=randsample(trs,sampa);
                tsset=[]; for tr=1:trs; if sum(avset==tr)<1; tsset=vertcat(tsset,tr); end; end;

                %mnvecTYP=mean(mean(effectsallTR1(:,ttt,avset,1),3),2);
                %mnvecDEV=mean(mean(effectsallTR1(:,ttt,avset,2),3),2);

                % mnvecTYP=mean(mean(effectsallTR1(:,time1:time2,avset,1),3),2);
                %mnvecDEV=mean(mean(effectsallTR1(:,time1:time2,avset,2),3),2);

                %fullvec
                mnvecTYP=mean(effectsallTR1(:,time1:time2,avset,1),3); mnvecTYP=mnvecTYP(:);
                mnvecDEV=mean(effectsallTR1(:,time1:time2,avset,2),3); mnvecDEV=mnvecDEV(:);

                %typ identify;
                idents=[]; idents1=[];
                for t=1:(trs-sampa);
                    %a=mean(effectsallTR1(:,time1:time2,tsset(t),1),2);% a=a./max(a);
                    %a=mean(effectsallTR1(:,ttt,tsset(t),1),2);
                    a=(effectsallTR1(:,time1:time2,tsset(t),1)); a=a(:);
                    b=mnvecTYP;%./max(vertcat(mnvecTYP));
                    c=mnvecDEV;%./max(vertcat(mnvecDEV));
                    d=a(randsample(size(a,1),size(a,1)));
                    % simil1=dot(a,b)./((dot(a,a)+dot(b,b))/2);
                    % simil2=dot(a,c)./((dot(a,a)+dot(c,c))/2);

                    % simil3=dot(d,b)./((dot(d,d)+dot(b,b))/2);
                    % simil4=dot(d,c)./((dot(d,d)+dot(c,c))/2);
                    r=corrcoef(a,b); simil1=r(2,1);
                    r=corrcoef(a,c); simil2=r(2,1);
                    r=corrcoef(d,b); simil3=r(2,1);
                    r=corrcoef(d,c); simil4=r(2,1);

                    if simil1>simil2;
                        idents(t)=1;
                    elseif simil1<simil2;
                        idents(t)=0;
                    else
                        idents(t)=0;
                    end


                    if simil3>simil4;
                        idents1(t)=0;
                    elseif simil3<simil4;
                        idents1(t)=1;
                    else
                        idents1(t)=0;
                    end
                end
                typacc=vertcat(typacc,idents');
                ranacc=vertcat(ranacc,idents1');

                %dev identify;
                idents=[]; idents1=[];
                for t=1:(trs-sampa);
                    %a=mean(effectsallTR1(:,time1:time2,tsset(t),2),2);% a=a./max(a);
                    %a=mean(effectsallTR1(:,ttt,tsset(t),2),2);% a=a./max(a);
                    a=(effectsallTR1(:,time1:time2,tsset(t),2)); a=a(:);
                    b=mnvecTYP;%./max(vertcat(mnvecTYP));
                    c=mnvecDEV;%./max(vertcat(mnvecDEV));
                    d=a(randsample(size(a,1),size(a,1)));
                    % simil1=dot(a,b)./((dot(a,a)+dot(b,b))/2);
                    % simil2=dot(a,c)./((dot(a,a)+dot(c,c))/2);

                    %  simil3=dot(d,b)./((dot(d,d)+dot(b,b))/2);
                    %  simil4=dot(d,c)./((dot(d,d)+dot(c,c))/2);
                    r=corrcoef(a,b); simil1=r(2,1);
                    r=corrcoef(a,c); simil2=r(2,1);
                    r=corrcoef(d,b); simil3=r(2,1);
                    r=corrcoef(d,c); simil4=r(2,1);

                    if simil1<simil2;
                        idents(t)=1;
                    elseif simil1>simil2;
                        idents(t)=0;
                    else
                        idents(t)=.0;
                    end

                    if simil3>simil4;
                        idents1(t)=0;
                    elseif simil3<simil4;
                        idents1(t)=1;
                    else
                        idents1(t)=0;
                    end
                end
                devacc=vertcat(devacc,idents');
                ranacc=vertcat(ranacc,idents1');
            end
            devvec(ttt,stim)=mean(devacc);
            typvec(ttt,stim)=mean(typacc);
            ranvec(ttt,stim)=mean(ranacc);
        end

    end

    aa=((mean(devvec,2)+mean(typvec,2))/2);%
    mouseavs(file)=mean(aa);
    mouseavsTIME(file,:)=aa';

    nummacells(file)=size(c,1)./((time2-time1)+1);

    aa=((mean(ranvec,2)+mean(ranvec,2))/2);
    mouseavs_ran(file)=mean(aa);
end

excl=nummacells-nummacells;
for m=1:size(nummacells,2);
    if nummacells(m)<20;
        excl(m)=1;
    end
end
nummacells=nummacells(excl<1);mouseavs=mouseavs(excl<1); mouseavsTIME=mouseavsTIME(excl<1,:);mouseavs_ran=mouseavs_ran(excl<1); %mousevec=mousevec(excl<1);
mouseavs1=[];

mouseavs
mean(mouseavs)

figure;errorbar_groups([mean(mouseavs) mean(mouseavs_ran)],[std(mouseavs)./sqrt(file) std(mouseavs_ran)./sqrt(file)], 'bar_colors',[.5 .5 .5]);hold on; scatter((mouseavs-mouseavs)+1,mouseavs,200,'k'); ylim([0 1]); hold on; scatter((mouseavs-mouseavs)+2,mouseavs_ran,200,'k'); ylim([0 1]);make_eps_saveable

