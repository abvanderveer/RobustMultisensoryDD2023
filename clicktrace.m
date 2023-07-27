clcl
%read the titles of each plot and the command line for instructions
%ask jordan if you have any questions.
res=256; framerate=10;
for loadfiles=1;
    
    file2='2023_0626_Y0052_Blue_PV_multi002'; %put filename here
    excludecellsfromhalo=1; nopca=0;
    
    load cmapGRAY
    scrsz=get(0,'ScreenSize');
    invv=0;
    
    for loadraw=1;
        fileb=strcat(file2,'.tif');
        file=fileb;
        ii=1;
        res1=res;
        res2=res;
        
        datraw=zeros(res1,res2,1200);
        done=0;
        
        datraw=bigread(file);
        av1=mean(datraw,3);
        i=size(datraw,3);
        datraw=single(datraw(:,:,1:i-1));
        if invv==1;
            datraw=(datraw-max(datraw(:))).*-1;
        end
    end
end

for focus_time_window=1;
    if focus_time_window==1;
        figure;  plot(squeeze(mean(mean(datraw,2),1))); title('select beginning and end of desired time window');
        [xap,yap]=ginput(2);
        begin1=round(xap(1)); end1=round(xap(2)); if begin1<1; begin1=1; end; if end1>size(datraw,3); end1=size(datraw,3); end;
        dat=datraw(:,:,begin1:end1);
    end
end

%clear scoresraw_8*
vars=zeros(res,res); load cmapGRAYred
vars1=zeros(res,res);
vars2=zeros(res,res); h=waitbar(0,'Computing template images');  z=randsample(size(dat,3),round(size(dat,3)/10));
%z=1:2:size(dat,3); 
for r=1:size(datraw,1);
    for c=1:size(datraw,2);
%        akak=(smooth(squeeze(dat(r,c,:)),10,'lowess')); 
        akak=((squeeze(dat(r,c,z)))); 
        vars2(r,c)=std(akak); 
        vars1(r,c)=max(akak)-min(akak);
        if isnan(vars1(r,c));
            vars1(r,c)=0;
        else if (vars1(r,c))<0;
                vars1(r,c)=0;
            end
        end
        if isnan(vars2(r,c));
            vars1(r,c)=0;
        else if (vars2(r,c))<0;
                vars1(r,c)=0;
            end
        end
    end
    waitbar(r/res,h);
end
close force

switchy=1; vars=vars1;
figure; imagesc(vars); set(gca,'Colormap',cmapGRAY); title('adjust colormap to accentuate cells');
cells=200;
graytmp=get(gca,'Colormap');
close all
blk2=vars-vars;
blkrej=vars-vars;
cellpixels=[];
cell=1; smoothy=0; 
 

finishedd=0; gbgb=0;
figure; sz=get(0,'Screensize');posful=sz;  set(gcf,'Position',posful);
subplot(3,2,[1 3]); imagesc(zeros(res,res)); subplot(3,2,5:6); plot(1:100,zeros(100));
while finishedd==0;
    yt=0; testo=1; 
    btn = uicontrol('Style', 'pushbutton', 'String', 'FINISHED','Position', [10 100 100 200],'Callback','finishedd=1; close; cells=cell-1;');
    btn3 = uicontrol('Style', 'pushbutton', 'String', 'RESCORE LAST CELL','Position', [10 350 100 100],'Callback','cell=cell-1');
    try
    while yt<testo;
        btn2 = uicontrol('Style', 'pushbutton', 'String', 'SWITCH IMAGE','Position', [10 500 100 200],'Callback','switchy=1+switchy; gbgb=1');
        if rem(switchy,2)==0; vars=vars2; else vars=vars1; end;
        subplot(3,2,[2 4]); imagesc(sign(blk2-blkrej));  colormap(cmapGRAYred); freezeColors; caxis([-1 1]); title(strcat('pick cell ',num2str(cell))); set(gca,'XTick',0:10:250,'YTick',0:10:250,'TickDir','out');% grid minor;
        %pos2=get(gcf,'Position');
        
        subplot(3,2,[1 3]);  imagesc(vars);  colormap(graytmp); set(gca,'XTick',0:10:250,'YTick',0:10:250,'TickDir','out'); title('click button once, then on gray part of figure next to button, then click below trace');
        %pos1=get(gcf,'Position');
        [y1,x1]=ginput(2); x(1)=round(min(x1)); x(2)=round(max(x1)); y(1)=round(min(y1)); y(2)=round(max(y1)); graytmp=get(gca,'Colormap');
        if x(1)<1; x(1)=1; end; if x(2)<1; x(2)=1; end; if y(1)<1; y(1)=1; end; if y(2)<1; y(2)=1; end; if x(1)>res1; x(1)=res1; end; if x(2)>res1; x(2)=res1; end; if y(1)>res1; y(1)=res1; end; if y(2)>res1; y(2)=res1; end;
        if (x(2)-x(1))>20; x(2)=x(1)+20; end; if (y(2)-y(1))>20; y(2)=y(1)+20; end;
        tmpwin=datraw(x(1):x(2),y(1):y(2),begin1:end1);
        sqr=[x(1) x(2) y(1) y(2)];
        posful=get(gcf,'Position');  set(gcf,'Position',posful);
        
        blk=vars-vars; blk(x(1):x(2),y(1):y(2))=std(tmpwin,0,3);
        subplot(3,2,[1 3]); imagesc(blk);  colormap(cmapGRAY); 
        clear tmpdat
        for r=1:size(tmpwin,1);
            for c=1:size(tmpwin,2);
                for z=1:size(tmpwin,3);
                    tmpdat(z,c+((r-1)*size(tmpwin,2)))=tmpwin(r,c,z);
                end
            end
        end
        
        
        if nopca==0;
        [coeff,scorea,latent,tsquared,explained,mu] = pca(tmpdat); %variables in columns
        else
            coeff=ones(size(tmpdat,2),size(tmpdat,2));
            abc=mean(tmpdat',2);
            prpr=prctile(abc,20);
            for t5=1:size(coeff,1);
                if abc(t5,1)<prpr;
                    coeff(t5,1)=0;
                end
            end
        end
        aa=max(coeff(:,1));
        for p=1:size(coeff,1);
            if coeff(p,1)<aa*.8;
                coeff(p,1)=0;
            end
        end
        for fixa=1;
            if nopca==0;
            if or(sum(coeff(:,1)>0)<6,sum(coeff(:,1)>0)>.75*size(coeff,1));
                [coeff,scorea,latent,tsquared,explained,mu] = pca(tmpdat); %variables in columns
                aaa=sort(coeff(:,1));
                aa=aaa(round(size(coeff,1)*.85),1);
                for p=1:size(coeff,1);
                    if coeff(p,1)<aa;
                        coeff(p,1)=0;
                    end
                end
            end
            end
        end
        comptmp=zeros(size(tmpwin,1),size(tmpwin,2));
       
        for r=1:size(tmpwin,1);
            comptmp(r,:)=coeff(((r-1)*size(tmpwin,2))+1:r*size(tmpwin,2),1)';
        end
        blk=vars-vars; blk(x(1):x(2),y(1):y(2))=comptmp;
        blktmp=blk; blktmp1=blktmp;
        kk=0; a=zeros(size(datraw,3),1);
        for xa=1:res1;
            for ya=1:res2;
                if abs(blk(xa,ya))>0;
                    kk=1+kk;
                    for t=1:size(datraw,3);
                        a(t,1)=datraw(xa,ya,t)+a(t,1);
                    end
                end
            end
        end
        a=a/kk;
        for makehaloFIRSTPASS=1;
            pixelstmp2=[];
            sr=sqr;
            if sqr(1)<3; sr(1)=3; end; if sqr(2)>res1-3; sr(2)=res1-3; end; if sqr(3)<3; sr(3)=3; end; if sqr(4)>res2-3; sr(4)=res2-3; end;
            for c=sr(3)-2:sr(4)+2;
                pixelstmp2=vertcat(pixelstmp2,[sr(1)-1 c]);
            end
            for c=sr(3)-2:sr(4)+2;
                pixelstmp2=vertcat(pixelstmp2,[sr(1)-2 c]);
            end
            for c=sr(3)-2:sr(4)+2;
                pixelstmp2=vertcat(pixelstmp2,[sr(2)+1 c]);
            end
            for c=sr(3)-2:sr(4)+2;
                pixelstmp2=vertcat(pixelstmp2,[sr(2)+2 c]);
            end
            for r=sr(1):sr(2);
                pixelstmp2=vertcat(pixelstmp2,[r sr(3)-1]);
            end
            for r=sr(1):sr(2);
                pixelstmp2=vertcat(pixelstmp2,[r sr(3)-2]);
            end
            for r=sr(1):sr(2);
                pixelstmp2=vertcat(pixelstmp2,[r sr(4)+1]);
            end
            for r=sr(1):sr(2);
                pixelstmp2=vertcat(pixelstmp2,[r sr(4)+2]);
            end
            haltmp=zeros(size(datraw,3),1);
            for cp=1:size(pixelstmp2,1);
                for t=1:size(datraw,3);
                    haltmp(t,1)=datraw(pixelstmp2(cp,1),pixelstmp2(cp,2),t)+haltmp(t,1);
                end
            end
            for cp=1:size(pixelstmp2,1);
                for xxx=1:res1;
                    for yyy=1:res2;
                        if and(pixelstmp2(cp,1)==xxx,pixelstmp2(cp,2)==yyy);
                            blktmp(xxx,yyy)=max(blk(:));
                        end
                    end
                end
            end
            
        end
        haltmp=haltmp/cp;
        subplot(3,2,[1 3]); imagesc(blktmp);   colormap(graytmp);  set(gca,'XTick',0:10:250,'YTick',0:10:250,'TickDir','out');
        atmpa=a-haltmp; %atmpa(1:end-1)=diff(smooth(atmpa,framerate*2,'lowess')); 
        ckk=max(atmpa(framerate*60+1:end-(framerate*60))); 
        ck1=sort(atmpa(framerate*60+1:end-(framerate*60))); ckk=ck1(end-5); 
        tm1=0; for t=1:size(atmpa,1); if (atmpa(t))==ckk; tm1=t; end; end; tm2=round(size(atmpa,1)/2);
        if smoothy==1;  %smooth or not
            subplot(3,6,[17 18]); plot(smooth(a-haltmp,framerate,'lowess')); tm2=tm2/framerate;  %xlim([tm2-200 tm2+200]);      
            subplot(3,6,13:16); plot(smooth(a-haltmp,framerate,'lowess')); tm1=tm1/framerate; %    xlim([tm2-400 tm2+400]);   
        else
            subplot(3,6,[17 18]); plot(a-haltmp); tm2=tm2/framerate; % xlim([tm2-200 tm2+200]);      
            subplot(3,6,13:16); plot(a-haltmp); tm1=tm1/framerate;  % xlim([tm2-400 tm2+400]);   
        end
        %xlim([tm2-200 tm2+200]); 
        title('click above trace if good, below if bad');  [xt,yt]=ginput(1);  xt=round(xt); 
        testo=a(xt,1)-haltmp(xt,1); 
        %if yt<0; blkrej=blkrej+blktmp1; end
        if yt<testo; blkrej=blkrej+blktmp1; end
    end
    disp(strcat('scored : ',num2str(cell)));
    
    for t=1:size(datraw,3);
        scoresraw(cell,t)=a(t);
    end
    compweights(:,:,cell)=blk;
    blk2(x(1):x(2),y(1):y(2))=comptmp;
    cx=round(mean(x)); cy=round(mean(y)); if cx<4; cx=4; end; if cy<4; cy=4; end; if cx>res1-3; cx=res1-3; end; if cy>res2-3; cy=res2-3; end;
    pixelstmp=[];
    for r=sqr(1):sqr(2);
        for c=sqr(3):sqr(4);
            pixelstmp=vertcat(pixelstmp,[r c]);
        end
    end
    
    
    loc{1,cell}=[cx-1 cy-1; cx-1 cy; cx-1 cy+1; cx cy-1; cx cy; cx cy+1; cx+1 cy-1; cx+1 cy; cx+1 cy+1];
    cellpixels=vertcat(cellpixels,pixelstmp);
    for makehalo=1;
        pixelstmp2=[];
        sr=sqr;
        if sqr(1)<3; sr(1)=3; end; if sqr(2)>res1-3; sr(2)=res1-3; end; if sqr(3)<3; sr(3)=3; end; if sqr(4)>res2-3; sr(4)=res2-3; end;
        for c=sr(3)-2:sr(4)+2;
            pixelstmp2=vertcat(pixelstmp2,[sr(1)-1 c]);
        end
        for c=sr(3)-2:sr(4)+2;
            pixelstmp2=vertcat(pixelstmp2,[sr(1)-2 c]);
        end
        for c=sr(3)-2:sr(4)+2;
            pixelstmp2=vertcat(pixelstmp2,[sr(2)+1 c]);
        end
        for c=sr(3)-2:sr(4)+2;
            pixelstmp2=vertcat(pixelstmp2,[sr(2)+2 c]);
        end
        for r=sr(1):sr(2);
            pixelstmp2=vertcat(pixelstmp2,[r sr(3)-1]);
        end
        for r=sr(1):sr(2);
            pixelstmp2=vertcat(pixelstmp2,[r sr(3)-2]);
        end
        for r=sr(1):sr(2);
            pixelstmp2=vertcat(pixelstmp2,[r sr(4)+1]);
        end
        for r=sr(1):sr(2);
            pixelstmp2=vertcat(pixelstmp2,[r sr(4)+2]);
        end
        
    end
    halos{1,cell}=pixelstmp2;
    zz=zeros(res,res); for r=1:size(pixelstmp,1); zz(pixelstmp(r,1),pixelstmp(r,2))=1; end; for r=1:size(pixelstmp2,1); zz(pixelstmp2(r,1),pixelstmp2(r,2))=.5; end; 
    cell=1+cell;    
   %save('tmpwrk', '-regexp', '^(?!(dat|datraw)$).'); %this will save your progress in case the program crashes. to start where you left off, you'll still need to run everything before this loop. then load this workspace.
    
    end
end%actually does the ROI selection

for visualizing_results=1;
    %     cellpixels=[];
    %     for cc=1:cells;
    %         blk=compweights(:,:,cc);
    %         for x1=1:res1;
    %             for y1=1:res2;
    %                 if blk(x1,y1)>0;
    %                     cellpixels=vertcat(cellpixels,[x1 y1]);
    %                 end
    %             end
    %         end
    %     end
    halosraw=scoresraw-scoresraw;
    for c=1:cells; %for making halos;
        tmpmat=halos{1,c};
        for aa=1:size(tmpmat,1);
            if tmpmat(aa,1)>res;
                tmpmat(aa,1)=res;
            end
            if tmpmat(aa,2)>res;
                tmpmat(aa,2)=res;
            end
            if tmpmat(aa,1)<1;
                tmpmat(aa,1)=1;
            end
            if tmpmat(aa,2)<1;
                tmpmat(aa,2)=1;
            end
        end
        tmppix=[];
        halotracetmp=scoresraw(c,:)-scoresraw(c,:);
        for c1=1:size(tmpmat,1);
            k=0;
            if excludecellsfromhalo==1; 
            for c2=1:size(cellpixels);
                if sum(tmpmat(c1,:)==cellpixels(c2,:))==2;
                    k=1;
                end
            end
            end
            if k==0;
                tmppix=vertcat(tmppix,tmpmat(c1,:));
                for t=1:size(halotracetmp,2);
                    halotracetmp(1,t)=halotracetmp(1,t)+datraw(tmpmat(c1,1),tmpmat(c1,2),t);
                end
            end
        end
        if size(tmppix,1)>0;
            halosraw(c,:)=halotracetmp./size(tmppix,1);
        end
        disp(size(tmppix,1));
    end
    for c=1:cells;
        clear scoretmp
        for t=1:size(halosraw,2);
            if (scoresraw(c,t)-halosraw(c,t))<0;
                scoretmp(1,t)=0;
            else scoretmp(1,t)=(scoresraw(c,t)-halosraw(c,t));
            end
        end
        m=max(scoretmp);
        scores(c,:)=(scoretmp)./m;
    end
    for c=1:cells;
        m=max(scoresraw(c,:));
        scores1(c,:)=(scoresraw(c,:))./m;
    end
    
    
    figure; imagesc(scores1); set(gca,'Colormap',cmapGRAY); caxis([0 max(scores1(:))]);
    figure; imagesc(scores); set(gca,'Colormap',cmapGRAY); caxis([0 max(scores(:))]); title('halosubtracted');
    
    for c=1:cells;
        tmp=loc{1,c}(:,:);
        centroids(c,:)=mean(tmp,1);
    end
    tmp=zeros(res1,res2);
    for cell=1:cells;
        tmp(round(centroids(cell,1)),round(centroids(cell,2)))=cell;
    end
    figure; imagesc(tmp); set(gca,'Colormap',cmapGRAY); caxis([0 1]);
    figure; imagesc(blk2); set(gca,'Colormap',cmapGRAY);
end %compiles all youve done.
dlmwrite(strcat('TracesRAW_',file2,'.csv'),scoresraw,',');
dlmwrite(strcat('halosRAW_',file2,'.csv'),halosraw,',');
save(strcat('CONTOURS_',file2),'loc');
save(strcat('compweights_',file2),'compweights');
save(strcat('halopixels_',file2),'halos');
save(strcat('cellboxpixels_',file2),'cellpixels');

for redtrace=0;
    if redtrace==0;
        file3=horzcat(file2,'_REDCHAN');
        datraw2=squeeze(mean(mean(datraw,1),2));
        dlmwrite(strcat('TracesRAW_',file3,'.csv'),datraw2',',');
    elseif redtrace==1;
        file3=horzcat(file2,'_REDCHAN');
        fileb=strcat(file3,'.tif');
        file=fileb;
        ii=1;
        
        datraw=zeros(res1,res2,1200);
        done=0;
        i=ii;
        datraw=bigread(file);
        av1=mean(datraw,3);x(1)=round(min(x1)); x(2)=round(max(x1)); y(1)=round(min(y1)); y(2)=round(max(y1));
        
        tmpwin=datraw(x(1):x(2),y(1):y(2),1:end);
        datraw=squeeze(mean(mean(tmpwin,1),2));
      
        
        load cmapGRAY
        figure; imagesc(av1); set(gca,'Colormap',cmapGRAY);
        
        figure; plot(datraw);
        dlmwrite(strcat('TracesRAW_',file3,'.csv'),datraw',',');
    end
end
for checkmakesurecellrescore=0;
    if checkmakesurecellrescore==1;
        
        cellwrong=44;
        figure; imagesc(mean(compweights(:,:,1:cellwrong-1),3));
        figure; imagesc(compweights(:,:,cellwrong));
    end
end

