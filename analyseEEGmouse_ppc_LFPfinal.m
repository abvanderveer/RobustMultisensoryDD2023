clear all 
close all
clc

load PPC_LFP_multisensory_workspace
sex=[0 0 1 1 1 1 1 1 0 0 0 0]; age= [159 153 91 92 92 84 142 142 142 142 108 108]; %%%0 is female, 1 is male!!!!!!!!
inc=[1 1 1 1 1 1 1 0 1 1 0 1]; subs2=sum(inc); itcs_dev=itcs_dev(:,:,:,inc==1);itcs_all=itcs_all(:,:,inc==1); age1=age(inc==1); sex1=sex(inc==1);
ersps_dev=ersps_dev(:,:,:,inc==1);ersps_all=ersps_all(:,:,inc==1); 
%%%fourth in ersps_dev is the deviant. 

for rebaseCorr=1;
    for f=1:subs2;
        for fr=1:size(ersps_dev,1);
            for cond=1:size(ersps_dev,3);
                 bb=mean(mean(ersps_dev(fr,(basetime2-4):time1-1,cond,f)));
                 ersps_dev(fr,:,cond,f)=ersps_dev(fr,:,cond,f)-bb;
            end
            bb=mean(mean(ersps_all(fr,basetime2:time1-1,f)));
            ersps_all(fr,:,f)=ersps_all(fr,:,f)-bb;
        end
    end
end



figure; imagesc(t_sG1,f_sG1,mean(itcs_all,3)); axis xy; colormap jet; %caxis([.05 .15]);
figure;
for tr=1:4;
    subplot(1,4,tr); imagesc(t_sG1,f_sG1,mean(itcs_dev(:,:,tr,:),4)); axis xy; colormap jet; caxis([0 .1]);
end

figure; imagesc(t_sG1,f_sG1,mean(ersps_all,3)); axis xy; colormap jet;%caxis([-2 1]);
figure;
for tr=1:4;
    subplot(1,4,tr); imagesc(mean(ersps_dev(:,:,tr,:),4)); axis xy; colormap jet; caxis([-1 1.5]);

end

for powercompare=1;


for fr=1:size(ersps_dev,1); %ttests
    for t=1:size(ersps_dev,2);
        [h,p,ci,stats]=ttest(squeeze(mean(ersps_dev(fr,t,3,:),3)),squeeze(ersps_dev(fr,t,4,:)));
        ppvals(fr,t)=1-p;
    end
end
figure; contourf(t_sG1,f_sG1,ppvals,110, 'linecolor', 'none'); axis xy; colormap gray; xlim([-50 550]); ylim([2 65]); set(gca,'yscale','log'); caxis([.95 1]);make_eps_saveable
figure; subplot(1,2,1); contourf(t_sG1,f_sG1,mean(mean(ersps_dev(:,:,3,:),3),4),40, 'linecolor', 'none'); axis xy; colormap jet; xlim([-100 835]); ylim([2 65]); set(gca,'yscale','log');
subplot(1,2,2); contourf(t_sG1,f_sG1,mean(mean(ersps_dev(:,:,4,:),3),4),40, 'linecolor', 'none'); axis xy; colormap jet; xlim([-100 835]); ylim([2 65]); set(gca,'yscale','log'); cc=caxis;
subplot(1,2,1); caxis(cc);title('typical');
subplot(1,2,2);  caxis(cc);title('mismatch');make_eps_saveable


% figure; imagesc(ppvals); axis xy; colormap gray;caxis([.9 1]);
% [xx,yy]=ginput(2); xx=round(xx); yy=round(yy);
% fs=min(yy):max(yy); ts=min(xx):max(xx);
% a=squeeze(mean(mean(ersps_dev(fs,ts,3,:),1),2));b=squeeze(mean(mean(ersps_dev(fs,ts,4,:),1),2));

%figure; scatter(sex1',b);xlim([-.5 1.5]);
figure; subplot(1,2,1); contourf(t_sG1,f_sG1,mean(mean(ersps_dev(:,:,4,sex1==0),3),4),40, 'linecolor', 'none'); axis xy; colormap jet; xlim([-100 835]); ylim([2 65]); set(gca,'yscale','log');
subplot(1,2,2); contourf(t_sG1,f_sG1,mean(mean(ersps_dev(:,:,4,sex1==1),3),4),40, 'linecolor', 'none'); axis xy; colormap jet; xlim([-100 835]); ylim([2 65]); set(gca,'yscale','log');cc=caxis;
subplot(1,2,1);  caxis(cc); title('female');
subplot(1,2,2);  caxis(cc);title('male');make_eps_saveable


%%%%% comparing sex
clear pows
for f=1:subs2;%make avs
    for tr=1:4;
        pows(:,f,tr)=mean(ersps_dev(:,time1:time2,tr,f),2);
    end
end
for fr=1:size(pows,1); %ttests
    [h,p,ci,stats]=ttest(mean(pows(fr,:,3),3),pows(fr,:,4));
    pvals(fr)=p;
end
figure; plot(pvals);
figure;
freqs=40:55; 
freqs=2:9; 
plot(f_sG1,mean(mean(pows(:,:,3),3),2)); hold on
plot(f_sG1,mean(pows(:,:,4),2)); hold on
figure; scatter(age1,mean(pows(freqs,:,4),1),'fill'); 
figure; scatter(sex1,mean(pows(freqs,:,4),1),'r','fill'); xlim([-.5 1.5]);hold on; scatter(sex1-.15,mean(pows(freqs,:,3),1),'k','fill');
clear tmp1 tmp2; tmp1=mean(pows(freqs,sex1==0,3),1); tmp2=mean(pows(freqs,sex1==0,4),1); 
for x=1:size(tmp1,2); line([-.15 0],[tmp1(x) tmp2(x)]); hold on; end
clear tmp1 tmp2; tmp1=mean(pows(freqs,sex1==1,3),1); tmp2=mean(pows(freqs,sex1==1,4),1); 
for x=1:size(tmp1,2); line([.85 1],[tmp1(x) tmp2(x)]); hold on; end
 [h,p,ci,stats]=ttest(mean(pows(freqs,sex1==0,4),1),mean(pows(freqs,sex1==1,4),1))
 [h,p,ci,stats]=ttest(mean(pows(freqs,:,3),1),mean(pows(freqs,:,4),1))
 anova_rm({horzcat(mean(pows(freqs,sex1==0,3),1)',mean(pows(freqs,sex1==0,4),1)') horzcat(mean(pows(freqs,sex1==1,3),1)',mean(pows(freqs,sex1==1,4),1)')});
figure; scatter((sex1*0)+1,mean(pows(freqs,:,4),1),'r','fill'); xlim([-.5 1.5]);hold on; scatter(sex1*0,mean(pows(freqs,:,3),1),'k','fill');
for x=1:subs2; line([0 1],[mean(pows(freqs,x,3),1) mean(pows(freqs,x,4),1)]); hold on; end;

% figure; subplot(1,2,1); contourf(t_sG1,f_sG1,mean(mean(ersps_dev(:,:,3,sex1==0),3),4),40, 'linecolor', 'none'); axis xy; colormap jet; xlim([-50 650]); ylim([2 90]); set(gca,'yscale','log');
% subplot(1,2,2); contourf(t_sG1,f_sG1,mean(mean(ersps_dev(:,:,3,sex1==1),3),4),40, 'linecolor', 'none'); axis xy; colormap jet; xlim([-50 650]); ylim([2 90]); set(gca,'yscale','log');cc=caxis;
% subplot(1,2,1);  caxis(cc);
% subplot(1,2,2);  caxis(cc);
end


for itccompare=1;
clear itcs
for f=1:subs2;%make avs
    for tr=1:4;
        itcs(:,f,tr)=mean(itcs_dev(:,37:50,tr,f),2);
    end
end
for fr=1:size(itcs,1); %ttests
    [h,p,ci,stats]=ttest(mean(itcs(fr,:,3),3),itcs(fr,:,4));
    pvals(fr)=p;
end
figure; plot(pvals);

for fr=1:size(itcs,1); %ttests
    for t=1:size(itcs_dev,2);
        [h,p,ci,stats]=ttest(squeeze(mean(itcs_dev(fr,t,3,:),3)),squeeze(itcs_dev(fr,t,4,:)));
        ppvals(fr,t)=1-p;
    end
end
figure; contourf(t_sG1,f_sG1,ppvals,110, 'linecolor', 'none'); axis xy; colormap gray; xlim([-50 550]); ylim([2 65]); set(gca,'yscale','log'); caxis([.95 1]);make_eps_saveable
figure; subplot(1,2,1); contourf(t_sG1,f_sG1,mean(mean(itcs_dev(:,:,3,:),3),4),40, 'linecolor', 'none'); axis xy; colormap jet; xlim([-100 850]); ylim([2 65]); set(gca,'yscale','log');
subplot(1,2,2); contourf(t_sG1,f_sG1,mean(mean(itcs_dev(:,:,4,:),3),4),40, 'linecolor', 'none'); axis xy; colormap jet; xlim([-100 850]); ylim([2 65]); set(gca,'yscale','log'); cc=caxis;
cc(1)=.03; subplot(1,2,1);  caxis(cc);
subplot(1,2,2);  caxis(cc); 
end

%for ERP analysis
for erpanalslys=1;
for s=1:subs;
    tmp1=erpDATall{1,s};
    stp1=erpDATall{2,s};
    ERPavs1(:,1,s)=mean(tmp1(:,:),1)';
    ERPavs1(:,2,s)=mean(tmp1(stp1==5,:),1)';    
end
ERPavs=ERPavs1(:,:,inc==1);
flpo=(inc-inc)+1; flpo(8)=-1; flpo(6)=(-1); 
for s=1:size(ERPavs,3); ERPavs(:,:,s)=ERPavs(:,:,s).*flpo(s); for cnd=1:2; ERPavs(:,cnd,s)=ERPavs(:,cnd,s)-0; end; end;
figure; plot(-474:975,mean(ERPavs(:,1,:),3),'LineWidth',3,'Color','k'); xlim([-100 600]); hold on
plot(-474:975,mean(ERPavs(:,2,:),3),'LineWidth',3,'Color','r'); xlim([-100 600]);

end


