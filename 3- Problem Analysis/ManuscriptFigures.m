clear all; close all

%To find where capacity stops being constraining
if true
  SCARCITY=[0.05 0.1 0.15]
  
  

load WorkspacePT_Gauss.mat
    scarcePS_D05=max(find(round(uD1s05+uD2s05,4)==round(MaxTreat(1)*SCARCITY(1),4)));
    scarcePS_D10=max(find(round(uD1s10+uD2s10,4)==round(MaxTreat(1)*SCARCITY(2),4)));
    scarcePS_D15=max(find(round(uD1s15+uD2s15,4)==round(MaxTreat(1)*SCARCITY(3),4)));

load WorkspacePN_Gauss.mat
    scarcePN_D05=max(find(round(uD1s05+uD2s05,4)==round(MaxTreat(1)*SCARCITY(1),4)));
    scarcePN_D10=max(find(round(uD1s10+uD2s10,4)==round(MaxTreat(1)*SCARCITY(2),4)));
    scarcePN_D15=max(find(round(uD1s15+uD2s15,4)==round(MaxTreat(1)*SCARCITY(3),4)));

load WorkspaceST_Gauss.mat
    scarceSS_D05=max(find(round(uD1s05+uD2s05,4)==round(MaxTreat(1)*SCARCITY(1),4)));
    scarceSS_D10=max(find(round(uD1s10+uD2s10,4)==round(MaxTreat(1)*SCARCITY(2),4)));
    scarceSS_D15=max(find(round(uD1s15+uD2s15,4)==round(MaxTreat(1)*SCARCITY(3),4)));

    
load WorkspaceSN_Gauss.mat
    scarceSN_D05=max(find(round(uD1s05+uD2s05,4)==round(MaxTreat(1)*SCARCITY(1),4)));
    scarceSN_D10=max(find(round(uD1s10+uD2s10,4)==round(MaxTreat(1)*SCARCITY(2),4)));
    scarceSN_D15=max(find(round(uD1s15+uD2s15,4)==round(MaxTreat(1)*SCARCITY(3),4)));

load WorkspacePT_Gauss.mat
    scarcePS_V05=max(find(round(uV1s05+uV2s05,4)==round(MaxTreat(2)*SCARCITY(1),4)));
    scarcePS_V10=max(find(round(uV1s10+uV2s10,4)==round(MaxTreat(2)*SCARCITY(2),4)));
    scarcePS_V15=max(find(round(uV1s15+uV2s15,4)==round(MaxTreat(2)*SCARCITY(3),4)));

load WorkspacePN_Gauss.mat
    scarcePN_V05=max(find(round(uV1s05+uV2s05,4)==round(MaxTreat(2)*SCARCITY(1),4)));
    scarcePN_V10=max(find(round(uV1s10+uV2s10,4)==round(MaxTreat(2)*SCARCITY(2),4)));
    scarcePN_V15=max(find(round(uV1s15+uV2s15,4)==round(MaxTreat(2)*SCARCITY(3),4)));

    
load WorkspaceST_Gauss.mat
    scarceSS_V05=max(find(round(uV1s05+uV2s05,4)==round(MaxTreat(2)*SCARCITY(1),4)));
    scarceSS_V10=max(find(round(uV1s10+uV2s10,4)==round(MaxTreat(2)*SCARCITY(2),4)));
    scarceSS_V15=max(find(round(uV1s15+uV2s15,4)==round(MaxTreat(2)*SCARCITY(3),4)));
    
load WorkspaceSN_Gauss.mat
    scarceSN_V05=max(find(round(uV1s05+uV2s05,4)==round(MaxTreat(2)*SCARCITY(1),4)));
    scarceSN_V10=max(find(round(uV1s10+uV2s10,4)==round(MaxTreat(2)*SCARCITY(2),4)));
    scarceSN_V15=max(find(round(uV1s15+uV2s15,4)==round(MaxTreat(2)*SCARCITY(3),4)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%   Figures for Optimal and Ad hoc Allocations  %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%   Figures for drugs   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%   Figures to see effect of TR   %%%%%%%%%%%%%%%%%%%%%%%%%  

%Figure A1: Drug deviation to see effect of TR with permanent immunity
if true

bl=[0, 0.4470, 0.7410];
re=[0.8500, 0.3250, 0.0980];

    
    
fig= figure

ts=ts05;

load WorkspacePT_Gauss.mat
    subplot1=subplot(2,2,1)
    p1=plot(ts,uD1s10,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,uD2s10,'LineWidth',3,'Color',re);hold on
    p1=plot(ts,uD1s10_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on
    p2=plot(ts,uD2s10_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    xlim([0 ts(scarcePS_D10)])
    ylim([0 SCARCITY(2)*MaxTreat(1)])
    title({'Compliance to Travel Restrictions','(A)'})
    ylabel({'Prop. of','Treated'}, 'FontSize', 16); 
 
    
    subplot(2,2,3)
    p1=plot(ts,I1s10_D,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,I2s10_D,'LineWidth',3,'Color',re);hold on
    p3=plot(ts,I1s10_Dah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on 
    p4=plot(ts,I2s10_Dah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    ylabel({'Prop. of','Infected'}, 'FontSize', 16); 
    xlim([0 ts(scarcePS_D10)])
    ylim([0 0.035])
    title({'(C)'})
    legend1=legend([p1 p3 p2 p4],{'State 1: Optimal~~~','State 1: Ad Hoc~~~','State 2: Optimal~~~','State 2: Ad Hoc'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
    'Position',[0.106626425988088 0.493992323950414 0.811608069283622 0.0352380951245623],...
    'Orientation','horizontal',...
    'Interpreter','latex');

%  'Position',[0.137876385116331 0.413039942998033 0.763393865312849 0.0352380951245626],...
load WorkspacePN_Gauss.mat
    subplot(2,2,2)
    p1=plot(ts,uD1s10,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,uD2s10,'LineWidth',3,'Color',re);hold on
    p1=plot(ts,uD1s10_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on
    p2=plot(ts,uD2s10_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    xlim([0 ts(scarcePN_D10)])
    ylim([0 SCARCITY(2)*MaxTreat(1)])
    title({'Noncompliance to Travel Restrictions','(B)'})
    ax = gca;
    ax.XRuler.Exponent = 0;
    
    
    subplot(2,2,4)
load WorkspacePN_Gauss.mat
    p2=plot(ts,I1s10_D,'LineWidth',3);hold on
    p2=plot(ts,I2s10_D,'LineWidth',3);hold on
    p2=plot(ts,I1s10_Dah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on 
    p2=plot(ts,I2s10_Dah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    xlim([0 ts(scarcePN_D10)])
    ylim([0 0.1])
   title({'(D)'})
    
    
    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
 %  ylabel(han,'Quantity of Treatment', 'FontSize', 16);
    xlabel(han,'Time (months)', 'FontSize', 16);
   


   % sgtitle({'Optimal Deviation from Ad hoc Drug Allocation and','Corresponding Infection Levels for when Immunity is Permanent'},'Color','Black', 'FontSize', 18);
   % sgtitle({'Optimal and Ad Hoc Drug Allocation with','Corresponding Infection Levels for When Immunity',' Is Permanent and Treatment Capacity Is 10%'},'Color','Black', 'FontSize', 18);  
    
    saveas(gcf,'DeviationDrug_EffectOfSIP.png'); hold off
      
end

%Figure 3: Drug deviation to see effect of TR with 6-mo. immunity
if true

bl=[0, 0.4470, 0.7410];
re=[0.8500, 0.3250, 0.0980];

    
    
fig= figure

ts=ts05;

load WorkspaceST_Gauss.mat
    subplot1=subplot(2,2,1)
    p1=plot(ts,uD1s10,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,uD2s10,'LineWidth',3,'Color',re);hold on
    p1=plot(ts,uD1s10_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on
    p2=plot(ts,uD2s10_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    xlim([0 ts(scarceSS_D10)])
    ylim([0 SCARCITY(2)*MaxTreat(1)])
    title({'Compliance to Travel Restrictions','(A)'})
    ylabel({'Prop. of','Treated'}, 'FontSize', 16); 
 
    
    subplot(2,2,3)
    p1=plot(ts,I1s10_D,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,I2s10_D,'LineWidth',3,'Color',re);hold on
    p3=plot(ts,I1s10_Dah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on 
    p4=plot(ts,I2s10_Dah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    ylabel({'Prop. of','Infected'}, 'FontSize', 16); 
    xlim([0 ts(scarceSS_D10)])
    ylim([0 0.035])
    title({'(C)'})
    legend1=legend([p1 p3 p2 p4],{'State 1: Optimal~~~','State 1: Ad Hoc~~~','State 2: Optimal~~~','State 2: Ad Hoc'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
    'Position',[0.106626425988088 0.493992323950414 0.811608069283622 0.0352380951245623],...
    'Orientation','horizontal',...
    'Interpreter','latex');




load WorkspaceSN_Gauss.mat
    subplot(2,2,2)
    p1=plot(ts,uD1s10,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,uD2s10,'LineWidth',3,'Color',re);hold on
    p1=plot(ts,uD1s10_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on
    p2=plot(ts,uD2s10_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    xlim([0 ts(scarceSN_D10)])
    ylim([0 SCARCITY(2)*MaxTreat(1)])
    title({'Noncompliance to Travel Restrictions','(B)'})
    
    
    subplot(2,2,4)
    p2=plot(ts,I1s10_D,'LineWidth',3);hold on
    p2=plot(ts,I2s10_D,'LineWidth',3);hold on
    p2=plot(ts,I1s10_Dah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on 
    p2=plot(ts,I2s10_Dah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    xlim([0 ts(scarceSN_D10)])
    ylim([0 0.1])
    title({'(D)'})
    
    
    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
 %  ylabel(han,'Quantity of Treatment', 'FontSize', 16);
    xlabel(han,'Time (months)', 'FontSize', 16);
        
   % sgtitle({'Optimal Deviation from Ad hoc Drug Allocation and','Corresponding Infection Levels for when Immunity is Permanent'},'Color','Black', 'FontSize', 18);
   %   sgtitle({'Optimal and Ad Hoc Drug Allocation with','Corresponding Infection Levels for When Immunity','Lasts 6 Months and Treatment Capacity Is 10%'},'Color','Black', 'FontSize', 18);  
    saveas(gcf,'DeviationDrug_EffectOfSIP_with6mos.png'); hold off
    
      
end


%%%%%%%%%%%%%%%%   Figures to see effect of Immunity   %%%%%%%%%%%%%%%%%%%%  

%Figure 2: Drug deviation to see effect immunity with TR
if true

bl=[0, 0.4470, 0.7410];
re=[0.8500, 0.3250, 0.0980];
    
    
fig= figure

ts=ts05;

load WorkspacePT_Gauss.mat
    subplot1=subplot(2,2,1)
    p1=plot(ts,uD1s10,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,uD2s10,'LineWidth',3,'Color',re);hold on
    p1=plot(ts,uD1s10_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on
    p2=plot(ts,uD2s10_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    xlim([0 ts(scarcePS_D10)])
    ylim([0 SCARCITY(2)*MaxTreat(1)])
    title({'Permanent Immunity','(A)'})
    ylabel({'Prop. of','Treated'}, 'FontSize', 16); 
 
    
    subplot(2,2,3)
    p1=plot(ts,I1s10_D,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,I2s10_D,'LineWidth',3,'Color',re);hold on
    p3=plot(ts,I1s10_Dah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on 
    p4=plot(ts,I2s10_Dah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    ylabel({'Prop. of','Infected'}, 'FontSize', 16); 
    xlim([0 ts(scarcePS_D10)])
    ylim([0 0.035])
        title({'(C)'})
     legend1=legend([p1 p3 p2 p4],{'State 1: Optimal~~~','State 1: Ad Hoc~~~','State 2: Optimal~~~','State 2: Ad Hoc'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
    'Position',[0.106626425988088 0.493992323950414 0.811608069283622 0.0352380951245623],...
    'Orientation','horizontal',...
    'Interpreter','latex');



load WorkspaceST_Gauss.mat
    subplot(2,2,2)
    p1=plot(ts,uD1s10,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,uD2s10,'LineWidth',3,'Color',re);hold on
    p1=plot(ts,uD1s10_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on
    p2=plot(ts,uD2s10_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    xlim([0 ts(scarceSS_D10)])
    ylim([0 SCARCITY(2)*MaxTreat(1)])
    title({'6-month Immunity','(B)'})
    
    subplot(2,2,4)
    p1=plot(ts,I1s10_D,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,I2s10_D,'LineWidth',3,'Color',re);hold on
    p3=plot(ts,I1s10_Dah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on 
    p4=plot(ts,I2s10_Dah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    xlim([0 ts(scarceSS_D10)])
    ylim([0 0.035])
    title({'(D)'})
    
    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
 
    xlabel(han,'Time (months)', 'FontSize', 16);
        
   % sgtitle({'Optimal Deviation from Ad hoc Drug Allocation and','Corresponding Infection Levels for when Immunity is Permanent'},'Color','Black', 'FontSize', 18);
   % sgtitle({'Optimal and Ad Hoc Drug Allocation with','Corresponding Infection Levels for When There Is Compliance','to a Shelter-in-Place Order and Treatment Capacity Is 10%'},'Color','Black', 'FontSize', 18);  
    
  saveas(gcf,'DeviationDrug_EffectOfImmunityWithSIP.png'); hold off
    
      
end

%Figure A7: Drug deviation to see effect immunity with no TR
if true
bl=[0, 0.4470, 0.7410];
re=[0.8500, 0.3250, 0.0980];
    
    
fig= figure

ts=ts05;

load WorkspacePN_Gauss.mat
    subplot1=subplot(2,2,1)
    p1=plot(ts,uD1s10,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,uD2s10,'LineWidth',3,'Color',re);hold on
    p1=plot(ts,uD1s10_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on
    p2=plot(ts,uD2s10_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    xlim([0 ts(scarcePN_D10)])
    ylim([0 SCARCITY(2)*MaxTreat(1)])
    title({'Permanent Immunity','(A)'})
    ylabel({'Prop. of','Treated'}, 'FontSize', 16); 
 
    
    subplot(2,2,3)
    p1=plot(ts,I1s10_D,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,I2s10_D,'LineWidth',3,'Color',re);hold on
    p3=plot(ts,I1s10_Dah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on 
    p4=plot(ts,I2s10_Dah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    ylabel({'Prop. of','Infected'}, 'FontSize', 16); 
    xlim([0 ts(scarcePS_D10)])
    ylim([0 0.1])
    title({'(C)'})
      legend1=legend([p1 p3 p2 p4],{'State 1: Optimal~~~','State 1: Ad Hoc~~~','State 2: Optimal~~~','State 2: Ad Hoc'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
       'Position',[0.106626425988088 0.493992323950414 0.811608069283622 0.0352380951245623],...
    'Orientation','horizontal',...
    'Interpreter','latex');




load WorkspaceSN_Gauss.mat
    subplot(2,2,2)
    p1=plot(ts,uD1s10,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,uD2s10,'LineWidth',3,'Color',re);hold on
    p1=plot(ts,uD1s10_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on
    p2=plot(ts,uD2s10_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    xlim([0 ts(scarceSN_D10)])
    ylim([0 SCARCITY(2)*MaxTreat(1)])
    title({'6-month Immunity','(B)'})

    
    subplot(2,2,4)
    p1=plot(ts,I1s10_D,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,I2s10_D,'LineWidth',3,'Color',re);hold on
    p3=plot(ts,I1s10_Dah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on 
    p4=plot(ts,I2s10_Dah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    xlim([0 ts(scarceSN_D10)])
    ylim([0 0.1])
    title({'(D)'})
    
    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
 
    xlabel(han,'Time (months)', 'FontSize', 16);
        
   % sgtitle({'Optimal Deviation from Ad hoc Drug Allocation and','Corresponding Infection Levels for when Immunity is Permanent'},'Color','Black', 'FontSize', 18);
   %sgtitle({'Optimal and Ad Hoc Drug Allocation with','Corresponding Infection Levels for When There Is No Compliance','to a Shelter-in-Place Order and Treatment Capacity Is 10%'},'Color','Black', 'FontSize', 18);  
    
   saveas(gcf,'DeviationDrug_EffectOfImmunityWithoutSIP.png'); hold off
      
      
end

%Figure A6: Drug deviation to see effect immunity with no TR 5%
if true
bl=[0, 0.4470, 0.7410];
re=[0.8500, 0.3250, 0.0980];
    
    
fig= figure

ts=ts05;

load WorkspacePN_Gauss.mat
    subplot1=subplot(2,2,1)
    p1=plot(ts,uD1s05,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,uD2s05,'LineWidth',3,'Color',re);hold on
    p1=plot(ts,uD1s05_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on
    p2=plot(ts,uD2s05_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    xlim([0 ts(scarcePN_D05)])
    ylim([0 SCARCITY(1)*MaxTreat(1)])
    title({'Permanent Immunity','(A)'})
    ylabel({'Prop. of','Treated'}, 'FontSize', 16);  
 
    
    subplot(2,2,3)
    p1=plot(ts,I1s05_D,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,I2s05_D,'LineWidth',3,'Color',re);hold on
    p3=plot(ts,I1s05_Dah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on 
    p4=plot(ts,I2s05_Dah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    ylabel({'Prop. of','Infected'}, 'FontSize', 16); 
    xlim([0 ts(scarcePS_D05)])
    ylim([0 0.1])
    title({'(C)'})
      legend1=legend([p1 p3 p2 p4],{'State 1: Optimal~~~','State 1: Ad Hoc~~~','State 2: Optimal~~~','State 2: Ad Hoc'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
       'Position',[0.106626425988088 0.493992323950414 0.811608069283622 0.0352380951245623],...
    'Orientation','horizontal',...
    'Interpreter','latex');




load WorkspaceSN_Gauss.mat
    subplot(2,2,2)
    p1=plot(ts,uD1s05,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,uD2s05,'LineWidth',3,'Color',re);hold on
    p1=plot(ts,uD1s05_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on
    p2=plot(ts,uD2s05_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    xlim([0 ts(scarceSN_D05)])
    ylim([0 SCARCITY(1)*MaxTreat(1)])
    title({'6-month Immunity','(B)'})
    
    subplot(2,2,4)
    p1=plot(ts,I1s05_D,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,I2s05_D,'LineWidth',3,'Color',re);hold on
    p3=plot(ts,I1s05_Dah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on 
    p4=plot(ts,I2s05_Dah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    xlim([0 ts(scarceSN_D05)])
    ylim([0 0.1])
    title({'(D)'})

    
    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
 
    xlabel(han,'Time (months)', 'FontSize', 16);
        
   % sgtitle({'Optimal Deviation from Ad hoc Drug Allocation and','Corresponding Infection Levels for when Immunity is Permanent'},'Color','Black', 'FontSize', 18);
   %sgtitle({'Optimal and Ad Hoc Drug Allocation with','Corresponding Infection Levels for When There Is No Compliance','to a Shelter-in-Place Order and Treatment Capacity Is 10%'},'Color','Black', 'FontSize', 18);  
    
   saveas(gcf,'DeviationDrug_EffectOfImmunityWithoutSIP_5percent.png'); hold off
      
      
end


%%%%%%%%%%%%%%%%   Figures to see effect of Constraint   %%%%%%%%%%%%%%%%%%  

 %Figure A2: Drug deviation to see effect of Constraint with TR and permanent immunity
if true

bl=[0, 0.4470, 0.7410];
re=[0.8500, 0.3250, 0.0980];
    
    
fig= figure

ts=ts05;

load WorkspacePT_Gauss.mat
    subplot1=subplot(2,3,1)
    p1=plot(ts,uD1s05,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,uD2s05,'LineWidth',3,'Color',re);hold on
    p1=plot(ts,uD1s05_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on
    p2=plot(ts,uD2s05_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    xlim([0 ts(scarcePS_D05)])
    ylim([0 SCARCITY(1)*MaxTreat(1)])
    ylabel({'Prop. of','Treated'}, 'FontSize', 16);
    title({'5% Drug Capacity','(A)'}) 
 
    
    subplot(2,3,4)
    p1=plot(ts,I1s05_D,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,I2s05_D,'LineWidth',3,'Color',re);hold on
    p3=plot(ts,I1s05_Dah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on 
    p4=plot(ts,I2s05_Dah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    ylabel({'Prop. of','Infected'}, 'FontSize', 16); 
    xlim([0 ts(scarcePS_D05)])
    ylim([0 0.035])
      title({'(D)'}) 
     legend1=legend([p1 p3 p2 p4],{'State 1: Optimal~~~','State 1: Ad Hoc~~~','State 2: Optimal~~~','State 2: Ad Hoc'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
        'Position',[0.106626425988088 0.493992323950414 0.811608069283622 0.0352380951245623],...
    'Orientation','horizontal',...
    'Interpreter','latex');



    subplot1=subplot(2,3,2)
    p1=plot(ts,uD1s10,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,uD2s10,'LineWidth',3,'Color',re);hold on
    p1=plot(ts,uD1s10_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on
    p2=plot(ts,uD2s10_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    xlim([0 ts(scarcePS_D10)])
    ylim([0 SCARCITY(2)*MaxTreat(1)])
   % ylabel({'Prop. of','Treated'}, 'FontSize', 16);
    title({'10% Drug Capacity','(B)'}) 
 
    
    subplot(2,3,5)
    p1=plot(ts,I1s10_D,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,I2s10_D,'LineWidth',3,'Color',re);hold on
    p3=plot(ts,I1s10_Dah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on 
    p4=plot(ts,I2s10_Dah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    %ylabel({'Prop. of','Infected'}, 'FontSize', 16); 
    xlim([0 ts(scarcePS_D10)])
    ylim([0 0.035])
    title({'(E)'}) 
   % legend1=legend([p1 p3 p2 p4],{'State 1: Optimal','State 1: Ad Hoc','State 2: Optimal','State 2: Ad Hoc'},'Interpreter','latex','Orientation','vertical','Location','northeast');



    subplot1=subplot(2,3,3)
    p1=plot(ts,uD1s15,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,uD2s15,'LineWidth',3,'Color',re);hold on
    p1=plot(ts,uD1s15_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on
    p2=plot(ts,uD2s15_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    xlim([0 ts(scarcePS_D15)])
    ylim([0 SCARCITY(3)*MaxTreat(1)])
    title({'15% Drug Capacity','(C)'}) 
    
    subplot(2,3,6)
    p1=plot(ts,I1s15_D,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,I2s15_D,'LineWidth',3,'Color',re);hold on
    p3=plot(ts,I1s15_Dah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on 
    p4=plot(ts,I2s15_Dah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    xlim([0 ts(scarcePS_D15)])
    ylim([0 0.035])
    title({'(F)'}) 

    
    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    xlabel(han,'Time (months)', 'FontSize', 16);
        
     % sgtitle({'Optimal and Ad Hoc Drug Allocation with','Corresponding Infection Levels for when Immunity Is Permanent','and There Is Compliance to a Shelter-in-Place Order'},'Color','Black', 'FontSize', 18);  
    
   saveas(gcf,'DeviationDrug_EffectOfCons_withPermAndSIP.png'); hold off
   
end

 %Figure A3: Drug deviation to see effect of Constraint with no TR and permanent immunity
if true
bl=[0, 0.4470, 0.7410];
re=[0.8500, 0.3250, 0.0980];
    
    
fig= figure

ts=ts05;

load WorkspacePN_Gauss.mat

  subplot1=subplot(2,3,1)
    p1=plot(ts,uD1s05,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,uD2s05,'LineWidth',3,'Color',re);hold on
    p1=plot(ts,uD1s05_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on
    p2=plot(ts,uD2s05_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    xlim([0 ts(scarcePN_D05)])
    ylim([0 SCARCITY(1)*MaxTreat(1)])
    ylabel({'Prop. of','Treated'}, 'FontSize', 16);
    title({'5% Drug Capacity','(A)'}) 
 
    
    subplot(2,3,4)
    p1=plot(ts,I1s05_D,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,I2s05_D,'LineWidth',3,'Color',re);hold on
    p3=plot(ts,I1s05_Dah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on 
    p4=plot(ts,I2s05_Dah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    ylabel({'Prop. of','Infected'}, 'FontSize', 16); 
    xlim([0 ts(scarcePN_D05)])
    ylim([0 0.1])
    title({'(D)'}) 
     legend1=legend([p1 p3 p2 p4],{'State 1: Optimal~~~','State 1: Ad Hoc~~~','State 2: Optimal~~~','State 2: Ad Hoc'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
        'Position',[0.106626425988088 0.493992323950414 0.811608069283622 0.0352380951245623],...
    'Orientation','horizontal',...
    'Interpreter','latex');


    subplot1=subplot(2,3,2)
    p1=plot(ts,uD1s10,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,uD2s10,'LineWidth',3,'Color',re);hold on
    p1=plot(ts,uD1s10_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on
    p2=plot(ts,uD2s10_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    xlim([0 ts(scarcePN_D10)])
    ylim([0 SCARCITY(2)*MaxTreat(1)])
   % ylabel({'Prop. of','Treated'}, 'FontSize', 16);
    title({'10% Drug Capacity','(B)'}) 
 
    
    subplot(2,3,5)
    p1=plot(ts,I1s10_D,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,I2s10_D,'LineWidth',3,'Color',re);hold on
    p3=plot(ts,I1s10_Dah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on 
    p4=plot(ts,I2s10_Dah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    %ylabel({'Prop. of','Infected'}, 'FontSize', 16); 
    xlim([0 ts(scarcePN_D10)])
    ylim([0 0.1])    
    title({'(E)'}) 
    %legend1=legend([p1 p3 p2 p4],{'State 1: Optimal','State 1: Ad Hoc','State 2: Optimal','State 2: Ad Hoc'},'Interpreter','latex','Orientation','vertical','Location','northeast');



    subplot1=subplot(2,3,3)
    p1=plot(ts,uD1s15,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,uD2s15,'LineWidth',3,'Color',re);hold on
    p1=plot(ts,uD1s15_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on
    p2=plot(ts,uD2s15_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    xlim([0 ts(scarcePN_D15)])
    ylim([0 SCARCITY(3)*MaxTreat(1)])
    title({'15% Drug Capacity','(C)'}) 
 
    
    subplot(2,3,6)
    p1=plot(ts,I1s15_D,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,I2s15_D,'LineWidth',3,'Color',re);hold on
    p3=plot(ts,I1s15_Dah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on 
    p4=plot(ts,I2s15_Dah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    xlim([0 ts(scarcePN_D15)])
    ylim([0 0.1])
    title({'(F)'}) 

    
    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    xlabel(han,'Time (months)', 'FontSize', 16);
        
    % sgtitle({'Optimal and Ad Hoc Drug Allocation with','Corresponding Infection Levels for when Immunity Is Permanent','and There Is No Compliance to a Shelter-in-Place Order'},'Color','Black', 'FontSize', 18);  
      
   saveas(gcf,'DeviationDrug_EffectOfCons_withPermAndnoSIP.png'); hold off
      
end

 %Figure A4: Drug deviation to see effect of Constraint with TR and 6-mo. immunity
if true
    
bl=[0, 0.4470, 0.7410];
re=[0.8500, 0.3250, 0.0980];
    
    
fig= figure

ts=ts05;

load WorkspaceST_Gauss.mat


    subplot1=subplot(2,3,1)
    p1=plot(ts,uD1s05,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,uD2s05,'LineWidth',3,'Color',re);hold on
    p1=plot(ts,uD1s05_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on
    p2=plot(ts,uD2s05_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    xlim([0 ts(scarceSS_D05)])
    ylim([0 SCARCITY(1)*MaxTreat(1)])
    ylabel({'Prop. of','Treated'}, 'FontSize', 16);
    title({'5% Drug Capacity','(A)'}) 
 
    
    subplot(2,3,4)
    p1=plot(ts,I1s05_D,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,I2s05_D,'LineWidth',3,'Color',re);hold on
    p3=plot(ts,I1s05_Dah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on 
    p4=plot(ts,I2s05_Dah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    ylabel({'Prop. of','Infected'}, 'FontSize', 16); 
    xlim([0 ts(scarceSS_D05)])
    ylim([0 0.035])
    title({'(D)'}) 
     legend1=legend([p1 p3 p2 p4],{'State 1: Optimal~~~','State 1: Ad Hoc~~~','State 2: Optimal~~~','State 2: Ad Hoc'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
        'Position',[0.106626425988088 0.493992323950414 0.811608069283622 0.0352380951245623],...
    'Orientation','horizontal',...
    'Interpreter','latex');



    subplot1=subplot(2,3,2)
    p1=plot(ts,uD1s10,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,uD2s10,'LineWidth',3,'Color',re);hold on
    p1=plot(ts,uD1s10_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on
    p2=plot(ts,uD2s10_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    xlim([0 ts(scarceSS_D10)])
    ylim([0 SCARCITY(2)*MaxTreat(1)])
    %ylabel({'Prop. of','Treated'}, 'FontSize', 16);
    title({'10% Drug Capacity','(B)'}) 
 
    
    subplot(2,3,5)
    p1=plot(ts,I1s10_D,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,I2s10_D,'LineWidth',3,'Color',re);hold on
    p3=plot(ts,I1s10_Dah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on 
    p4=plot(ts,I2s10_Dah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    %ylabel({'Prop. of','Infected'}, 'FontSize', 16); 
    xlim([0 ts(scarceSS_D10)])
    ylim([0 0.035])
    %legend1=legend([p1 p3 p2 p4],{'State 1: Optimal','State 1: Ad Hoc','State 2: Optimal','State 2: Ad Hoc'},'Interpreter','latex','Orientation','vertical','Location','northeast');
    title({'(E)'}) 


    subplot1=subplot(2,3,3)
    p1=plot(ts,uD1s15,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,uD2s15,'LineWidth',3,'Color',re);hold on
    p1=plot(ts,uD1s15_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on
    p2=plot(ts,uD2s15_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    xlim([0 ts(scarceSS_D15)])
    ylim([0 SCARCITY(3)*MaxTreat(1)])
    title({'15% Drug Capacity','(C)'}) 
 
    
    subplot(2,3,6)
    p1=plot(ts,I1s15_D,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,I2s15_D,'LineWidth',3,'Color',re);hold on
    p3=plot(ts,I1s15_Dah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on 
    p4=plot(ts,I2s15_Dah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    xlim([0 ts(scarceSS_D15)])
    ylim([0 0.035])
    title({'(F)'}) 

    
    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    xlabel(han,'Time (months)', 'FontSize', 16);
        
   %sgtitle({'Optimal and Ad Hoc Drug Allocation with','Corresponding Infection Levels for when Immunity Lasts 6 Months','and There Is Compliance to a Shelter-in-Place Order'},'Color','Black', 'FontSize', 18);  
      
   saveas(gcf,'DeviationDrug_EffectOfCons_with6moAndSIP.png'); hold off
end

 %Figure A5: Drug deviation to see effect of Constraint with no TR and 6-mo. immunity
if true

bl=[0, 0.4470, 0.7410];
re=[0.8500, 0.3250, 0.0980];
    
    
fig= figure

ts=ts05;

load WorkspaceSN_Gauss.mat


    subplot1=subplot(2,3,1)
    p1=plot(ts,uD1s05,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,uD2s05,'LineWidth',3,'Color',re);hold on
    p1=plot(ts,uD1s05_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on
    p2=plot(ts,uD2s05_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    xlim([0 ts(scarceSN_D05)])
    ylim([0 SCARCITY(1)*MaxTreat(1)])
    ylabel({'Prop. of','Treated'}, 'FontSize', 16);
    title({'5% Drug Capacity','(A)'}) 
 
    
    subplot(2,3,4)
    p1=plot(ts,I1s05_D,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,I2s05_D,'LineWidth',3,'Color',re);hold on
    p3=plot(ts,I1s05_Dah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on 
    p4=plot(ts,I2s05_Dah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    ylabel({'Prop. of','Infected'}, 'FontSize', 16); 
    xlim([0 ts(scarceSN_D05)])
    ylim([0 0.1])    
    title({'(D)'}) 
     legend1=legend([p1 p3 p2 p4],{'State 1: Optimal~~~','State 1: Ad Hoc~~~','State 2: Optimal~~~','State 2: Ad Hoc'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
        'Position',[0.106626425988088 0.493992323950414 0.811608069283622 0.0352380951245623],...
    'Orientation','horizontal',...
    'Interpreter','latex');


    subplot1=subplot(2,3,2)
    p1=plot(ts,uD1s10,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,uD2s10,'LineWidth',3,'Color',re);hold on
    p1=plot(ts,uD1s10_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on
    p2=plot(ts,uD2s10_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    xlim([0 ts(scarceSN_D10)])
    ylim([0 SCARCITY(2)*MaxTreat(1)])
    %ylabel({'Prop. of','Treated'}, 'FontSize', 16);
    title({'10% Drug Capacity','(B)'}) 
 
    
    subplot(2,3,5)
    p1=plot(ts,I1s10_D,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,I2s10_D,'LineWidth',3,'Color',re);hold on
    p3=plot(ts,I1s10_Dah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on 
    p4=plot(ts,I2s10_Dah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
   % ylabel({'Prop. of','Infected'}, 'FontSize', 16); 
    xlim([0 ts(scarceSN_D10)])
    ylim([0 0.1])
    %legend1=legend([p1 p3 p2 p4],{'State 1: Optimal','State 1: Ad Hoc','State 2: Optimal','State 2: Ad Hoc'},'Interpreter','latex','Orientation','vertical','Location','northeast');
    title({'(E)'}) 


    subplot1=subplot(2,3,3)
    p1=plot(ts,uD1s15,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,uD2s15,'LineWidth',3,'Color',re);hold on
    p1=plot(ts,uD1s15_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on
    p2=plot(ts,uD2s15_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    xlim([0 ts(scarceSN_D15)])
    ylim([0 SCARCITY(3)*MaxTreat(1)])
    title({'15% Drug Capacity','(C)'}) 
 
    
    subplot(2,3,6)
    p1=plot(ts,I1s15_D,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,I2s15_D,'LineWidth',3,'Color',re);hold on
    p3=plot(ts,I1s15_Dah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on 
    p4=plot(ts,I2s15_Dah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    xlim([0 ts(scarceSN_D15)])
    ylim([0 0.1])
     title({'(F)'})    

    
    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    xlabel(han,'Time (months)', 'FontSize', 16);
        
   % sgtitle({'Optimal and Ad Hoc Drug Allocation with','Corresponding Infection Levels for when Immunity Lasts 6 Months','and There Is No Compliance to a Shelter-in-Place Order'},'Color','Black', 'FontSize', 18);  
      
   saveas(gcf,'DeviationDrug_EffectOfCons_with6moAndNoSIP.png'); hold off
      
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%   Figures for vaccines   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%   Figures to see effect of TR   %%%%%%%%%%%%%%%%%%%%%%%%%  

%Figure 5: Vaccine deviation to see effect of TR with Perm. immunity
if true

bl=[0, 0.4470, 0.7410];
re=[0.8500, 0.3250, 0.0980];

    
    
fig= figure

ts=ts05;

load WorkspacePT_Gauss.mat
    subplot1=subplot(2,2,1)
    p1=plot(ts,uV1s10,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,uV2s10,'LineWidth',3,'Color',re);hold on
    p1=plot(ts,uV1s10_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on
    p2=plot(ts,uV2s10_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    xlim([0 ts(scarcePS_V10)])
    ylim([0 SCARCITY(2)*MaxTreat(2)])
    title({'Compliance to Travel Restrictions','(A)'})
    ylabel({'Prop. of','Treated'}, 'FontSize', 16); 
 
    
    subplot(2,2,3)
    p1=plot(ts,I1s10,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,I2s10,'LineWidth',3,'Color',re);hold on
    p3=plot(ts,I1s10_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on 
    p4=plot(ts,I2s10_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    ylabel({'Prop. of','Infected'}, 'FontSize', 16); 
    xlim([0 ts(scarcePS_V10)])
    ylim([0 .03])
     title({'(C)'})
      legend1=legend([p1 p3 p2 p4],{'State 1: Optimal~~~','State 1: Ad Hoc~~~','State 2: Optimal~~~','State 2: Ad Hoc'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
       'Position',[0.106626425988088 0.493992323950414 0.811608069283622 0.0352380951245623],...
    'Orientation','horizontal',...
    'Interpreter','latex');




load WorkspacePN_Gauss.mat
    subplot(2,2,2)
    p1=plot(ts,uV1s10,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,uV2s10,'LineWidth',3,'Color',re);hold on
    p1=plot(ts,uV1s10_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on
    p2=plot(ts,uV2s10_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    xlim([0 ts(scarcePN_V10)])
    ylim([0 SCARCITY(2)*MaxTreat(2)])
    title({'Noncompliance to Travel Restrictions','(B)'})
    
    
    subplot(2,2,4)
    p2=plot(ts,I1s10,'LineWidth',3);hold on
    p2=plot(ts,I2s10,'LineWidth',3);hold on
    p2=plot(ts,I1s10_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on 
    p2=plot(ts,I2s10_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    xlim([0 ts(scarcePN_V10)])
    ylim([0 .09])
     title({'(D)'})
    
    
    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
 %  ylabel(han,'Quantity of Treatment', 'FontSize', 16);
    xlabel(han,'Time (months)', 'FontSize', 16);
    
   % sgtitle({'Optimal and Ad Hoc Vaccine Allocation with','Corresponding Infection Levels for When Immunity',' Is Permanent and Treatment Capacity Is 10%'},'Color','Black', 'FontSize', 18);  
    
    saveas(gcf,'DeviationVaccine_EffectOfSIP.png'); hold off
    
    
      
end

%Figure A15: Vaccine deviation to see effect of TR with 6-mo. immunity
if true

bl=[0, 0.4470, 0.7410];
re=[0.8500, 0.3250, 0.0980];

    
    
fig= figure

ts=ts05;

load WorkspaceST_Gauss.mat
    subplot1=subplot(2,2,1)
    p1=plot(ts,uV1s10,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,uV2s10,'LineWidth',3,'Color',re);hold on
    p1=plot(ts,uV1s10_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on
    p2=plot(ts,uV2s10_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    xlim([0 ts(scarceSS_V10)])
    ylim([0 SCARCITY(2)*MaxTreat(2)])
    title({'Compliance to Travel Restrictions','(A)'})
    ylabel({'Prop. of','Treated'}, 'FontSize', 16); 
 
    
    subplot(2,2,3)
    p1=plot(ts,I1s10,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,I2s10,'LineWidth',3,'Color',re);hold on
    p3=plot(ts,I1s10_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on 
    p4=plot(ts,I2s10_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    ylabel({'Prop. of','Infected'}, 'FontSize', 16); 
    xlim([0 ts(scarceSS_V10)])
    ylim([0 .03])
    title({'(C)'})
     legend1=legend([p1 p3 p2 p4],{'State 1: Optimal~~~','State 1: Ad Hoc~~~','State 2: Optimal~~~','State 2: Ad Hoc'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
       'Position',[0.106626425988088 0.493992323950414 0.811608069283622 0.0352380951245623],...
    'Orientation','horizontal',...
    'Interpreter','latex');





load WorkspaceSN_Gauss.mat
    subplot(2,2,2)
    p1=plot(ts,uV1s10,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,uV2s10,'LineWidth',3,'Color',re);hold on
    p1=plot(ts,uV1s10_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on
    p2=plot(ts,uV2s10_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    xlim([0 ts(scarceSN_V10)])
    ylim([0 SCARCITY(2)*MaxTreat(2)])
    title({'Noncompliance to Travel Restrictions','(B)'})
    
    
    subplot(2,2,4)
    p2=plot(ts,I1s10,'LineWidth',3);hold on
    p2=plot(ts,I2s10,'LineWidth',3);hold on
    p2=plot(ts,I1s10_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on 
    p2=plot(ts,I2s10_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    xlim([0 ts(scarceSN_V10)])
    ylim([0 .09])
    title({'(D)'})
    
    
    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
 %  ylabel(han,'Quantity of Treatment', 'FontSize', 16);
    xlabel(han,'Time (months)', 'FontSize', 16);
        
   %sgtitle({'Optimal and Ad Hoc Vaccine Allocation with','Corresponding Recovery Levels for when Immunity','Lasts 6 months and Treatment Capacity Is 10%'},'Color','Black', 'FontSize', 18);  
    
  saveas(gcf,'DeviationVaccine_EffectOfSIP_with6mos.png'); hold off
    
    
    
      
end

%%%%%%%%%%%%%%%%   Figures to see effect of Immunity   %%%%%%%%%%%%%%%%%%%%  

 %Figure A19: Vaccine deviation to see effect of Immunity with TR
 if true

bl=[0, 0.4470, 0.7410];
re=[0.8500, 0.3250, 0.0980];
    
    
fig= figure

ts=ts05;

load WorkspacePT_Gauss.mat
    subplot1=subplot(2,2,1)
    p1=plot(ts,uV1s10,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,uV2s10,'LineWidth',3,'Color',re);hold on
    p1=plot(ts,uV1s10_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on
    p2=plot(ts,uV2s10_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    xlim([0 ts(scarcePS_V10)])
    ylim([0 SCARCITY(2)*MaxTreat(2)])
    title({'Permanent Immunity','(A)'})
    ylabel({'Prop. of','Treated'}, 'FontSize', 16); 
 
    
    subplot(2,2,3)
    p1=plot(ts,I1s10,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,I2s10,'LineWidth',3,'Color',re);hold on
    p3=plot(ts,I1s10_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on 
    p4=plot(ts,I2s10_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    ylabel({'Prop. of','Infected'}, 'FontSize', 16); 
    xlim([0 ts(scarcePS_V10)])
    ylim([0 .03])
    title({'(C)'}) 
       legend1=legend([p1 p3 p2 p4],{'State 1: Optimal~~~','State 1: Ad Hoc~~~','State 2: Optimal~~~','State 2: Ad Hoc'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
      'Position',[0.106626425988088 0.493992323950414 0.811608069283622 0.0352380951245623],...
    'Orientation','horizontal',...
    'Interpreter','latex');





load WorkspaceST_Gauss.mat
    subplot(2,2,2)
    p1=plot(ts,uV1s10,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,uV2s10,'LineWidth',3,'Color',re);hold on
    p1=plot(ts,uV1s10_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on
    p2=plot(ts,uV2s10_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    xlim([0 ts(scarceSS_V10)])
    ylim([0 SCARCITY(2)*MaxTreat(2)])
    ylabel({'Prop. of','Treated'}, 'FontSize', 16); 
     title({'6-month Immunity','(B)'})
    subplot(2,2,4)
    p1=plot(ts,I1s10,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,I2s10,'LineWidth',3,'Color',re);hold on
    p3=plot(ts,I1s10_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on 
    p4=plot(ts,I2s10_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    ylabel({'Prop. of','Infected'}, 'FontSize', 16); 
    xlim([0 ts(scarceSS_V10)])
    ylim([0 .03])
    title({'(D)'}) 
    
    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
 
    xlabel(han,'Time (months)', 'FontSize', 16);
        
   % sgtitle({'Optimal Deviation from Ad hoc Drug Allocation and','Corresponding Infection Levels for when Immunity is Permanent'},'Color','Black', 'FontSize', 18);
 % sgtitle({'Optimal and Ad Hoc Vaccine Allocation with','Corresponding Recovery Levels for When There Is Compliance','to a Selter-in-Place Order and Treatment Capacity Is 10%'},'Color','Black', 'FontSize', 18);  
    
   saveas(gcf,'DeviationVaccine_EffectOfImmunityWithSIP.png'); hold off
    
    
    
      
 end

 %Figure A20: Vaccine deviation to see effect of Immunity without TR
 if true

bl=[0, 0.4470, 0.7410];
re=[0.8500, 0.3250, 0.0980];
    
    
fig= figure

ts=ts05;

load WorkspacePN_Gauss.mat
    subplot1=subplot(2,2,1)
    p1=plot(ts,uV1s10,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,uV2s10,'LineWidth',3,'Color',re);hold on
    p1=plot(ts,uV1s10_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on
    p2=plot(ts,uV2s10_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    xlim([0 ts(scarcePN_V10)])
    ylim([0 SCARCITY(2)*MaxTreat(2)])
    title({'Permanent Immunity','(A)'})
    ylabel({'Prop. of','Treated'}, 'FontSize', 16); 
 
    
    subplot(2,2,3)
    p1=plot(ts,I1s10,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,I2s10,'LineWidth',3,'Color',re);hold on
    p3=plot(ts,I1s10_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on 
    p4=plot(ts,I2s10_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    ylabel({'Prop. of','Infected'}, 'FontSize', 16); 
    xlim([0 ts(scarcePN_V10)])
    ylim([0 .09])
        title({'(C)'}) 
        legend1=legend([p1 p3 p2 p4],{'State 1: Optimal~~~','State 1: Ad Hoc~~~','State 2: Optimal~~~','State 2: Ad Hoc'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
     'Position',[0.106626425988088 0.493992323950414 0.811608069283622 0.0352380951245623],...
    'Orientation','horizontal',...
    'Interpreter','latex');




load WorkspaceSN_Gauss.mat
    subplot(2,2,2)
    p1=plot(ts,uV1s10,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,uV2s10,'LineWidth',3,'Color',re);hold on
    p1=plot(ts,uV1s10_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on
    p2=plot(ts,uV2s10_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    xlim([0 ts(scarceSN_V10)])
    ylim([0 SCARCITY(2)*MaxTreat(2)])
    ylabel({'Prop. of','Treated'}, 'FontSize', 16); 
    title({'6-month Immunity','(B)'})
    
    subplot(2,2,4)
    p1=plot(ts,I1s10,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,I2s10,'LineWidth',3,'Color',re);hold on
    p3=plot(ts,I1s10_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on 
    p4=plot(ts,I2s10_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    ylabel({'Prop. of','Infected'}, 'FontSize', 16); 
    xlim([0 ts(scarceSN_V10)])
    ylim([0 .09])
    title({'(D)'}) 
    
    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
 
    xlabel(han,'Time (months)', 'FontSize', 16);
        
% sgtitle({'Optimal and Ad Hoc Vaccine Allocation with','Corresponding Recovery Levels for When There Is No Compliance','to a Shelter-in-Place Order and Treatment Capacity Is 10%'},'Color','Black', 'FontSize', 18);  
   
   saveas(gcf,'DeviationVaccine_EffectOfImmunityWithoutSIP.png'); hold off
    
    
    
      
 end
 
 
 %%%%%%%%%%%%%%%%   Figures to see effect of Constraint   %%%%%%%%%%%%%%%%%%  

 %Figure A16: Vaccine deviation to see effect of Constraint with TR and permanent immunity
if true

bl=[0, 0.4470, 0.7410];
re=[0.8500, 0.3250, 0.0980];

    
    
fig= figure

ts=ts05;

load WorkspacePT_Gauss.mat

    subplot1=subplot(2,3,1)
    p1=plot(ts,uV1s05,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,uV2s05,'LineWidth',3,'Color',re);hold on
    p1=plot(ts,uV1s05_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on
    p2=plot(ts,uV2s05_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    xlim([0 ts(scarcePS_V05)])
    ylim([0 SCARCITY(1)*MaxTreat(2)])
    title({'5% Vaccine Capacity','(A)'}) 
    ylabel({'Prop. of','Treated'}, 'FontSize', 16); 
 
    
    subplot(2,3,4)
    p1=plot(ts,I1s05,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,I2s05,'LineWidth',3,'Color',re);hold on
    p3=plot(ts,I1s05_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on 
    p4=plot(ts,I2s05_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    ylabel({'Prop. of','Infected'}, 'FontSize', 16); 
    xlim([0 ts(scarcePS_V05)])
    ylim([0 .031])
    title({'(D)'}) 
     legend1=legend([p1 p3 p2 p4],{'State 1: Optimal~~~','State 1: Ad Hoc~~~','State 2: Optimal~~~','State 2: Ad Hoc'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
    'Position',[0.106626425988088 0.493992323950414 0.811608069283622 0.0352380951245623],...
    'Orientation','horizontal',...
    'Interpreter','latex');


    subplot1=subplot(2,3,2)
    p1=plot(ts,uV1s10,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,uV2s10,'LineWidth',3,'Color',re);hold on
    p1=plot(ts,uV1s10_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on
    p2=plot(ts,uV2s10_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    xlim([0 ts(scarcePS_V10)])
    ylim([0 SCARCITY(2)*MaxTreat(2)])
    title({'10% Vaccine Capacity','(B)'}) 
    %ylabel({'Prop. of','Treated'}, 'FontSize', 16); 

    subplot(2,3,5)
    p1=plot(ts,I1s10,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,I2s10,'LineWidth',3,'Color',re);hold on
    p3=plot(ts,I1s10_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on 
    p4=plot(ts,I2s10_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    %ylabel({'Prop. of','Recovered'}, 'FontSize', 16); 
    xlim([0 ts(scarcePS_V10)])
    ylim([0 .031])
    title({'(E)'}) 
    %legend1=legend([p1 p3 p2 p4],{'State 1: Optimal','State 1: Ad Hoc','State 2: Optimal','State 2: Ad Hoc'},'Interpreter','latex','Orientation','vertical','Location','southeast');


    subplot(2,3,3)
    p1=plot(ts,uV1s15,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,uV2s15,'LineWidth',3,'Color',re);hold on
    p1=plot(ts,uV1s15_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on
    p2=plot(ts,uV2s15_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    xlim([0 ts(scarcePS_V15)])
    ylim([0 SCARCITY(3)*MaxTreat(2)])
    title({'15% Vaccine Capacity','(C)'}) 
    
    
    subplot(2,3,6)
    p2=plot(ts,I1s15,'LineWidth',3);hold on
    p2=plot(ts,I2s15,'LineWidth',3);hold on
    p2=plot(ts,I1s15_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on 
    p2=plot(ts,I2s15_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    xlim([0 ts(scarcePS_V15)])
    ylim([0 .031])
    title({'(F)'})

    
    
    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    xlabel(han,'Time (months)', 'FontSize', 16);
        
 %   sgtitle({'Optimal and Ad Hoc Vaccine Allocation with','Corresponding Infection Levels for when Immunity Is Permanent','and There Is Compliance to a Shelter-in-Place Order'},'Color','Black', 'FontSize', 18);  
    
   saveas(gcf,'DeviationVaccine_EffectOfCons_withPermAndSIP.png'); hold off
    
      
end 

%Figure 6: Vaccine deviation to see effect of Constraint with no TR and permanent immunity
if true

bl=[0, 0.4470, 0.7410];
re=[0.8500, 0.3250, 0.0980];

    
    
fig= figure

ts=ts05;

load WorkspacePN_Gauss.mat

    subplot1=subplot(2,3,1)
    p1=plot(ts,uV1s05,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,uV2s05,'LineWidth',3,'Color',re);hold on
    p1=plot(ts,uV1s05_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on
    p2=plot(ts,uV2s05_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    xlim([0 ts(scarcePN_V05)])
    ylim([0 SCARCITY(1)*MaxTreat(2)])
    title({'5% Vaccine Capacity','(A)'}) 
    ylabel({'Prop. of','Treated'}, 'FontSize', 16); 
 
    
    subplot(2,3,4)
    p1=plot(ts,I1s05,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,I2s05,'LineWidth',3,'Color',re);hold on
    p3=plot(ts,I1s05_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on 
    p4=plot(ts,I2s05_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    ylabel({'Prop. of','Infected'}, 'FontSize', 16); 
    xlim([0 ts(scarcePN_V05)])
    ylim([0 .09])
     title({'(D)'})
     legend1=legend([p1 p3 p2 p4],{'State 1: Optimal~~~','State 1: Ad Hoc~~~','State 2: Optimal~~~','State 2: Ad Hoc'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
        'Position',[0.106626425988088 0.493992323950414 0.811608069283622 0.0352380951245623],...
    'Orientation','horizontal',...
    'Interpreter','latex');




    subplot1=subplot(2,3,2)
    p1=plot(ts,uV1s10,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,uV2s10,'LineWidth',3,'Color',re);hold on
    p1=plot(ts,uV1s10_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on
    p2=plot(ts,uV2s10_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    xlim([0 ts(scarcePN_V10)])
    ylim([0 SCARCITY(2)*MaxTreat(2)])
    title({'10% Vaccine Capacity','(B)'}) 
    %ylabel({'Prop. of','Treated'}, 'FontSize', 16); 
 
    
    subplot(2,3,5)
    p1=plot(ts,I1s10,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,I2s10,'LineWidth',3,'Color',re);hold on
    p3=plot(ts,I1s10_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on 
    p4=plot(ts,I2s10_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    %ylabel({'Prop. of','Recovered'}, 'FontSize', 16); 
    xlim([0 ts(scarcePN_V10)])
    ylim([0 .09])
    title({'(E)'})
    %legend1=legend([p1 p3 p2 p4],{'State 1: Optimal','State 1: Ad Hoc','State 2: Optimal','State 2: Ad Hoc'},'Interpreter','latex','Orientation','vertical','Location','northwest');



    subplot(2,3,3)
    p1=plot(ts,uV1s15,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,uV2s15,'LineWidth',3,'Color',re);hold on
    p1=plot(ts,uV1s15_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on
    p2=plot(ts,uV2s15_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    xlim([0 ts(scarcePN_V15)])
    ylim([0 SCARCITY(3)*MaxTreat(2)])
    title({'15% Vaccine Capacity','(C)'}) 
    
    
    subplot(2,3,6)
    p2=plot(ts,I1s15,'LineWidth',3);hold on
    p2=plot(ts,I2s15,'LineWidth',3);hold on
    p2=plot(ts,I1s15_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on 
    p2=plot(ts,I2s15_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    xlim([0 ts(scarcePN_V15)])
    ylim([0 .09])
    title({'(F)'})
    
    
    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    xlabel(han,'Time (months)', 'FontSize', 16);
        
  % sgtitle({'Optimal and Ad Hoc Vaccine Allocation with','Corresponding Infection Levels for when Immunity Is Permanent','and There Is No Compliance to a Shelter-in-Place Order'},'Color','Black', 'FontSize', 18);  
   
      saveas(gcf,'DeviationVaccine_EffectOfCons_withPermAndnoSIP.png'); hold off
    
    
    
end 

%Figure A17: Vaccine deviation to see effect of Constraint with TR and 6-mo. immunity
if true

bl=[0, 0.4470, 0.7410];
re=[0.8500, 0.3250, 0.0980];

    
    
fig= figure

ts=ts05;

load WorkspaceST_Gauss.mat

    subplot1=subplot(2,3,1)
    p1=plot(ts,uV1s05,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,uV2s05,'LineWidth',3,'Color',re);hold on
    p1=plot(ts,uV1s05_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on
    p2=plot(ts,uV2s05_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    xlim([0 ts(scarceSS_V05)])
    ylim([0 SCARCITY(1)*MaxTreat(2)])
    title({'5% Vaccine Capacity','(A)'}) 
    ylabel({'Prop. of','Treated'}, 'FontSize', 16); 
 
    
    subplot(2,3,4)
    p1=plot(ts,I1s05,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,I2s05,'LineWidth',3,'Color',re);hold on
    p3=plot(ts,I1s05_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on 
    p4=plot(ts,I2s05_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    ylabel({'Prop. of','Infectied'}, 'FontSize', 16); 
    xlim([0 ts(scarceSS_V05)])
    ylim([0 .033])
         title({'(D)'})
     legend1=legend([p1 p3 p2 p4],{'State 1: Optimal~~~','State 1: Ad Hoc~~~','State 2: Optimal~~~','State 2: Ad Hoc'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
    'Position',[0.106626425988088 0.493992323950414 0.811608069283622 0.0352380951245623],...
    'Orientation','horizontal',...
    'Interpreter','latex');



    subplot1=subplot(2,3,2)
    p1=plot(ts,uV1s10,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,uV2s10,'LineWidth',3,'Color',re);hold on
    p1=plot(ts,uV1s10_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on
    p2=plot(ts,uV2s10_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    xlim([0 ts(scarceSS_V10)])
    ylim([0 SCARCITY(2)*MaxTreat(2)])
    title({'10% Vaccine Capacity','(B)'}) 
    %ylabel({'Prop. of','Treated'}, 'FontSize', 16); 
 
    
    subplot(2,3,5)
    p1=plot(ts,I1s10,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,I2s10,'LineWidth',3,'Color',re);hold on
    p3=plot(ts,I1s10_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on 
    p4=plot(ts,I2s10_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
   % ylabel({'Prop. of','Recovered'}, 'FontSize', 16); 
    xlim([0 ts(scarceSS_V10)])
    ylim([0 .033])
    %legend1=legend([p1 p3 p2 p4],{'State 1: Optimal','State 1: Ad Hoc','State 2: Optimal','State 2: Ad Hoc'},'Interpreter','latex','Orientation','vertical','Location','southeast');
     title({'(E)'})


load WorkspaceST_Gauss.mat
    subplot(2,3,3)
    p1=plot(ts,uV1s15,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,uV2s15,'LineWidth',3,'Color',re);hold on
    p1=plot(ts,uV1s15_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on
    p2=plot(ts,uV2s15_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    xlim([0 ts(scarceSS_V15)])
    ylim([0 SCARCITY(3)*MaxTreat(2)])
    title({'15% Vaccine Capacity','(C)'}) 
    
    
    subplot(2,3,6)
    p2=plot(ts,I1s15,'LineWidth',3);hold on
    p2=plot(ts,I2s15,'LineWidth',3);hold on
    p2=plot(ts,I1s15_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on 
    p2=plot(ts,I2s15_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    xlim([0 ts(scarceSS_V15)])
    ylim([0 .033])
     title({'(F)'})
    
    
    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    xlabel(han,'Time (months)', 'FontSize', 16);
        
%  sgtitle({'Optimal and Ad Hoc Vaccine Allocation with','Corresponding Infection Levels for when Immunity Lasts 6 Months','and There Is Compliance to a Shelter-in-Place Order'},'Color','Black', 'FontSize', 18);  
      
   saveas(gcf,'DeviationVaccine_EffectOfCons_with6moAndSIP.png'); hold off
    
    
    
      
end 

%Figure A18: Vaccine deviation to see effect of Constraint with no TR and 6-mo. immunity
if true

bl=[0, 0.4470, 0.7410];
re=[0.8500, 0.3250, 0.0980];

    
    
fig= figure

ts=ts05;

load WorkspaceSN_Gauss.mat

    subplot1=subplot(2,3,1)
    p1=plot(ts,uV1s05,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,uV2s05,'LineWidth',3,'Color',re);hold on
    p1=plot(ts,uV1s05_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on
    p2=plot(ts,uV2s05_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    xlim([0 ts(scarceSN_V05)])
    ylim([0 SCARCITY(1)*MaxTreat(2)])
    title({'5% Vaccine Capacity','(A)'}) 
    ylabel({'Prop. of','Treated'}, 'FontSize', 16); 
 
    
    subplot(2,3,4)
    p1=plot(ts,I1s05,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,I2s05,'LineWidth',3,'Color',re);hold on
    p3=plot(ts,I1s05_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on 
    p4=plot(ts,I2s05_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    ylabel({'Prop. of','Infected'}, 'FontSize', 16); 
    xlim([0 ts(scarceSN_V05)])
    ylim([0 .1])
    title({'(D)'})
      legend1=legend([p1 p3 p2 p4],{'State 1: Optimal~~~','State 1: Ad Hoc~~~','State 2: Optimal~~~','State 2: Ad Hoc'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
     'Position',[0.106626425988088 0.493992323950414 0.811608069283622 0.0352380951245623],...
    'Orientation','horizontal',...
    'Interpreter','latex');


    %%%
    subplot1=subplot(2,3,2)
    p1=plot(ts,uV1s10,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,uV2s10,'LineWidth',3,'Color',re);hold on
    p1=plot(ts,uV1s10_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on
    p2=plot(ts,uV2s10_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    xlim([0 ts(scarceSN_V10)])
    ylim([0 SCARCITY(2)*MaxTreat(2)])
    title({'10% Vaccine Capacity','(B)'}) 
   % ylabel({'Prop. of','Treated'}, 'FontSize', 16); 
 
    
    subplot(2,3,5)
    p1=plot(ts,I1s10,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,I2s10,'LineWidth',3,'Color',re);hold on
    p3=plot(ts,I1s10_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on 
    p4=plot(ts,I2s10_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    %ylabel({'Prop. of','Recovered'}, 'FontSize', 16); 
    xlim([0 ts(scarceSN_V10)])
    ylim([0 .1])
    %legend1=legend([p1 p3 p2 p4],{'State 1: Optimal','State 1: Ad Hoc','State 2: Optimal','State 2: Ad Hoc'},'Interpreter','latex','Orientation','vertical','Location','northwest');
     title({'(E)'})


load WorkspaceSN_Gauss.mat
    subplot(2,3,3)
    p1=plot(ts,uV1s15,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,uV2s15,'LineWidth',3,'Color',re);hold on
    p1=plot(ts,uV1s15_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on
    p2=plot(ts,uV2s15_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    xlim([0 ts(scarceSN_V15)])
    ylim([0 SCARCITY(3)*MaxTreat(2)])
    title({'15% Vaccine Capacity','(C)'}) 
    
    
    subplot(2,3,6)
    p2=plot(ts,I1s15,'LineWidth',3);hold on
    p2=plot(ts,I2s15,'LineWidth',3);hold on
    p2=plot(ts,I1s15_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on 
    p2=plot(ts,I2s15_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    xlim([0 ts(scarceSN_V15)])
    ylim([0 .1])
     title({'(F)'})    
    
    
    
    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    xlabel(han,'Time (months)', 'FontSize', 16);
        
  % sgtitle({'Optimal and Ad Hoc Vaccine Allocation with','Corresponding Infection Levels for when Immunity Lasts 6 Months','and There Is No Compliance Shelter-in-Place Order'},'Color','Black', 'FontSize', 18);  
     
   saveas(gcf,'DeviationVaccine_EffectOfCons_with6moAndNoSIP.png'); hold off
    
    
    
    
end 

 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%   Figures for cumulative number of cases  %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Drugs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    
%Figure A8: Cumulative Infections for Drugs for 5% capacity per state
if true
  
pop=1000000;


ts=ts05;

%Type of interpolation
fitApprox = fittype('pchipinterp');

load WorkspacePT_Gauss.mat
gI1s=fit(ts,(I1s05_Dah),fitApprox);
gI1s_opt=fit(ts,(I1s05_D),fitApprox);
gI2s=fit(ts,(I2s05_Dah),fitApprox);
gI2s_opt=fit(ts,(I2s05_D),fitApprox);

load workspaceNoControls.mat
gI1s_no=fit(ts,I1s_Dno1PS,fitApprox);
gI2s_no=fit(ts,I2s_Dno1PS,fitApprox);


fig=figure
subplot(2,4,1)
p1=plot(ts,integrate(gI1s_opt,ts,0)./integrate(gI1s_no,ts,0)-1,'LineWidth',3,'Color',bl); hold on
p3=plot(ts,integrate(gI1s,ts,0)./integrate(gI1s_no,ts,0)-1,'LineWidth',3,'Color',bl,'LineStyle',':'); hold on
p2=plot(ts,integrate(gI2s_opt,ts,0)./integrate(gI2s_no,ts,0)-1,'LineWidth',3,'Color',re); hold on
p4=plot(ts,integrate(gI2s,ts,0)./integrate(gI2s_no,ts,0)-1,'LineWidth',3,'Color',re,'LineStyle',':'); hold on
ylim([-0.16 0])
title({'Perm. Immunity &','Compliance to TR','(A)'})
ylabel({'Relative','Difference'}, 'FontSize', 16); 
    legend1=legend([p1 p3 p2 p4],{'State 1: Optimal~~~','State 1: Ad Hoc~~~','State 2: Optimal~~~','State 2: Ad Hoc'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
    'Position',[0.117340711702374 0.479706609664701 0.811608069283622 0.0352380951245625],...
    'Orientation','horizontal',...
    'Interpreter','latex');


subplot(2,4,5)
p1=plot(ts,(integrate(gI1s_opt,ts,0)-integrate(gI1s_no,ts,0))*pop,'LineWidth',3,'Color',bl); hold on
p3=plot(ts,(integrate(gI1s,ts,0)-integrate(gI1s_no,ts,0))*pop,'LineWidth',3,'Color',bl,'LineStyle',':'); hold on
p2=plot(ts,(integrate(gI2s_opt,ts,0)-integrate(gI2s_no,ts,0))*pop,'LineWidth',3,'Color',re); hold on
p4=plot(ts,(integrate(gI2s,ts,0)-integrate(gI2s_no,ts,0))*pop,'LineWidth',3,'Color',re,'LineStyle',':'); hold on
ylim([-13000 0])
ylabel({'Absolute Difference','per 1M People'}, 'FontSize', 16); 
title({'(E)'})



%%%
load WorkspacePN_Gauss.mat
gI1s=fit(ts,(I1s05_Dah),fitApprox);
gI1s_opt=fit(ts,(I1s05_D),fitApprox);
gI2s=fit(ts,(I2s05_Dah),fitApprox);
gI2s_opt=fit(ts,(I2s05_D),fitApprox);

load workspaceNoControls.mat
gI1s_no=fit(ts,I1s_Dno3PN,fitApprox);
gI2s_no=fit(ts,I2s_Dno3PN,fitApprox);

subplot(2,4,2)
p1=plot(ts,integrate(gI1s_opt,ts,0)./integrate(gI1s_no,ts,0)-1,'LineWidth',3,'Color',bl); hold on
p3=plot(ts,integrate(gI1s,ts,0)./integrate(gI1s_no,ts,0)-1,'LineWidth',3,'Color',bl,'LineStyle',':'); hold on
p2=plot(ts,integrate(gI2s_opt,ts,0)./integrate(gI2s_no,ts,0)-1,'LineWidth',3,'Color',re); hold on
p4=plot(ts,integrate(gI2s,ts,0)./integrate(gI2s_no,ts,0)-1,'LineWidth',3,'Color',re,'LineStyle',':'); hold on
title({'Perm. Immunity &','Noncompliance to TR','(B)'})
ylim([-0.16 0])

subplot(2,4,6)
p1=plot(ts,(integrate(gI1s_opt,ts,0)-integrate(gI1s_no,ts,0))*pop,'LineWidth',3,'Color',bl); hold on
p3=plot(ts,(integrate(gI1s,ts,0)-integrate(gI1s_no,ts,0))*pop,'LineWidth',3,'Color',bl,'LineStyle',':'); hold on
p2=plot(ts,(integrate(gI2s_opt,ts,0)-integrate(gI2s_no,ts,0))*pop,'LineWidth',3,'Color',re); hold on
p4=plot(ts,(integrate(gI2s,ts,0)-integrate(gI2s_no,ts,0))*pop,'LineWidth',3,'Color',re,'LineStyle',':'); hold on
ylim([-13000 0])
title({'(F)'})

%%%
load WorkspaceST_Gauss.mat
gI1s=fit(ts,(I1s05_Dah),fitApprox);
gI1s_opt=fit(ts,(I1s05_D),fitApprox);
gI2s=fit(ts,(I2s05_Dah),fitApprox);
gI2s_opt=fit(ts,(I2s05_D),fitApprox);

load workspaceNoControls.mat
gI1s_no=fit(ts,I1s_Dno2SS,fitApprox);
gI2s_no=fit(ts,I2s_Dno2SS,fitApprox);

subplot(2,4,3)
p1=plot(ts,integrate(gI1s_opt,ts,0)./integrate(gI1s_no,ts,0)-1,'LineWidth',3,'Color',bl); hold on
p3=plot(ts,integrate(gI1s,ts,0)./integrate(gI1s_no,ts,0)-1,'LineWidth',3,'Color',bl,'LineStyle',':'); hold on
p2=plot(ts,integrate(gI2s_opt,ts,0)./integrate(gI2s_no,ts,0)-1,'LineWidth',3,'Color',re); hold on
p4=plot(ts,integrate(gI2s,ts,0)./integrate(gI2s_no,ts,0)-1,'LineWidth',3,'Color',re,'LineStyle',':'); hold on
title({'6-mo. Immunity &','Compliance to TR','(C)'})
ylim([-0.16 0])

subplot(2,4,7)
p1=plot(ts,(integrate(gI1s_opt,ts,0)-integrate(gI1s_no,ts,0))*pop,'LineWidth',3,'Color',bl); hold on
p3=plot(ts,(integrate(gI1s,ts,0)-integrate(gI1s_no,ts,0))*pop,'LineWidth',3,'Color',bl,'LineStyle',':'); hold on
p2=plot(ts,(integrate(gI2s_opt,ts,0)-integrate(gI2s_no,ts,0))*pop,'LineWidth',3,'Color',re); hold on
p4=plot(ts,(integrate(gI2s,ts,0)-integrate(gI2s_no,ts,0))*pop,'LineWidth',3,'Color',re,'LineStyle',':'); hold on
ylim([-13000 0])
title({'(G)'})
%%%
load WorkspaceSN_Gauss.mat
gI1s=fit(ts,(I1s05_Dah),fitApprox);
gI1s_opt=fit(ts,(I1s05_D),fitApprox);
gI2s=fit(ts,(I2s05_Dah),fitApprox);
gI2s_opt=fit(ts,(I2s05_D),fitApprox);

load workspaceNoControls.mat
gI1s_no=fit(ts,I1s_Dno4SN,fitApprox);
gI2s_no=fit(ts,I2s_Dno4SN,fitApprox);

subplot(2,4,4)
p1=plot(ts,integrate(gI1s_opt,ts,0)./integrate(gI1s_no,ts,0)-1,'LineWidth',3,'Color',bl); hold on
p3=plot(ts,integrate(gI1s,ts,0)./integrate(gI1s_no,ts,0)-1,'LineWidth',3,'Color',bl,'LineStyle',':'); hold on
p2=plot(ts,integrate(gI2s_opt,ts,0)./integrate(gI2s_no,ts,0)-1,'LineWidth',3,'Color',re); hold on
p4=plot(ts,integrate(gI2s,ts,0)./integrate(gI2s_no,ts,0)-1,'LineWidth',3,'Color',re,'LineStyle',':'); hold on
title({'6-mo. Immunity &','Noncompliance to TR','(D)'})
ylim([-0.16 0])


subplot(2,4,8)
p1=plot(ts,(integrate(gI1s_opt,ts,0)-integrate(gI1s_no,ts,0))*pop,'LineWidth',3,'Color',bl); hold on
p3=plot(ts,(integrate(gI1s,ts,0)-integrate(gI1s_no,ts,0))*pop,'LineWidth',3,'Color',bl,'LineStyle',':'); hold on
p2=plot(ts,(integrate(gI2s_opt,ts,0)-integrate(gI2s_no,ts,0))*pop,'LineWidth',3,'Color',re); hold on
p4=plot(ts,(integrate(gI2s,ts,0)-integrate(gI2s_no,ts,0))*pop,'LineWidth',3,'Color',re,'LineStyle',':'); hold on
ylim([-13000 0])
title({'(H)'})

han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
%ylabel(han,'Quantity of Treatment', 'FontSize', 16);

xlabel(han,'Time (months)', 'FontSize', 16);

%sgtitle({'Cumulative Difference in the Number of Infections','Compared to the No-Treatment Scenario for Different','Drug Allocation Rules when Capacity is 5%'},'Color','Black', 'FontSize', 18);
saveas(gcf,'CumInfectionsDrug_5percent.png'); hold off

end

%Figure 4: Cumulative Infections for Drugs for 10% capacity per state
if true
  
pop=1000000;


ts=ts05;

%Type of interpolation
fitApprox = fittype('pchipinterp');

load WorkspacePT_Gauss.mat
gI1s=fit(ts,(I1s10_Dah),fitApprox);
gI1s_opt=fit(ts,(I1s10_D),fitApprox);
gI2s=fit(ts,(I2s10_Dah),fitApprox);
gI2s_opt=fit(ts,(I2s10_D),fitApprox);

load workspaceNoControls.mat
gI1s_no=fit(ts,I1s_Dno1PS,fitApprox);
gI2s_no=fit(ts,I2s_Dno1PS,fitApprox);


fig=figure
subplot(2,4,1)
p1=plot(ts,integrate(gI1s_opt,ts,0)./integrate(gI1s_no,ts,0)-1,'LineWidth',3,'Color',bl); hold on
p3=plot(ts,integrate(gI1s,ts,0)./integrate(gI1s_no,ts,0)-1,'LineWidth',3,'Color',bl,'LineStyle',':'); hold on
p2=plot(ts,integrate(gI2s_opt,ts,0)./integrate(gI2s_no,ts,0)-1,'LineWidth',3,'Color',re); hold on
p4=plot(ts,integrate(gI2s,ts,0)./integrate(gI2s_no,ts,0)-1,'LineWidth',3,'Color',re,'LineStyle',':'); hold on
ylim([-0.3 0])
title({'Perm. Immunity &','Compliance to TR','(A)'})
ylabel({'Relative','Difference'}, 'FontSize', 16); 
    legend1=legend([p1 p3 p2 p4],{'State 1: Optimal~~~','State 1: Ad Hoc~~~','State 2: Optimal~~~','State 2: Ad Hoc'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
    'Position',[0.117340711702374 0.479706609664701 0.811608069283622 0.0352380951245625],...
    'Orientation','horizontal',...
    'Interpreter','latex');


subplot(2,4,5)
p1=plot(ts,(integrate(gI1s_opt,ts,0)-integrate(gI1s_no,ts,0))*pop,'LineWidth',3,'Color',bl); hold on
p3=plot(ts,(integrate(gI1s,ts,0)-integrate(gI1s_no,ts,0))*pop,'LineWidth',3,'Color',bl,'LineStyle',':'); hold on
p2=plot(ts,(integrate(gI2s_opt,ts,0)-integrate(gI2s_no,ts,0))*pop,'LineWidth',3,'Color',re); hold on
p4=plot(ts,(integrate(gI2s,ts,0)-integrate(gI2s_no,ts,0))*pop,'LineWidth',3,'Color',re,'LineStyle',':'); hold on
ylim([-23000 0])
ylabel({'Absolute Difference','per 1M People'}, 'FontSize', 16); 
title({'(E)'})



%%%
load WorkspacePN_Gauss.mat
gI1s=fit(ts,(I1s10_Dah),fitApprox);
gI1s_opt=fit(ts,(I1s10_D),fitApprox);
gI2s=fit(ts,(I2s10_Dah),fitApprox);
gI2s_opt=fit(ts,(I2s10_D),fitApprox);

load workspaceNoControls.mat
gI1s_no=fit(ts,I1s_Dno3PN,fitApprox);
gI2s_no=fit(ts,I2s_Dno3PN,fitApprox);

subplot(2,4,2)
p1=plot(ts,integrate(gI1s_opt,ts,0)./integrate(gI1s_no,ts,0)-1,'LineWidth',3,'Color',bl); hold on
p3=plot(ts,integrate(gI1s,ts,0)./integrate(gI1s_no,ts,0)-1,'LineWidth',3,'Color',bl,'LineStyle',':'); hold on
p2=plot(ts,integrate(gI2s_opt,ts,0)./integrate(gI2s_no,ts,0)-1,'LineWidth',3,'Color',re); hold on
p4=plot(ts,integrate(gI2s,ts,0)./integrate(gI2s_no,ts,0)-1,'LineWidth',3,'Color',re,'LineStyle',':'); hold on
title({'Perm. Immunity &','Noncompliance to TR'})
ylim([-0.3 0])
title({'Perm. Immunity &','Noncompliance to TR','(B)'})

subplot(2,4,6)
p1=plot(ts,(integrate(gI1s_opt,ts,0)-integrate(gI1s_no,ts,0))*pop,'LineWidth',3,'Color',bl); hold on
p3=plot(ts,(integrate(gI1s,ts,0)-integrate(gI1s_no,ts,0))*pop,'LineWidth',3,'Color',bl,'LineStyle',':'); hold on
p2=plot(ts,(integrate(gI2s_opt,ts,0)-integrate(gI2s_no,ts,0))*pop,'LineWidth',3,'Color',re); hold on
p4=plot(ts,(integrate(gI2s,ts,0)-integrate(gI2s_no,ts,0))*pop,'LineWidth',3,'Color',re,'LineStyle',':'); hold on
ylim([-23000 0])
title({'(F)'})

%%%
load WorkspaceST_Gauss.mat
gI1s=fit(ts,(I1s10_Dah),fitApprox);
gI1s_opt=fit(ts,(I1s10_D),fitApprox);
gI2s=fit(ts,(I2s10_Dah),fitApprox);
gI2s_opt=fit(ts,(I2s10_D),fitApprox);

load workspaceNoControls.mat
gI1s_no=fit(ts,I1s_Dno2SS,fitApprox);
gI2s_no=fit(ts,I2s_Dno2SS,fitApprox);

subplot(2,4,3)
p1=plot(ts,integrate(gI1s_opt,ts,0)./integrate(gI1s_no,ts,0)-1,'LineWidth',3,'Color',bl); hold on
p3=plot(ts,integrate(gI1s,ts,0)./integrate(gI1s_no,ts,0)-1,'LineWidth',3,'Color',bl,'LineStyle',':'); hold on
p2=plot(ts,integrate(gI2s_opt,ts,0)./integrate(gI2s_no,ts,0)-1,'LineWidth',3,'Color',re); hold on
p4=plot(ts,integrate(gI2s,ts,0)./integrate(gI2s_no,ts,0)-1,'LineWidth',3,'Color',re,'LineStyle',':'); hold on
title({'6-mo. Immunity &','Compliance to TR','(C)'})
ylim([-0.3 0])

subplot(2,4,7)
p1=plot(ts,(integrate(gI1s_opt,ts,0)-integrate(gI1s_no,ts,0))*pop,'LineWidth',3,'Color',bl); hold on
p3=plot(ts,(integrate(gI1s,ts,0)-integrate(gI1s_no,ts,0))*pop,'LineWidth',3,'Color',bl,'LineStyle',':'); hold on
p2=plot(ts,(integrate(gI2s_opt,ts,0)-integrate(gI2s_no,ts,0))*pop,'LineWidth',3,'Color',re); hold on
p4=plot(ts,(integrate(gI2s,ts,0)-integrate(gI2s_no,ts,0))*pop,'LineWidth',3,'Color',re,'LineStyle',':'); hold on
ylim([-23000 0])
title({'(G)'})

%%%
load WorkspaceSN_Gauss.mat
gI1s=fit(ts,(I1s10_Dah),fitApprox);
gI1s_opt=fit(ts,(I1s10_D),fitApprox);
gI2s=fit(ts,(I2s10_Dah),fitApprox);
gI2s_opt=fit(ts,(I2s10_D),fitApprox);

load workspaceNoControls.mat
gI1s_no=fit(ts,I1s_Dno4SN,fitApprox);
gI2s_no=fit(ts,I2s_Dno4SN,fitApprox);

subplot(2,4,4)
p1=plot(ts,integrate(gI1s_opt,ts,0)./integrate(gI1s_no,ts,0)-1,'LineWidth',3,'Color',bl); hold on
p3=plot(ts,integrate(gI1s,ts,0)./integrate(gI1s_no,ts,0)-1,'LineWidth',3,'Color',bl,'LineStyle',':'); hold on
p2=plot(ts,integrate(gI2s_opt,ts,0)./integrate(gI2s_no,ts,0)-1,'LineWidth',3,'Color',re); hold on
p4=plot(ts,integrate(gI2s,ts,0)./integrate(gI2s_no,ts,0)-1,'LineWidth',3,'Color',re,'LineStyle',':'); hold on
title({'6-mo. Immunity &','Noncompliance to TR','(D)'})
ylim([-0.3 0])
%    legend1=legend([p1 p3 p2 p4],{'State 1: Optimal','State 1: Ad Hoc','State 2: Optimal','State 2: Ad Hoc'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
%set(legend1,...
%    'Position',[0.728571142469134 0.0853001236695748 0.192857428959438 0.11595238049825],...
%    'Orientation','vertical',...
%    'Interpreter','latex');
%%'Position',[0.137876385116331 0.413039942998033 0.763393865312849 0.0352380951245626],...


subplot(2,4,8)
p1=plot(ts,(integrate(gI1s_opt,ts,0)-integrate(gI1s_no,ts,0))*pop,'LineWidth',3,'Color',bl); hold on
p3=plot(ts,(integrate(gI1s,ts,0)-integrate(gI1s_no,ts,0))*pop,'LineWidth',3,'Color',bl,'LineStyle',':'); hold on
p2=plot(ts,(integrate(gI2s_opt,ts,0)-integrate(gI2s_no,ts,0))*pop,'LineWidth',3,'Color',re); hold on
p4=plot(ts,(integrate(gI2s,ts,0)-integrate(gI2s_no,ts,0))*pop,'LineWidth',3,'Color',re,'LineStyle',':'); hold on
ylim([-23000 0])
title({'(H)'})

han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
%ylabel(han,'Quantity of Treatment', 'FontSize', 16);

xlabel(han,'Time (months)', 'FontSize', 16);

%sgtitle({'Cumulative Difference in the Number of Infections','Compared to the No-Treatment Scenario for Different','Drug Allocation Rules when Capacity is 10%'},'Color','Black', 'FontSize', 18);
saveas(gcf,'CumInfectionsDrug_10percent.png'); hold off

end

%Figure A9: Cumulative Infections for Drugs for 15% capacity per state
if true
  
pop=1000000;


ts=ts05;

%Type of interpolation
fitApprox = fittype('pchipinterp');

load WorkspacePT_Gauss.mat
gI1s=fit(ts,(I1s15_Dah),fitApprox);
gI1s_opt=fit(ts,(I1s15_D),fitApprox);
gI2s=fit(ts,(I2s15_Dah),fitApprox);
gI2s_opt=fit(ts,(I2s15_D),fitApprox);

load workspaceNoControls.mat
gI1s_no=fit(ts,I1s_Dno1PS,fitApprox);
gI2s_no=fit(ts,I2s_Dno1PS,fitApprox);


fig=figure
subplot(2,4,1)
p1=plot(ts,integrate(gI1s_opt,ts,0)./integrate(gI1s_no,ts,0)-1,'LineWidth',3,'Color',bl); hold on
p3=plot(ts,integrate(gI1s,ts,0)./integrate(gI1s_no,ts,0)-1,'LineWidth',3,'Color',bl,'LineStyle',':'); hold on
p2=plot(ts,integrate(gI2s_opt,ts,0)./integrate(gI2s_no,ts,0)-1,'LineWidth',3,'Color',re); hold on
p4=plot(ts,integrate(gI2s,ts,0)./integrate(gI2s_no,ts,0)-1,'LineWidth',3,'Color',re,'LineStyle',':'); hold on
ylim([-0.3 0])
title({'Perm. Immunity &','Compliance to TR','(A)'})
ylabel({'Relative','Difference'}, 'FontSize', 16); 
    legend1=legend([p1 p3 p2 p4],{'State 1: Optimal~~~','State 1: Ad Hoc~~~','State 2: Optimal~~~','State 2: Ad Hoc'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
    'Position',[0.117340711702374 0.479706609664701 0.811608069283622 0.0352380951245625],...
    'Orientation','horizontal',...
    'Interpreter','latex');


subplot(2,4,5)
p1=plot(ts,(integrate(gI1s_opt,ts,0)-integrate(gI1s_no,ts,0))*pop,'LineWidth',3,'Color',bl); hold on
p3=plot(ts,(integrate(gI1s,ts,0)-integrate(gI1s_no,ts,0))*pop,'LineWidth',3,'Color',bl,'LineStyle',':'); hold on
p2=plot(ts,(integrate(gI2s_opt,ts,0)-integrate(gI2s_no,ts,0))*pop,'LineWidth',3,'Color',re); hold on
p4=plot(ts,(integrate(gI2s,ts,0)-integrate(gI2s_no,ts,0))*pop,'LineWidth',3,'Color',re,'LineStyle',':'); hold on
ylim([-25000 0])
ylabel({'Absolute Difference','per 1M People'}, 'FontSize', 16); 
title({'(E)'})



%%%
load WorkspacePN_Gauss.mat
gI1s=fit(ts,(I1s15_Dah),fitApprox);
gI1s_opt=fit(ts,(I1s15_D),fitApprox);
gI2s=fit(ts,(I2s15_Dah),fitApprox);
gI2s_opt=fit(ts,(I2s15_D),fitApprox);

load workspaceNoControls.mat
gI1s_no=fit(ts,I1s_Dno3PN,fitApprox);
gI2s_no=fit(ts,I2s_Dno3PN,fitApprox);

subplot(2,4,2)
p1=plot(ts,integrate(gI1s_opt,ts,0)./integrate(gI1s_no,ts,0)-1,'LineWidth',3,'Color',bl); hold on
p3=plot(ts,integrate(gI1s,ts,0)./integrate(gI1s_no,ts,0)-1,'LineWidth',3,'Color',bl,'LineStyle',':'); hold on
p2=plot(ts,integrate(gI2s_opt,ts,0)./integrate(gI2s_no,ts,0)-1,'LineWidth',3,'Color',re); hold on
p4=plot(ts,integrate(gI2s,ts,0)./integrate(gI2s_no,ts,0)-1,'LineWidth',3,'Color',re,'LineStyle',':'); hold on
title({'Perm. Immunity &','Noncompliance to TR','(B)'})
ylim([-0.3 0])

subplot(2,4,6)
p1=plot(ts,(integrate(gI1s_opt,ts,0)-integrate(gI1s_no,ts,0))*pop,'LineWidth',3,'Color',bl); hold on
p3=plot(ts,(integrate(gI1s,ts,0)-integrate(gI1s_no,ts,0))*pop,'LineWidth',3,'Color',bl,'LineStyle',':'); hold on
p2=plot(ts,(integrate(gI2s_opt,ts,0)-integrate(gI2s_no,ts,0))*pop,'LineWidth',3,'Color',re); hold on
p4=plot(ts,(integrate(gI2s,ts,0)-integrate(gI2s_no,ts,0))*pop,'LineWidth',3,'Color',re,'LineStyle',':'); hold on
ylim([-25000 0])
title({'(F)'})

%%%
load WorkspaceST_Gauss.mat
gI1s=fit(ts,(I1s15_Dah),fitApprox);
gI1s_opt=fit(ts,(I1s15_D),fitApprox);
gI2s=fit(ts,(I2s15_Dah),fitApprox);
gI2s_opt=fit(ts,(I2s15_D),fitApprox);

load workspaceNoControls.mat
gI1s_no=fit(ts,I1s_Dno2SS,fitApprox);
gI2s_no=fit(ts,I2s_Dno2SS,fitApprox);

subplot(2,4,3)
p1=plot(ts,integrate(gI1s_opt,ts,0)./integrate(gI1s_no,ts,0)-1,'LineWidth',3,'Color',bl); hold on
p3=plot(ts,integrate(gI1s,ts,0)./integrate(gI1s_no,ts,0)-1,'LineWidth',3,'Color',bl,'LineStyle',':'); hold on
p2=plot(ts,integrate(gI2s_opt,ts,0)./integrate(gI2s_no,ts,0)-1,'LineWidth',3,'Color',re); hold on
p4=plot(ts,integrate(gI2s,ts,0)./integrate(gI2s_no,ts,0)-1,'LineWidth',3,'Color',re,'LineStyle',':'); hold on
title({'6-mo. Immunity &','Compliance to TR','(C)'})
ylim([-0.3 0])

subplot(2,4,7)
p1=plot(ts,(integrate(gI1s_opt,ts,0)-integrate(gI1s_no,ts,0))*pop,'LineWidth',3,'Color',bl); hold on
p3=plot(ts,(integrate(gI1s,ts,0)-integrate(gI1s_no,ts,0))*pop,'LineWidth',3,'Color',bl,'LineStyle',':'); hold on
p2=plot(ts,(integrate(gI2s_opt,ts,0)-integrate(gI2s_no,ts,0))*pop,'LineWidth',3,'Color',re); hold on
p4=plot(ts,(integrate(gI2s,ts,0)-integrate(gI2s_no,ts,0))*pop,'LineWidth',3,'Color',re,'LineStyle',':'); hold on
ylim([-25000 0])
title({'(G)'})
%%%
load WorkspaceSN_Gauss.mat
gI1s=fit(ts,(I1s15_Dah),fitApprox);
gI1s_opt=fit(ts,(I1s15_D),fitApprox);
gI2s=fit(ts,(I2s15_Dah),fitApprox);
gI2s_opt=fit(ts,(I2s15_D),fitApprox);

load workspaceNoControls.mat
gI1s_no=fit(ts,I1s_Dno4SN,fitApprox);
gI2s_no=fit(ts,I2s_Dno4SN,fitApprox);

subplot(2,4,4)
p1=plot(ts,integrate(gI1s_opt,ts,0)./integrate(gI1s_no,ts,0)-1,'LineWidth',3,'Color',bl); hold on
p3=plot(ts,integrate(gI1s,ts,0)./integrate(gI1s_no,ts,0)-1,'LineWidth',3,'Color',bl,'LineStyle',':'); hold on
p2=plot(ts,integrate(gI2s_opt,ts,0)./integrate(gI2s_no,ts,0)-1,'LineWidth',3,'Color',re); hold on
p4=plot(ts,integrate(gI2s,ts,0)./integrate(gI2s_no,ts,0)-1,'LineWidth',3,'Color',re,'LineStyle',':'); hold on
title({'6-mo. Immunity &','Noncompliance to TR','(D)'})
ylim([-0.3 0])

subplot(2,4,8)
p1=plot(ts,(integrate(gI1s_opt,ts,0)-integrate(gI1s_no,ts,0))*pop,'LineWidth',3,'Color',bl); hold on
p3=plot(ts,(integrate(gI1s,ts,0)-integrate(gI1s_no,ts,0))*pop,'LineWidth',3,'Color',bl,'LineStyle',':'); hold on
p2=plot(ts,(integrate(gI2s_opt,ts,0)-integrate(gI2s_no,ts,0))*pop,'LineWidth',3,'Color',re); hold on
p4=plot(ts,(integrate(gI2s,ts,0)-integrate(gI2s_no,ts,0))*pop,'LineWidth',3,'Color',re,'LineStyle',':'); hold on
ylim([-25000 0])
title({'(H)'})

han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
%ylabel(han,'Quantity of Treatment', 'FontSize', 16);

xlabel(han,'Time (months)', 'FontSize', 16);

%sgtitle({'Cumulative Difference in the Number of Infections','Compared to the No-Treatment Scenario for Different','Drug Allocation Rules when Capacity is 15%'},'Color','Black', 'FontSize', 18);
saveas(gcf,'CumInfectionsDrug_15percent.png'); hold off

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Vaccines %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 %Figure A21: Cumulative Infections for Vaccines for 5% capacity per state
if true
  
pop=1000000;


ts=ts05;

%Type of interpolation
fitApprox = fittype('pchipinterp');

load WorkspacePT_Gauss.mat
gI1s=    fit(ts,(I1s05_ah),fitApprox);
gI1s_opt=fit(ts,(I1s05),fitApprox);
gI2s=    fit(ts,(I2s05_ah),fitApprox);
gI2s_opt=fit(ts,(I2s05),fitApprox);

load workspaceNoControls.mat
gI1s_no=fit(ts,I1s_no1PS,fitApprox);
gI2s_no=fit(ts,I2s_no1PS,fitApprox);


fig=figure
subplot(2,4,1)
p1=plot(ts,integrate(gI1s_opt,ts,0)./integrate(gI1s_no,ts,0)-1,'LineWidth',3,'Color',bl); hold on
p3=plot(ts,integrate(gI1s,ts,0)./integrate(gI1s_no,ts,0)-1,'LineWidth',3,'Color',bl,'LineStyle',':'); hold on
p2=plot(ts,integrate(gI2s_opt,ts,0)./integrate(gI2s_no,ts,0)-1,'LineWidth',3,'Color',re); hold on
p4=plot(ts,integrate(gI2s,ts,0)./integrate(gI2s_no,ts,0)-1,'LineWidth',3,'Color',re,'LineStyle',':'); hold on
ylim([-0.2 0])
title({'Perm. Immunity &','Compliance to TR','(A)'})
ylabel({'Relative','Difference'}, 'FontSize', 16); 
    legend1=legend([p1 p3 p2 p4],{'State 1: Optimal~~~','State 1: Ad Hoc~~~','State 2: Optimal~~~','State 2: Ad Hoc'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
    'Position',[0.117340711702374 0.479706609664701 0.811608069283622 0.0352380951245625],...
    'Orientation','horizontal',...
    'Interpreter','latex');


subplot(2,4,5)
p1=plot(ts,(integrate(gI1s_opt,ts,0)-integrate(gI1s_no,ts,0))*pop,'LineWidth',3,'Color',bl); hold on
p3=plot(ts,(integrate(gI1s,ts,0)-integrate(gI1s_no,ts,0))*pop,'LineWidth',3,'Color',bl,'LineStyle',':'); hold on
p2=plot(ts,(integrate(gI2s_opt,ts,0)-integrate(gI2s_no,ts,0))*pop,'LineWidth',3,'Color',re); hold on
p4=plot(ts,(integrate(gI2s,ts,0)-integrate(gI2s_no,ts,0))*pop,'LineWidth',3,'Color',re,'LineStyle',':'); hold on
ylim([-17000 0])
ylabel({'Absolute Difference','per 1M People'}, 'FontSize', 16); 
title({'(E)'})



%%%
load WorkspacePN_Gauss.mat
gI1s=fit(ts,(I1s05_ah),fitApprox);
gI1s_opt=fit(ts,(I1s05),fitApprox);
gI2s=fit(ts,(I2s05_ah),fitApprox);
gI2s_opt=fit(ts,(I2s05),fitApprox);

load workspaceNoControls.mat
gI1s_no=fit(ts,I1s_no3PN,fitApprox);
gI2s_no=fit(ts,I2s_no3PN,fitApprox);

subplot(2,4,2)
p1=plot(ts,integrate(gI1s_opt,ts,0)./integrate(gI1s_no,ts,0)-1,'LineWidth',3,'Color',bl); hold on
p3=plot(ts,integrate(gI1s,ts,0)./integrate(gI1s_no,ts,0)-1,'LineWidth',3,'Color',bl,'LineStyle',':'); hold on
p2=plot(ts,integrate(gI2s_opt,ts,0)./integrate(gI2s_no,ts,0)-1,'LineWidth',3,'Color',re); hold on
p4=plot(ts,integrate(gI2s,ts,0)./integrate(gI2s_no,ts,0)-1,'LineWidth',3,'Color',re,'LineStyle',':'); hold on
title({'Perm. Immunity &','Noncompliance to TR','(B)'})
ylim([-0.2 0])


subplot(2,4,6)
p1=plot(ts,(integrate(gI1s_opt,ts,0)-integrate(gI1s_no,ts,0))*pop,'LineWidth',3,'Color',bl); hold on
p3=plot(ts,(integrate(gI1s,ts,0)-integrate(gI1s_no,ts,0))*pop,'LineWidth',3,'Color',bl,'LineStyle',':'); hold on
p2=plot(ts,(integrate(gI2s_opt,ts,0)-integrate(gI2s_no,ts,0))*pop,'LineWidth',3,'Color',re); hold on
p4=plot(ts,(integrate(gI2s,ts,0)-integrate(gI2s_no,ts,0))*pop,'LineWidth',3,'Color',re,'LineStyle',':'); hold on
ylim([-17000 0])
title({'(F)'})

%%%
load WorkspaceST_Gauss.mat
gI1s=fit(ts,(I1s05_ah),fitApprox);
gI1s_opt=fit(ts,(I1s05),fitApprox);
gI2s=fit(ts,(I2s05_ah),fitApprox);
gI2s_opt=fit(ts,(I2s05),fitApprox);

load workspaceNoControls.mat
gI1s_no=fit(ts,I1s_no2SS,fitApprox);
gI2s_no=fit(ts,I2s_no2SS,fitApprox);

subplot(2,4,3)
p1=plot(ts,integrate(gI1s_opt,ts,0)./integrate(gI1s_no,ts,0)-1,'LineWidth',3,'Color',bl); hold on
p3=plot(ts,integrate(gI1s,ts,0)./integrate(gI1s_no,ts,0)-1,'LineWidth',3,'Color',bl,'LineStyle',':'); hold on
p2=plot(ts,integrate(gI2s_opt,ts,0)./integrate(gI2s_no,ts,0)-1,'LineWidth',3,'Color',re); hold on
p4=plot(ts,integrate(gI2s,ts,0)./integrate(gI2s_no,ts,0)-1,'LineWidth',3,'Color',re,'LineStyle',':'); hold on
title({'6-mo. Immunity &','Compliance to TR','(C)'})
ylim([-0.2 0])

subplot(2,4,7)
p1=plot(ts,(integrate(gI1s_opt,ts,0)-integrate(gI1s_no,ts,0))*pop,'LineWidth',3,'Color',bl); hold on
p3=plot(ts,(integrate(gI1s,ts,0)-integrate(gI1s_no,ts,0))*pop,'LineWidth',3,'Color',bl,'LineStyle',':'); hold on
p2=plot(ts,(integrate(gI2s_opt,ts,0)-integrate(gI2s_no,ts,0))*pop,'LineWidth',3,'Color',re); hold on
p4=plot(ts,(integrate(gI2s,ts,0)-integrate(gI2s_no,ts,0))*pop,'LineWidth',3,'Color',re,'LineStyle',':'); hold on
ylim([-17000 0])
title({'(G)'})
%%%
load WorkspaceSN_Gauss.mat
gI1s=fit(ts,(I1s05_ah),fitApprox);
gI1s_opt=fit(ts,(I1s05),fitApprox);
gI2s=fit(ts,(I2s05_ah),fitApprox);
gI2s_opt=fit(ts,(I2s05),fitApprox);

load workspaceNoControls.mat
gI1s_no=fit(ts,I1s_no4SN,fitApprox);
gI2s_no=fit(ts,I2s_no4SN,fitApprox);

subplot(2,4,4)
p1=plot(ts,integrate(gI1s_opt,ts,0)./integrate(gI1s_no,ts,0)-1,'LineWidth',3,'Color',bl); hold on
p3=plot(ts,integrate(gI1s,ts,0)./integrate(gI1s_no,ts,0)-1,'LineWidth',3,'Color',bl,'LineStyle',':'); hold on
p2=plot(ts,integrate(gI2s_opt,ts,0)./integrate(gI2s_no,ts,0)-1,'LineWidth',3,'Color',re); hold on
p4=plot(ts,integrate(gI2s,ts,0)./integrate(gI2s_no,ts,0)-1,'LineWidth',3,'Color',re,'LineStyle',':'); hold on
title({'6-mo. Immunity &','Noncompliance to TR','(D)'})
ylim([-0.2 0])


subplot(2,4,8)
p1=plot(ts,(integrate(gI1s_opt,ts,0)-integrate(gI1s_no,ts,0))*pop,'LineWidth',3,'Color',bl); hold on
p3=plot(ts,(integrate(gI1s,ts,0)-integrate(gI1s_no,ts,0))*pop,'LineWidth',3,'Color',bl,'LineStyle',':'); hold on
p2=plot(ts,(integrate(gI2s_opt,ts,0)-integrate(gI2s_no,ts,0))*pop,'LineWidth',3,'Color',re); hold on
p4=plot(ts,(integrate(gI2s,ts,0)-integrate(gI2s_no,ts,0))*pop,'LineWidth',3,'Color',re,'LineStyle',':'); hold on
ylim([-17000 0])
title({'(H)'})

han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
%ylabel(han,'Quantity of Treatment', 'FontSize', 16);

xlabel(han,'Time (months)', 'FontSize', 16);

%sgtitle({'Cumulative Difference in the Number of Infections','Compared to the No-Treatment Scenario for Different','Vaccine Allocation Rules when Capacity is 5%'},'Color','Black', 'FontSize', 18);
saveas(gcf,'CumInfectionsVaccine_5percent.png'); hold off

end

%Figure 7: Cumulative Infections for Vaccines for 15% capacity per state
if true
  
pop=1000000;




%Type of interpolation
fitApprox = fittype('pchipinterp');

load WorkspacePT_Gauss.mat
ts=ts05;
gI1s=    fit(ts,(I1s10_ah),fitApprox);
gI1s_opt=fit(ts,(I1s10),fitApprox);
gI2s=    fit(ts,(I2s10_ah),fitApprox);
gI2s_opt=fit(ts,(I2s10),fitApprox);

load workspaceNoControls.mat
gI1s_no=fit(ts,I1s_no1PS,fitApprox);
gI2s_no=fit(ts,I2s_no1PS,fitApprox);


fig=figure
subplot(2,4,1)
p1=plot(ts,integrate(gI1s_opt,ts,0)./integrate(gI1s_no,ts,0)-1,'LineWidth',3,'Color',bl); hold on
p3=plot(ts,integrate(gI1s,ts,0)./integrate(gI1s_no,ts,0)-1,'LineWidth',3,'Color',bl,'LineStyle',':'); hold on
p2=plot(ts,integrate(gI2s_opt,ts,0)./integrate(gI2s_no,ts,0)-1,'LineWidth',3,'Color',re); hold on
p4=plot(ts,integrate(gI2s,ts,0)./integrate(gI2s_no,ts,0)-1,'LineWidth',3,'Color',re,'LineStyle',':'); hold on
ylim([-0.3 0])
title({'Perm. Immunity &','Compliance to TR','(A)'})
ylabel({'Relative','Difference'}, 'FontSize', 16); 
    legend1=legend([p1 p3 p2 p4],{'State 1: Optimal~~~','State 1: Ad Hoc~~~','State 2: Optimal~~~','State 2: Ad Hoc'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
    'Position',[0.117340711702374 0.479706609664701 0.811608069283622 0.0352380951245625],...
    'Orientation','horizontal',...
    'Interpreter','latex');


subplot(2,4,5)
p1=plot(ts,(integrate(gI1s_opt,ts,0)-integrate(gI1s_no,ts,0))*pop,'LineWidth',3,'Color',bl); hold on
p3=plot(ts,(integrate(gI1s,ts,0)-integrate(gI1s_no,ts,0))*pop,'LineWidth',3,'Color',bl,'LineStyle',':'); hold on
p2=plot(ts,(integrate(gI2s_opt,ts,0)-integrate(gI2s_no,ts,0))*pop,'LineWidth',3,'Color',re); hold on
p4=plot(ts,(integrate(gI2s,ts,0)-integrate(gI2s_no,ts,0))*pop,'LineWidth',3,'Color',re,'LineStyle',':'); hold on
ylim([-26000 0])
ylabel({'Absolute Difference','per 1M People'}, 'FontSize', 16); 
title({'(E)'})



%%%
load WorkspacePN_Gauss.mat
gI1s=    fit(ts,(I1s10_ah),fitApprox);
gI1s_opt=fit(ts,(I1s10),fitApprox);
gI2s=    fit(ts,(I2s10_ah),fitApprox);
gI2s_opt=fit(ts,(I2s10),fitApprox);

load workspaceNoControls.mat
gI1s_no=fit(ts,I1s_no3PN,fitApprox);
gI2s_no=fit(ts,I2s_no3PN,fitApprox);

subplot(2,4,2)
p1=plot(ts,integrate(gI1s_opt,ts,0)./integrate(gI1s_no,ts,0)-1,'LineWidth',3,'Color',bl); hold on
p3=plot(ts,integrate(gI1s,ts,0)./integrate(gI1s_no,ts,0)-1,'LineWidth',3,'Color',bl,'LineStyle',':'); hold on
p2=plot(ts,integrate(gI2s_opt,ts,0)./integrate(gI2s_no,ts,0)-1,'LineWidth',3,'Color',re); hold on
p4=plot(ts,integrate(gI2s,ts,0)./integrate(gI2s_no,ts,0)-1,'LineWidth',3,'Color',re,'LineStyle',':'); hold on
title({'Perm. Immunity &','Noncompliance to TR','(B)'})
ylim([-0.3 0])

subplot(2,4,6)
p1=plot(ts,(integrate(gI1s_opt,ts,0)-integrate(gI1s_no,ts,0))*pop,'LineWidth',3,'Color',bl); hold on
p3=plot(ts,(integrate(gI1s,ts,0)-integrate(gI1s_no,ts,0))*pop,'LineWidth',3,'Color',bl,'LineStyle',':'); hold on
p2=plot(ts,(integrate(gI2s_opt,ts,0)-integrate(gI2s_no,ts,0))*pop,'LineWidth',3,'Color',re); hold on
p4=plot(ts,(integrate(gI2s,ts,0)-integrate(gI2s_no,ts,0))*pop,'LineWidth',3,'Color',re,'LineStyle',':'); hold on
ylim([-26000 0])
title({'(F)'})

%%%
load WorkspaceST_Gauss.mat
gI1s=    fit(ts,(I1s10_ah),fitApprox);
gI1s_opt=fit(ts,(I1s10),fitApprox);
gI2s=    fit(ts,(I2s10_ah),fitApprox);
gI2s_opt=fit(ts,(I2s10),fitApprox);

load workspaceNoControls.mat
gI1s_no=fit(ts,I1s_no2SS,fitApprox);
gI2s_no=fit(ts,I2s_no2SS,fitApprox);

subplot(2,4,3)
p1=plot(ts,integrate(gI1s_opt,ts,0)./integrate(gI1s_no,ts,0)-1,'LineWidth',3,'Color',bl); hold on
p3=plot(ts,integrate(gI1s,ts,0)./integrate(gI1s_no,ts,0)-1,'LineWidth',3,'Color',bl,'LineStyle',':'); hold on
p2=plot(ts,integrate(gI2s_opt,ts,0)./integrate(gI2s_no,ts,0)-1,'LineWidth',3,'Color',re); hold on
p4=plot(ts,integrate(gI2s,ts,0)./integrate(gI2s_no,ts,0)-1,'LineWidth',3,'Color',re,'LineStyle',':'); hold on
title({'6-mo. Immunity &','Compliance to TR','(C)'})
ylim([-0.3 0])

subplot(2,4,7)
p1=plot(ts,(integrate(gI1s_opt,ts,0)-integrate(gI1s_no,ts,0))*pop,'LineWidth',3,'Color',bl); hold on
p3=plot(ts,(integrate(gI1s,ts,0)-integrate(gI1s_no,ts,0))*pop,'LineWidth',3,'Color',bl,'LineStyle',':'); hold on
p2=plot(ts,(integrate(gI2s_opt,ts,0)-integrate(gI2s_no,ts,0))*pop,'LineWidth',3,'Color',re); hold on
p4=plot(ts,(integrate(gI2s,ts,0)-integrate(gI2s_no,ts,0))*pop,'LineWidth',3,'Color',re,'LineStyle',':'); hold on
ylim([-26000 0])
title({'(G)'})

%%%
load WorkspaceSN_Gauss.mat
gI1s=    fit(ts,(I1s10_ah),fitApprox);
gI1s_opt=fit(ts,(I1s10),fitApprox);
gI2s=    fit(ts,(I2s10_ah),fitApprox);
gI2s_opt=fit(ts,(I2s10),fitApprox);

load workspaceNoControls.mat
gI1s_no=fit(ts,I1s_no4SN,fitApprox);
gI2s_no=fit(ts,I2s_no4SN,fitApprox);

subplot(2,4,4)
p1=plot(ts,integrate(gI1s_opt,ts,0)./integrate(gI1s_no,ts,0)-1,'LineWidth',3,'Color',bl); hold on
p3=plot(ts,integrate(gI1s,ts,0)./integrate(gI1s_no,ts,0)-1,'LineWidth',3,'Color',bl,'LineStyle',':'); hold on
p2=plot(ts,integrate(gI2s_opt,ts,0)./integrate(gI2s_no,ts,0)-1,'LineWidth',3,'Color',re); hold on
p4=plot(ts,integrate(gI2s,ts,0)./integrate(gI2s_no,ts,0)-1,'LineWidth',3,'Color',re,'LineStyle',':'); hold on
title({'6-mo. Immunity &','Noncompliance to TR','(D)'})
ylim([-0.3 0])


subplot(2,4,8)
p1=plot(ts,(integrate(gI1s_opt,ts,0)-integrate(gI1s_no,ts,0))*pop,'LineWidth',3,'Color',bl); hold on
p3=plot(ts,(integrate(gI1s,ts,0)-integrate(gI1s_no,ts,0))*pop,'LineWidth',3,'Color',bl,'LineStyle',':'); hold on
p2=plot(ts,(integrate(gI2s_opt,ts,0)-integrate(gI2s_no,ts,0))*pop,'LineWidth',3,'Color',re); hold on
p4=plot(ts,(integrate(gI2s,ts,0)-integrate(gI2s_no,ts,0))*pop,'LineWidth',3,'Color',re,'LineStyle',':'); hold on
ylim([-26000 0])
title({'(H)'})

han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
%ylabel(han,'Quantity of Treatment', 'FontSize', 16);

xlabel(han,'Time (months)', 'FontSize', 16);

%sgtitle({'Cumulative Difference in the Number of Infections','Compared to the No-Treatment Scenario for Different','Vaccine Allocation Rules when Capacity is 10%'},'Color','Black', 'FontSize', 18);
saveas(gcf,'CumInfectionsVaccine_10percent.png'); hold off

end

 %Figure A22: Cumulative Infections for Vaccines for 10% capacity per state
if true
  
pop=1000000;


ts=ts05;

%Type of interpolation
fitApprox = fittype('pchipinterp');

load WorkspacePT_Gauss.mat
gI1s=    fit(ts,(I1s15_ah),fitApprox);
gI1s_opt=fit(ts,(I1s15),fitApprox);
gI2s=    fit(ts,(I2s15_ah),fitApprox);
gI2s_opt=fit(ts,(I2s15),fitApprox);

load workspaceNoControls.mat
gI1s_no=fit(ts,I1s_no1PS,fitApprox);
gI2s_no=fit(ts,I2s_no1PS,fitApprox);


fig=figure
subplot(2,4,1)
p1=plot(ts,integrate(gI1s_opt,ts,0)./integrate(gI1s_no,ts,0)-1,'LineWidth',3,'Color',bl); hold on
p3=plot(ts,integrate(gI1s,ts,0)./integrate(gI1s_no,ts,0)-1,'LineWidth',3,'Color',bl,'LineStyle',':'); hold on
p2=plot(ts,integrate(gI2s_opt,ts,0)./integrate(gI2s_no,ts,0)-1,'LineWidth',3,'Color',re); hold on
p4=plot(ts,integrate(gI2s,ts,0)./integrate(gI2s_no,ts,0)-1,'LineWidth',3,'Color',re,'LineStyle',':'); hold on
ylim([-0.4 0])
title({'Perm. Immunity &','Compliance to TR','(A)'})
ylabel({'Relative','Difference'}, 'FontSize', 16); 
    legend1=legend([p1 p3 p2 p4],{'State 1: Optimal~~~','State 1: Ad Hoc~~~','State 2: Optimal~~~','State 2: Ad Hoc'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
    'Position',[0.117340711702374 0.479706609664701 0.811608069283622 0.0352380951245625],...
    'Orientation','horizontal',...
    'Interpreter','latex');


subplot(2,4,5)
p1=plot(ts,(integrate(gI1s_opt,ts,0)-integrate(gI1s_no,ts,0))*pop,'LineWidth',3,'Color',bl); hold on
p3=plot(ts,(integrate(gI1s,ts,0)-integrate(gI1s_no,ts,0))*pop,'LineWidth',3,'Color',bl,'LineStyle',':'); hold on
p2=plot(ts,(integrate(gI2s_opt,ts,0)-integrate(gI2s_no,ts,0))*pop,'LineWidth',3,'Color',re); hold on
p4=plot(ts,(integrate(gI2s,ts,0)-integrate(gI2s_no,ts,0))*pop,'LineWidth',3,'Color',re,'LineStyle',':'); hold on
ylim([-34000 0])
ylabel({'Absolute Difference','per 1M People'}, 'FontSize', 16); 
title({'(E)'})



%%%
load WorkspacePN_Gauss.mat
gI1s=    fit(ts,(I1s15_ah),fitApprox);
gI1s_opt=fit(ts,(I1s15),fitApprox);
gI2s=    fit(ts,(I2s15_ah),fitApprox);
gI2s_opt=fit(ts,(I2s15),fitApprox);

load workspaceNoControls.mat
gI1s_no=fit(ts,I1s_no3PN,fitApprox);
gI2s_no=fit(ts,I2s_no3PN,fitApprox);

subplot(2,4,2)
p1=plot(ts,integrate(gI1s_opt,ts,0)./integrate(gI1s_no,ts,0)-1,'LineWidth',3,'Color',bl); hold on
p3=plot(ts,integrate(gI1s,ts,0)./integrate(gI1s_no,ts,0)-1,'LineWidth',3,'Color',bl,'LineStyle',':'); hold on
p2=plot(ts,integrate(gI2s_opt,ts,0)./integrate(gI2s_no,ts,0)-1,'LineWidth',3,'Color',re); hold on
p4=plot(ts,integrate(gI2s,ts,0)./integrate(gI2s_no,ts,0)-1,'LineWidth',3,'Color',re,'LineStyle',':'); hold on
title({'Perm. Immunity &','Noncompliance to TR','(B)'})
ylim([-0.4 0])

subplot(2,4,6)
p1=plot(ts,(integrate(gI1s_opt,ts,0)-integrate(gI1s_no,ts,0))*pop,'LineWidth',3,'Color',bl); hold on
p3=plot(ts,(integrate(gI1s,ts,0)-integrate(gI1s_no,ts,0))*pop,'LineWidth',3,'Color',bl,'LineStyle',':'); hold on
p2=plot(ts,(integrate(gI2s_opt,ts,0)-integrate(gI2s_no,ts,0))*pop,'LineWidth',3,'Color',re); hold on
p4=plot(ts,(integrate(gI2s,ts,0)-integrate(gI2s_no,ts,0))*pop,'LineWidth',3,'Color',re,'LineStyle',':'); hold on
ylim([-34000 0])
title({'(F)'})

%%%
load WorkspaceST_Gauss.mat
gI1s=    fit(ts,(I1s15_ah),fitApprox);
gI1s_opt=fit(ts,(I1s15),fitApprox);
gI2s=    fit(ts,(I2s15_ah),fitApprox);
gI2s_opt=fit(ts,(I2s15),fitApprox);

load workspaceNoControls.mat
gI1s_no=fit(ts,I1s_no2SS,fitApprox);
gI2s_no=fit(ts,I2s_no2SS,fitApprox);

subplot(2,4,3)
p1=plot(ts,integrate(gI1s_opt,ts,0)./integrate(gI1s_no,ts,0)-1,'LineWidth',3,'Color',bl); hold on
p3=plot(ts,integrate(gI1s,ts,0)./integrate(gI1s_no,ts,0)-1,'LineWidth',3,'Color',bl,'LineStyle',':'); hold on
p2=plot(ts,integrate(gI2s_opt,ts,0)./integrate(gI2s_no,ts,0)-1,'LineWidth',3,'Color',re); hold on
p4=plot(ts,integrate(gI2s,ts,0)./integrate(gI2s_no,ts,0)-1,'LineWidth',3,'Color',re,'LineStyle',':'); hold on
title({'6-mo. Immunity &','Compliance to TR','(C)'})
ylim([-0.4 0])

subplot(2,4,7)
p1=plot(ts,(integrate(gI1s_opt,ts,0)-integrate(gI1s_no,ts,0))*pop,'LineWidth',3,'Color',bl); hold on
p3=plot(ts,(integrate(gI1s,ts,0)-integrate(gI1s_no,ts,0))*pop,'LineWidth',3,'Color',bl,'LineStyle',':'); hold on
p2=plot(ts,(integrate(gI2s_opt,ts,0)-integrate(gI2s_no,ts,0))*pop,'LineWidth',3,'Color',re); hold on
p4=plot(ts,(integrate(gI2s,ts,0)-integrate(gI2s_no,ts,0))*pop,'LineWidth',3,'Color',re,'LineStyle',':'); hold on
ylim([-34000 0])
title({'(G)'})
%%%
load WorkspaceSN_Gauss.mat
gI1s=    fit(ts,(I1s15_ah),fitApprox);
gI1s_opt=fit(ts,(I1s15),fitApprox);
gI2s=    fit(ts,(I2s15_ah),fitApprox);
gI2s_opt=fit(ts,(I2s15),fitApprox);

load workspaceNoControls.mat
gI1s_no=fit(ts,I1s_no4SN,fitApprox);
gI2s_no=fit(ts,I2s_no4SN,fitApprox);

subplot(2,4,4)
p1=plot(ts,integrate(gI1s_opt,ts,0)./integrate(gI1s_no,ts,0)-1,'LineWidth',3,'Color',bl); hold on
p3=plot(ts,integrate(gI1s,ts,0)./integrate(gI1s_no,ts,0)-1,'LineWidth',3,'Color',bl,'LineStyle',':'); hold on
p2=plot(ts,integrate(gI2s_opt,ts,0)./integrate(gI2s_no,ts,0)-1,'LineWidth',3,'Color',re); hold on
p4=plot(ts,integrate(gI2s,ts,0)./integrate(gI2s_no,ts,0)-1,'LineWidth',3,'Color',re,'LineStyle',':'); hold on
title({'6-mo. Immunity &','Noncompliance to TR','(D)'})
ylim([-0.4 0])


subplot(2,4,8)
p1=plot(ts,(integrate(gI1s_opt,ts,0)-integrate(gI1s_no,ts,0))*pop,'LineWidth',3,'Color',bl); hold on
p3=plot(ts,(integrate(gI1s,ts,0)-integrate(gI1s_no,ts,0))*pop,'LineWidth',3,'Color',bl,'LineStyle',':'); hold on
p2=plot(ts,(integrate(gI2s_opt,ts,0)-integrate(gI2s_no,ts,0))*pop,'LineWidth',3,'Color',re); hold on
p4=plot(ts,(integrate(gI2s,ts,0)-integrate(gI2s_no,ts,0))*pop,'LineWidth',3,'Color',re,'LineStyle',':'); hold on
ylim([-34000 0])
title({'(H)'})

han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
%ylabel(han,'Quantity of Treatment', 'FontSize', 16);

xlabel(han,'Time (months)', 'FontSize', 16);

%sgtitle({'Cumulative Difference in the Number of Infections','Compared to the No-Treatment Scenario for Different','Vaccine Allocation Rules when Capacity is 15%'},'Color','Black', 'FontSize', 18);
saveas(gcf,'CumInfectionsVaccine_15percent.png'); hold off

end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%   Sensitivity Analyses  %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%   Figures for workability cost parameter  %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Figure A13: Sensitivity of drug results to workability cost parameter
if true

    
    fig=figure
    
load Workspace_Cadj_PT_Drug.mat CADJ DrugProp DrugInfection CostProp
    
subplot(241)
yyaxis left
    plot(log(CADJ),DrugProp,'LineWidth',3);hold on
    xticks([log(CADJ(1)),log(1e5),log(CADJ(8))]);
    xticklabels({'1e^1','1e^5','1e^9'});
    xlim([log(CADJ(1)) log(CADJ(end))])
    ylabel({'Variance of the Optimal ','Deviation From the Ad Hoc'},'Interpreter','tex','FontSize', 14);     
    ylim([0 0.000012])
    plot([log(1e4) log(1e4)],[0 100],'k:','LineWidth',2); 
    plot([log(1e7) log(1e7)],[0 100],'k-','LineWidth',2); 
    title({'Perm. Immunity &','Compliance to TR','(A)',''})
yyaxis right
    plot(log(CADJ),DrugInfection*100,'LineWidth',3);hold on   
    ylim([-2.5 0])
subplot(245)
    plot(log(CADJ),CostProp,'LineWidth',3);hold on
    xticks([log(CADJ(1)),log(1e5),log(CADJ(8))]);
    xticklabels({'1e^1','1e^5','1e^9'});
    ylabel({'$\frac{\textrm{Total Workability}}{\textrm{Total Vaccine}}$'},'Interpreter','latex', 'FontSize', 14);
    ylim([0 25])
    plot([log(1e4) log(1e4)],[0 100],'k:','LineWidth',2); 
    plot([log(1e7) log(1e7)],[0 100],'k-','LineWidth',2); 
    title({'(E)'})
    xlim([log(CADJ(1)) log(CADJ(end))])
    load Workspace_Cadj_PN_Drug.mat CADJ DrugProp DrugInfection CostProp
    
subplot(242)
yyaxis left
    plot(log(CADJ),DrugProp,'LineWidth',3);hold on
    xticks([log(CADJ(1)),log(CADJ(5)),log(CADJ(8))]);
    xticklabels({'1e^1','1e^5','1e^9'});
    xlim([log(CADJ(1)) log(CADJ(end))])
    plot([log(1e4) log(1e4)],[0 100],'k:','LineWidth',2); 
    plot([log(1e7) log(1e7)],[0 100],'k-','LineWidth',2); 
    ylim([0 0.000012])
    title({'Perm. Immunity &','Noncompliance to TR','(B)',''})
        
yyaxis right
    plot(log(CADJ),DrugInfection*100,'LineWidth',3);hold on  
    ylim([-0.025 0])
subplot(246)
    plot(log(CADJ),CostProp,'LineWidth',3);hold on
    xticks([log(CADJ(1)),log(CADJ(5)),log(CADJ(8))]);
    xticklabels({'1e^1','1e^5','1e^9'});
    ylim([0 0.25])  
    plot([log(1e4) log(1e4)],[0 100],'k:','LineWidth',2); 
    plot([log(1e7) log(1e7)],[0 100],'k-','LineWidth',2); 
    title({'(F)'})
    xlim([log(CADJ(1)) log(CADJ(end))])    
load Workspace_Cadj_ST_Drug.mat CADJ DrugProp DrugInfection CostProp

subplot(243)
yyaxis left
    plot(log(CADJ),DrugProp,'LineWidth',3);hold on
    xticks([log(CADJ(1)),log(1e5),log(CADJ(8))]);
    xticklabels({'1e^1','1e^5','1e^9'});
    xlim([log(CADJ(1)) log(CADJ(end))])
    plot([log(1e4) log(1e4)],[0 100],'k:','LineWidth',2); 
    plot([log(1e7) log(1e7)],[0 100],'k-','LineWidth',2);  
    ylim([0 0.000012])
    title({'6-mo. Immunity &','Compliance to TR','(C)',''})
yyaxis right
    plot(log(CADJ),DrugInfection*100,'LineWidth',3);hold on  
    ylim([-2.5 0])
subplot(247)
    plot(log(CADJ),CostProp,'LineWidth',3);hold on
    xticks([log(CADJ(1)),log(CADJ(5)),log(CADJ(8))]);
    xticklabels({'1e^1','1e^5','1e^9'});
    ylim([0 25])
    plot([log(1e4) log(1e4)],[0 100],'k:','LineWidth',2); 
    plot([log(1e7) log(1e7)],[0 100],'k-','LineWidth',2); 
    title({'(G)'})    
    xlim([log(CADJ(1)) log(CADJ(end))])
        
 load Workspace_Cadj_SN_Drug.mat CADJ DrugProp DrugInfection CostProp

subplot(244)
yyaxis left
    plot(log(CADJ),DrugProp,'LineWidth',3);hold on
    xticks([log(CADJ(1)),log(1e5),log(CADJ(8))]);
    xticklabels({'1e^1','1e^5','1e^9'});
    xlim([log(CADJ(1)) log(CADJ(end))])
    ylim([0 0.000012])
    plot([log(1e4) log(1e4)],[0 100],'k:','LineWidth',2); 
    plot([log(1e7) log(1e7)],[0 100],'k-','LineWidth',2);  
    title({'6-mo. Immunity &','Noncompliance to TR','(D)',''})
yyaxis right
    plot(log(CADJ),DrugInfection*100,'LineWidth',3);hold on
    ylabel({'Difference in','Cumulative Cases'},'Interpreter','tex','FontSize', 14); 
    ylim([-0.025 0])
subplot(248)
    plot(log(CADJ),CostProp,'LineWidth',3);hold on
    xticks([log(CADJ(1)),log(CADJ(5)),log(CADJ(8))]);
    xticklabels({'1e^1','1e^5','1e^9'});
    ylim([0 0.25])
    plot([log(1e4) log(1e4)],[0 100],'k:','LineWidth',2); 
    plot([log(1e7) log(1e7)],[0 100],'k-','LineWidth',2); 
    title({'(H)'})
    xlim([log(CADJ(1)) log(CADJ(end))])
    
    
    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    xlabel(han,'Workability Cost Parameter', 'FontSize', 16,'Interpreter','tex');
    
    saveas(gcf,'SA_Workability_Drug.png'); hold off
    
end


%Figure A25: Sensitivity of vaccine results to workability cost parameter
if true

    
    fig=figure
    
load Workspace_Cadj_PT.mat CADJ VaccineProp VaccineInfection CostProp
    
subplot(241)
yyaxis left
    plot(log(CADJ),VaccineProp*100,'LineWidth',3);hold on
    xticks([log(CADJ(1)),log(CADJ(5)),log(CADJ(9))]);
    xticklabels({'1e^1','1e^5','1e^9'});
    xlim([log(CADJ(1)) log(CADJ(end))])
    ylabel({'Variance of the Optimal ','Deviation From the Ad Hoc'},'Interpreter','tex','FontSize', 14);     
    plot([log(CADJ(4)) log(CADJ(4))],[0 100],'k:','LineWidth',2); 
    plot([log(CADJ(7)) log(CADJ(7))],[0 100],'k-','LineWidth',2); 
    ylim([0 100])
    title({'Perm. Immunity &','Compliance to TR','(A)'})
yyaxis right
    plot(log(CADJ),VaccineInfection*100,'LineWidth',3);hold on
    ylim([-0.4 0])
subplot(245)
    plot(log(CADJ),CostProp,'LineWidth',3);hold on
    plot([log(CADJ(4)) log(CADJ(4))],[0 6],'k:','LineWidth',2);
    plot([log(CADJ(7)) log(CADJ(7))],[0 6],'k-','LineWidth',2);
    xticks([log(CADJ(1)),log(CADJ(5)),log(CADJ(9))]);
    xticklabels({'1e^1','1e^5','1e^9'});
    ylabel({'$\frac{\textrm{Total Workability}}{\textrm{Total Vaccine}}$'},'Interpreter','latex', 'FontSize', 14);
    ylim([0 6])
    title({'(E)'})
    
    xlim([log(CADJ(1)) log(CADJ(end))])
    load Workspace_Cadj_PN.mat CADJ VaccineProp VaccineInfection CostProp
    
subplot(242)
yyaxis left
    plot(log(CADJ),VaccineProp*100,'LineWidth',3);hold on
    xticks([log(CADJ(1)),log(CADJ(5)),log(CADJ(9))]);
    xticklabels({'1e^1','1e^5','1e^9'});
    xlim([log(CADJ(1)) log(CADJ(end))])
    plot([log(CADJ(4)) log(CADJ(4))],[0 100],'k:','LineWidth',2); 
    plot([log(CADJ(7)) log(CADJ(7))],[0 100],'k-','LineWidth',2); 
    ylim([0 100])
    title({'Perm. Immunity &','Noncompliance to TR','(B)'})
yyaxis right
    plot(log(CADJ),VaccineInfection*100,'LineWidth',3);hold on
    ylim([-0.4 0])
subplot(246)
    plot(log(CADJ),CostProp,'LineWidth',3);hold on
    plot([log(CADJ(4)) log(CADJ(4))],[0 6],'k:','LineWidth',2);
    plot([log(CADJ(7)) log(CADJ(7))],[0 6],'k-','LineWidth',2);
    xticks([log(CADJ(1)),log(CADJ(5)),log(CADJ(9))]);
    xticklabels({'1e^1','1e^5','1e^9'});
    ylim([0 6])  
    title({'(F)'})
    xlim([log(CADJ(1)) log(CADJ(end))])
    
load Workspace_Cadj_ST.mat CADJ VaccineProp VaccineInfection CostProp

subplot(243)
yyaxis left
    plot(log(CADJ),VaccineProp*100,'LineWidth',3);hold on
    xticks([log(CADJ(1)),log(CADJ(5)),log(CADJ(9))]);
    xticklabels({'1e^1','1e^5','1e^9'});
    xlim([log(CADJ(1)) log(CADJ(end))])
    plot([log(CADJ(4)) log(CADJ(4))],[0 100],'k:','LineWidth',2); 
    plot([log(CADJ(7)) log(CADJ(7))],[0 100],'k-','LineWidth',2); 
    ylim([0 100])
    title({'6-mo. Immunity &','Compliance to TR','(C)'})
yyaxis right
    plot(log(CADJ),VaccineInfection*100,'LineWidth',3);hold on
    ylim([-0.4 0])
subplot(247)
    plot(log(CADJ),CostProp,'LineWidth',3);hold on
    plot([log(CADJ(4)) log(CADJ(4))],[0 6],'k:','LineWidth',2);
    plot([log(CADJ(7)) log(CADJ(7))],[0 6],'k-','LineWidth',2);
    plot([log(CADJ(4)) log(CADJ(4))],[0 100],'k:','LineWidth',2); 
    plot([log(CADJ(7)) log(CADJ(7))],[0 100],'k-','LineWidth',2); 
    xticklabels({'1e^1','1e^5','1e^9'});
    ylim([0 6])
    title({'(G)'})
    xlim([log(CADJ(1)) log(CADJ(end))])
    
 load Workspace_Cadj_SN.mat CADJ VaccineProp VaccineInfection CostProp

subplot(244)
yyaxis left
    plot(log(CADJ),VaccineProp*100,'LineWidth',3);hold on
    xticks([log(CADJ(1)),log(CADJ(5)),log(CADJ(9))]);
    xticklabels({'1e^1','1e^5','1e^9'});
    xlim([log(CADJ(1)) log(CADJ(end))])
    plot([log(CADJ(4)) log(CADJ(4))],[0 100],'k:','LineWidth',2); 
    plot([log(CADJ(7)) log(CADJ(7))],[0 100],'k-','LineWidth',2); 
    ylim([0 100])
    title({'6-mo. Immunity &','Noncompliance to TR','(D)'})
yyaxis right
    plot(log(CADJ),VaccineInfection*100,'LineWidth',3);hold on
    ylabel({'Difference in','Cumulative Cases'},'Interpreter','tex','FontSize', 14);
    ylim([-0.4 0])
subplot(248)
    plot(log(CADJ),CostProp,'LineWidth',3);hold on
    plot([log(CADJ(4)) log(CADJ(4))],[0 6],'k:','LineWidth',2);
    plot([log(CADJ(7)) log(CADJ(7))],[0 6],'k-','LineWidth',2);
    xticks([log(CADJ(1)),log(CADJ(5)),log(CADJ(9))]);
    xticklabels({'1e^1','1e^5','1e^9'});
    ylim([0 6])
    title({'(H)'})
    xlim([log(CADJ(1)) log(CADJ(end))])
    
    
    
    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    xlabel(han,'Workability Cost Parameter', 'FontSize', 16,'Interpreter','tex');
    
    saveas(gcf,'SA_Workability.png'); hold off
    
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%   Figures for effectiveness of controls  %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Figure A14: Sensitivity of results to drug effectiveness 
if true

    
    fig=figure
    
load Workspace_qD_PT.mat QD DrugProp DrugInfection CostProp DrugProp2
    
subplot(241)
yyaxis left
    plot(QD,DrugProp,'LineWidth',3);hold on
    xticks([QD(1),QD(3),QD(5)]);
    xlim([QD(1) QD(end)])
    ylabel({'Variance of the Optimal ','Deviation From the Ad Hoc'},'Interpreter','tex','FontSize', 14);     
    plot([QD(2) QD(2)],[0 40],'k:','LineWidth',2); 
    ylim([0 1.5e-5])
    xlim([0.55 0.95])
    title({'Perm. Immunity &','Compliance to TR','(A)',''})
yyaxis right
    plot(QD,DrugInfection*100,'LineWidth',3);hold on 
    ylim([-5 0])
subplot(245)
    plot(QD,CostProp,'LineWidth',3);hold on
    plot([QD(2) QD(2)],[0 4],'k:','LineWidth',2);
    xticks([QD(1),QD(3),QD(5)]);
    ylabel({'$\frac{\textrm{Total Workability}}{\textrm{Total Vaccine}}$'},'Interpreter','latex', 'FontSize', 14);
    ylim([0 0.03])
    xlim([0.55 0.95])
    title({'(E)'})

    load Workspace_qD_PN.mat CADJ DrugProp DrugInfection CostProp
    
subplot(242)
yyaxis left
    plot(QD,DrugProp,'LineWidth',3);hold on
    xticks([QD(1),QD(3),QD(5)]);
    xticklabels({'0.55','0.75','0.95'});
    xlim([QD(1) QD(end)])
    plot([QD(2) QD(2)],[0 40],'k:','LineWidth',2); 
    ylim([0 1.5e-5])
    title({'Perm. Immunity &','Noncompliance to TR','(B)',''})
yyaxis right
    plot(QD,DrugInfection*100,'LineWidth',3);hold on
    ylim([-0.05 0])
subplot(246)
    plot(QD,CostProp,'LineWidth',3);hold on
    plot([QD(2) QD(2)],[0 4],'k:','LineWidth',2);
    xticks([QD(1),QD(3),QD(5)]);
    xticklabels({'0.55','0.75','0.95'});
    ylim([0 0.03])  
    xlim([0.55 0.95])
    title({'(F)'})
        
load Workspace_qD_ST.mat CADJ DrugProp DrugInfection CostProp

subplot(243)
yyaxis left
    plot(QD,DrugProp,'LineWidth',3);hold on
    xticks([QD(1),QD(3),QD(5)]);
    xticklabels({'0.55','0.75','0.95'});
    xlim([QD(1) QD(end)])
    plot([QD(2) QD(2)],[0 40],'k:','LineWidth',2); 
    ylim([0 1.5e-5])
    title({'6-mo. Immunity &','Compliance to TR','(C)',''})   
yyaxis right
    plot(QD,DrugInfection*100,'LineWidth',3);hold on 
    ylim([-5 0])
subplot(247)
    plot(QD,CostProp,'LineWidth',3);hold on
    plot([QD(2) QD(2)],[0 4],'k:','LineWidth',2);
    xticks([QD(1),QD(3),QD(5)]);
    xticklabels({'0.55','0.75','0.95'});
    ylim([0 0.03])
    xlim([0.55 0.95])
    title({'(G)'})    
 load Workspace_qD_SN.mat CADJ DrugProp DrugInfection CostProp

subplot(244)
yyaxis left
    plot(QD,DrugProp,'LineWidth',3);hold on
    xticks([QD(1),QD(3),QD(5)]);
    xticklabels({'0.55','0.75','0.95'});
    xlim([QD(1) QD(end)])
    plot([QD(2) QD(2)],[0 40],'k:','LineWidth',2); 
    ylim([0 1.5e-5])
    title({'6-mo. Immunity &','Noncompliance to TR','(D)',''})
yyaxis right
    plot(QD,DrugInfection*100,'LineWidth',3);hold on
    ylabel({'Difference in','Cumulative Cases'},'Interpreter','tex','FontSize', 14);
    ylim([-0.05 0])
subplot(248)
    plot(QD,CostProp,'LineWidth',3);hold on
    plot([QD(2) QD(2)],[0 4],'k:','LineWidth',2);
    xticks([QD(1),QD(3),QD(5)]);
    xticklabels({'0.55','0.75','0.95'});
    ylim([0 0.03])
    xlim([0.55 0.95])
    title({'(H)'})
    
    
    
    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    xlabel(han,'Effectiveness of Drug', 'FontSize', 16,'Interpreter','tex');
    
    saveas(gcf,'SA_Effectivenss_Drug.png'); hold off
    
end


%Figure A26: Sensitivity of results to vaccine effectiveness 
if true

    
    fig=figure
    
load Workspace_qV_PT.mat QV VaccineProp VaccineInfection CostProp
    
subplot(241)
yyaxis left
    plot(QV,VaccineProp*100,'LineWidth',3);hold on
    xticks([QV(1),QV(3),QV(5)]);
    xlim([QV(1) QV(end)])
    ylabel({'Variance of the Optimal ','Deviation From the Ad Hoc'},'Interpreter','tex','FontSize', 14);     
    plot([QV(2) QV(2)],[0 40],'k:','LineWidth',2); 
    ylim([0 40])
    title({'Perm. Immunity &','Compliance to TR','(A)'})
yyaxis right
    plot(QV,VaccineInfection*100,'LineWidth',3);hold on 
    ylim([-0.4 0])
subplot(245)
    plot(QV,CostProp,'LineWidth',3);hold on
    plot([QV(2) QV(2)],[0 4],'k:','LineWidth',2);
    xticks([QV(1),QV(3),QV(5)]);
    ylabel({'$\frac{\textrm{Total Workability}}{\textrm{Total Vaccine}}$'},'Interpreter','latex', 'FontSize', 14);
    ylim([0 4])
    xlim([0.55 0.95])
    title({'(E)'})

    load Workspace_qV_PN.mat CADJ VaccineProp VaccineInfection CostProp
    
subplot(242)
yyaxis left
    plot(QV,VaccineProp*100,'LineWidth',3);hold on
    xticks([QV(1),QV(3),QV(5)]);
    xticklabels({'0.55','0.75','0.95'});
    xlim([QV(1) QV(end)])
    plot([QV(2) QV(2)],[0 40],'k:','LineWidth',2); 
    ylim([0 40])
    title({'Perm. Immunity &','Noncompliance to TR','(B)'})
yyaxis right
    plot(QV,VaccineInfection*100,'LineWidth',3);hold on
    ylim([-0.4 0])
subplot(246)
    plot(QV,CostProp,'LineWidth',3);hold on
    plot([QV(2) QV(2)],[0 4],'k:','LineWidth',2);
    xticks([QV(1),QV(3),QV(5)]);
    xticklabels({'0.55','0.75','0.95'});
    ylim([0 4])  
    xlim([0.55 0.95])
    title({'(F)'})
        
load Workspace_qV_ST.mat CADJ VaccineProp VaccineInfection CostProp

subplot(243)
yyaxis left
    plot(QV,VaccineProp*100,'LineWidth',3);hold on
    xticks([QV(1),QV(3),QV(5)]);
    xticklabels({'0.55','0.75','0.95'});
    xlim([QV(1) QV(end)])
    plot([QV(2) QV(2)],[0 40],'k:','LineWidth',2); 
    ylim([0 40])
    title({'6-mo. Immunity &','Compliance to TR','(C)'})   
yyaxis right
    plot(QV,VaccineInfection*100,'LineWidth',3);hold on 
    ylim([-0.4 0])
subplot(247)
    plot(QV,CostProp,'LineWidth',3);hold on
    plot([QV(2) QV(2)],[0 4],'k:','LineWidth',2);
    xticks([QV(1),QV(3),QV(5)]);
    xticklabels({'0.55','0.75','0.95'});
    ylim([0 4])
    xlim([0.55 0.95])
    title({'(G)'})    
 load Workspace_qV_SN.mat CADJ VaccineProp VaccineInfection CostProp

subplot(244)
yyaxis left
    plot(QV,VaccineProp*100,'LineWidth',3);hold on
    xticks([QV(1),QV(3),QV(5)]);
    xticklabels({'0.55','0.75','0.95'});
    xlim([QV(1) QV(end)])
    plot([QV(2) QV(2)],[0 40],'k:','LineWidth',2); 
    ylim([0 40])
    title({'6-mo. Immunity &','Noncompliance to TR','(D)'})
yyaxis right
    plot(QV,VaccineInfection*100,'LineWidth',3);hold on
    ylabel({'Difference in','Cumulative Cases'},'Interpreter','tex','FontSize', 14);
    ylim([-0.4 0])
subplot(248)
    plot(QV,CostProp,'LineWidth',3);hold on
    plot([QV(2) QV(2)],[0 4],'k:','LineWidth',2);
    xticks([QV(1),QV(3),QV(5)]);
    xticklabels({'0.55','0.75','0.95'});
    ylim([0 4])
    xlim([0.55 0.95])
    title({'(H)'})
    
    
    
    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    xlabel(han,'Effectiveness of Vaccine', 'FontSize', 16,'Interpreter','tex');
    
    saveas(gcf,'SA_Effectivenss_Vaccine.png'); hold off
    
end



