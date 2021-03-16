clear all; close all

%To find where capacity stops being constraining
if true
  SCARCITY=[0.05 0.1 0.15]
  
 
load WorkspacePT_Death.mat
    scarcePT_V10_Death=max(find(round(uV1s10+uV2s10,4)==round(MaxTreat(2)*SCARCITY(2),4)));

load WorkspacePN_Death.mat
    scarcePN_V10_Death=max(find(round(uV1s10+uV2s10,4)==round(MaxTreat(2)*SCARCITY(2),4)));
    
load WorkspaceST_Death.mat
    scarceST_V10_Death=max(find(round(uV1s10+uV2s10,4)==round(MaxTreat(2)*SCARCITY(2),4)));
    
load WorkspaceSN_Death.mat
    scarceSN_V10_Death=max(find(round(uV1s10+uV2s10,4)==round(MaxTreat(2)*SCARCITY(2),4)));
 
    
 load WorkspacePT_Contact.mat
    scarcePT_V10_Contact=max(find(round(uV1s10+uV2s10,4)==round(MaxTreat(2)*SCARCITY(2),4)));

load WorkspacePN_Contact.mat
    scarcePN_V10_Contact=max(find(round(uV1s10+uV2s10,4)==round(MaxTreat(2)*SCARCITY(2),4)));
    
load WorkspaceST_Contact.mat
    scarceST_V10_Contact=max(find(round(uV1s10+uV2s10,4)==round(MaxTreat(2)*SCARCITY(2),4)));
    
load WorkspaceSN_Contact.mat
    scarceSN_V10_Contact=max(find(round(uV1s10+uV2s10,4)==round(MaxTreat(2)*SCARCITY(2),4)));
    
end


%Different Contact Rates across jurisdictions with Permanent Immunity
if false
  
bl=[0, 0.4470, 0.7410];
re=[0.8500, 0.3250, 0.0980];

    
    
fig= figure

ts=ts10;

load WorkspacePT_Contact.mat
    subplot1=subplot(2,2,1)
    p1=plot(ts,uV1s10,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,uV2s10,'LineWidth',3,'Color',re);hold on
    p1=plot(ts,uV1s10_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on
    p2=plot(ts,uV2s10_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    xlim([0 ts(scarcePT_V10_Contact)])
    ylim([0 SCARCITY(2)*MaxTreat(2)])
    title({'Perm. Immunity &','Compliance to TR','(A)'})
    ylabel({'Prop. of','Treated'}, 'FontSize', 16); 
 
    
    subplot(2,2,3)
    p1=plot(ts,I1s10,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,I2s10,'LineWidth',3,'Color',re);hold on
    p3=plot(ts,I1s10_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on 
    p4=plot(ts,I2s10_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    ylabel({'Prop. of','Infected'}, 'FontSize', 16); 
    xlim([0 ts(scarcePT_V10_Contact)])
    ylim([0 0.035])
    title({'(C)'})
    legend1=legend([p1 p3 p2 p4],{'State 1: Optimal~~~','State 1: Ad Hoc~~~','State 2: Optimal~~~','State 2: Ad Hoc'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
    'Position',[0.106626425988088 0.493992323950414 0.811608069283622 0.0352380951245623],...
    'Orientation','horizontal',...
    'Interpreter','latex');

load WorkspacePN_Contact.mat
    subplot(2,2,2)
    p1=plot(ts,uV1s10,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,uV2s10,'LineWidth',3,'Color',re);hold on
    p1=plot(ts,uV1s10_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on
    p2=plot(ts,uV2s10_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    xlim([0 ts(scarcePN_V10_Contact)])
    ylim([0 SCARCITY(2)*MaxTreat(2)])
    title({'Perm. Immunity &','Noncompliance to TR','(B)'})
    ylabel({'Prop. of','Treated'}, 'FontSize', 16); 
 
    
    subplot(2,2,4)
    p1=plot(ts,I1s10,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,I2s10,'LineWidth',3,'Color',re);hold on
    p3=plot(ts,I1s10_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on 
    p4=plot(ts,I2s10_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    ylabel({'Prop. of','Infected'}, 'FontSize', 16); 
    xlim([0 ts(scarcePN_V10_Contact)])
    ylim([0 0.1])
    title({'(D)'})
    
    

han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
%ylabel(han,'Quantity of Treatment', 'FontSize', 16);

xlabel(han,'Time (months)', 'FontSize', 16);

%sgtitle({'Cumulative Difference in the Number of Infections','Compared to the No-Treatment Scenario for Different','Drug Allocation Rules when Capacity is 10%'},'Color','Black', 'FontSize', 18);
saveas(gcf,'PermImm_Contact_10percent.png'); hold off

end


%Different Contact Rates across jurisdictions with Temporary Immunity
if false
  
bl=[0, 0.4470, 0.7410];
re=[0.8500, 0.3250, 0.0980];

    
    
fig= figure

ts=ts10;

load WorkspaceST_Contact.mat
    subplot1=subplot(2,2,1)
    p1=plot(ts,uV1s10,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,uV2s10,'LineWidth',3,'Color',re);hold on
    p1=plot(ts,uV1s10_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on
    p2=plot(ts,uV2s10_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    xlim([0 ts(scarceST_V10_Contact)])
    ylim([0 SCARCITY(2)*MaxTreat(2)])
    title({'6-mo. Immunity &','Compliance to TR','(A)'})
    ylabel({'Prop. of','Treated'}, 'FontSize', 16); 
 
    
    subplot(2,2,3)
    p1=plot(ts,I1s10,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,I2s10,'LineWidth',3,'Color',re);hold on
    p3=plot(ts,I1s10_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on 
    p4=plot(ts,I2s10_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    ylabel({'Prop. of','Infected'}, 'FontSize', 16); 
    xlim([0 ts(scarceST_V10_Contact)])
    ylim([0 0.035])
    title({'(C)'})
    legend1=legend([p1 p3 p2 p4],{'State 1: Optimal~~~','State 1: Ad Hoc~~~','State 2: Optimal~~~','State 2: Ad Hoc'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
    'Position',[0.106626425988088 0.493992323950414 0.811608069283622 0.0352380951245623],...
    'Orientation','horizontal',...
    'Interpreter','latex');

load WorkspaceSN_Contact.mat
    subplot(2,2,2)
    p1=plot(ts,uV1s10,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,uV2s10,'LineWidth',3,'Color',re);hold on
    p1=plot(ts,uV1s10_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on
    p2=plot(ts,uV2s10_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    xlim([0 ts(scarceSN_V10_Contact)])
    ylim([0 SCARCITY(2)*MaxTreat(2)])
    title({'6-mo. Immunity &','Noncompliance to TR','(B)'})
    ylabel({'Prop. of','Treated'}, 'FontSize', 16); 
 
    
    subplot(2,2,4)
    p1=plot(ts,I1s10,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,I2s10,'LineWidth',3,'Color',re);hold on
    p3=plot(ts,I1s10_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on 
    p4=plot(ts,I2s10_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    ylabel({'Prop. of','Infected'}, 'FontSize', 16); 
    xlim([0 ts(scarceSN_V10_Contact)])
    ylim([0 0.1])
    title({'(D)'})
    
    

han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
%ylabel(han,'Quantity of Treatment', 'FontSize', 16);

xlabel(han,'Time (months)', 'FontSize', 16);

%sgtitle({'Cumulative Difference in the Number of Infections','Compared to the No-Treatment Scenario for Different','Drug Allocation Rules when Capacity is 10%'},'Color','Black', 'FontSize', 18);
saveas(gcf,'TempImm_Contact_10percent.png'); hold off

end



%Different Death Rates across jurisdictions with Permanent Immunity
if true
  
bl=[0, 0.4470, 0.7410];
re=[0.8500, 0.3250, 0.0980];

    
    
fig= figure

ts=ts10;

load WorkspacePT_Death.mat
    subplot1=subplot(2,2,1)
    p1=plot(ts,uV1s10,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,uV2s10,'LineWidth',3,'Color',re);hold on
    p1=plot(ts,uV1s10_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on
    p2=plot(ts,uV2s10_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    xlim([0 ts(scarcePT_V10_Death)])
    ylim([0 SCARCITY(2)*MaxTreat(2)])
    title({'Perm. Immunity &','Compliance to TR','(A)'})
    ylabel({'Prop. of','Treated'}, 'FontSize', 16); 
 
    
    subplot(2,2,3)
    p1=plot(ts,I1s10,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,I2s10,'LineWidth',3,'Color',re);hold on
    p3=plot(ts,I1s10_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on 
    p4=plot(ts,I2s10_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    ylabel({'Prop. of','Infected'}, 'FontSize', 16); 
    xlim([0 ts(scarcePT_V10_Death)])
    ylim([0 0.035])
    title({'(C)'})
    legend1=legend([p1 p3 p2 p4],{'State 1: Optimal~~~','State 1: Ad Hoc~~~','State 2: Optimal~~~','State 2: Ad Hoc'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
    'Position',[0.106626425988088 0.493992323950414 0.811608069283622 0.0352380951245623],...
    'Orientation','horizontal',...
    'Interpreter','latex');

load WorkspacePN_Death.mat
    subplot(2,2,2)
    p1=plot(ts,uV1s10,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,uV2s10,'LineWidth',3,'Color',re);hold on
    p1=plot(ts,uV1s10_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on
    p2=plot(ts,uV2s10_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    xlim([0 ts(scarcePN_V10_Death)])
    ylim([0 SCARCITY(2)*MaxTreat(2)])
    title({'Perm. Immunity &','Noncompliance to TR','(B)'})
    ylabel({'Prop. of','Treated'}, 'FontSize', 16); 
 
    
    subplot(2,2,4)
    p1=plot(ts,I1s10,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,I2s10,'LineWidth',3,'Color',re);hold on
    p3=plot(ts,I1s10_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on 
    p4=plot(ts,I2s10_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    ylabel({'Prop. of','Infected'}, 'FontSize', 16); 
    xlim([0 ts(scarcePN_V10_Death)])
    ylim([0 0.1])
    title({'(D)'})
    
    

han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
%ylabel(han,'Quantity of Treatment', 'FontSize', 16);

xlabel(han,'Time (months)', 'FontSize', 16);

%sgtitle({'Cumulative Difference in the Number of Infections','Compared to the No-Treatment Scenario for Different','Drug Allocation Rules when Capacity is 10%'},'Color','Black', 'FontSize', 18);
saveas(gcf,'PermImm_Death_10percent.png'); hold off

end


%Different Death Rates across jurisdictions with Temporary Immunity
if true
  
bl=[0, 0.4470, 0.7410];
re=[0.8500, 0.3250, 0.0980];

    
    
fig= figure

ts=ts10;

load WorkspaceST_Death.mat
    subplot1=subplot(2,2,1)
    p1=plot(ts,uV1s10,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,uV2s10,'LineWidth',3,'Color',re);hold on
    p1=plot(ts,uV1s10_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on
    p2=plot(ts,uV2s10_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    xlim([0 ts(scarceST_V10_Death)])
    ylim([0 SCARCITY(2)*MaxTreat(2)])
    title({'6-mo. Immunity &','Compliance to TR','(A)'})
    ylabel({'Prop. of','Treated'}, 'FontSize', 16); 
 
    
    subplot(2,2,3)
    p1=plot(ts,I1s10,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,I2s10,'LineWidth',3,'Color',re);hold on
    p3=plot(ts,I1s10_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on 
    p4=plot(ts,I2s10_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    ylabel({'Prop. of','Infected'}, 'FontSize', 16); 
    xlim([0 ts(scarceST_V10_Death)])
    ylim([0 0.035])
    title({'(C)'})
    legend1=legend([p1 p3 p2 p4],{'State 1: Optimal~~~','State 1: Ad Hoc~~~','State 2: Optimal~~~','State 2: Ad Hoc'},'Interpreter','latex','Orientation','horizontal','Location','northeast');
set(legend1,...
    'Position',[0.106626425988088 0.493992323950414 0.811608069283622 0.0352380951245623],...
    'Orientation','horizontal',...
    'Interpreter','latex');

load WorkspaceSN_Death.mat
    subplot(2,2,2)
    p1=plot(ts,uV1s10,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,uV2s10,'LineWidth',3,'Color',re);hold on
    p1=plot(ts,uV1s10_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on
    p2=plot(ts,uV2s10_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    xlim([0 ts(scarceSN_V10_Death)])
    ylim([0 SCARCITY(2)*MaxTreat(2)])
    title({'6-mo. Immunity &','Noncompliance to TR','(B)'})
    ylabel({'Prop. of','Treated'}, 'FontSize', 16); 
 
    
    subplot(2,2,4)
    p1=plot(ts,I1s10,'LineWidth',3,'Color',bl);hold on
    p2=plot(ts,I2s10,'LineWidth',3,'Color',re);hold on
    p3=plot(ts,I1s10_ah,'LineWidth',3,'Color',bl,'LineStyle',':');hold on 
    p4=plot(ts,I2s10_ah,'LineWidth',3,'Color',re,'LineStyle',':');hold on
    ylabel({'Prop. of','Infected'}, 'FontSize', 16); 
    xlim([0 ts(scarceSN_V10_Death)])
    ylim([0 0.1])
    title({'(D)'})
    
    

han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
%ylabel(han,'Quantity of Treatment', 'FontSize', 16);

xlabel(han,'Time (months)', 'FontSize', 16);

%sgtitle({'Cumulative Difference in the Number of Infections','Compared to the No-Treatment Scenario for Different','Drug Allocation Rules when Capacity is 10%'},'Color','Black', 'FontSize', 18);
saveas(gcf,'TempImm_Death_10percent.png'); hold off

end
