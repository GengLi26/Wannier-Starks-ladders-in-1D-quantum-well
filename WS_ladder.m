% WS ladder _ 10 period
% 1st column E field (kV/cm)
% 2nd column \Delta E
WS=[5, 0.006;
10, 0.0075;
15, 0.0106;
20, 0.014;
25, 0.0175;
30, 0.021;
35, 0.0245;
40, 0.0280;
45, 0.0315;
50, 0.0350];

%%
scatter(WS(:,1),WS(:,2),100,'k')
hold on 
%% deltaE = eFd;
F_field_fit=[0:0.0001:55];
Period_fit=7e-9;
deltaE=F_field_fit.*1e5.*Period_fit;
plot(F_field_fit,deltaE,'color','k','linewidth',1.5)
hold on 
xlim([0,55])
ylim([0,0.04])
set(gca,'fontsize',28)
xlabel(['F (kV/cm)'],'FontSize',28)
ylabel(['\DeltaE (eV)'],'FontSize',28)
set(gcf,'Position',[500 300 800 600])
legend({' \DeltaE',' \DeltaE=eFd'},'FontSize')
legend('boxoff') 
box on 
set(gcf,'PaperOrientation','landscape')
print(gcf, 'WS_eFD.pdf', '-dpdf','-r0','-bestfit')

