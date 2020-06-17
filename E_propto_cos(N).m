
%% Fit e field at zero to cos (N)
% E=0;
Band_remove_trivial=Band(1:Stop); % get eigen energies;
Band_zero=Band_remove_trivial-min(Band_remove_trivial);

scatter([0:1:Stop-1],Band_zero,100,'k')
hold on

x_cos=[-2.1:0.1:39];
E_cos_fit=1-cos(x_cos./(2*pi)/2);
plot(x_cos+2.1,E_cos_fit./max(E_cos_fit)*max(Band_zero),'color','k','linewidth',2)

box on
xlim([0,40])
ylim([0,0.06])
set(gca,'fontsize',28)
xlabel(['N'],'FontSize',28)
ylabel(['E (eV)'],'FontSize',28)
set(gcf,'Position',[500 300 800 600])

legend({' Eigen energies',' E $\propto$ cos(N)'},'FontSize',26,'Interpreter','latex')
legend('boxoff') 
   
set(gcf,'PaperOrientation','landscape')
print(gcf, 'WS_zero_field_E_vs_cosN.pdf', '-dpdf','-r0','-bestfit')
