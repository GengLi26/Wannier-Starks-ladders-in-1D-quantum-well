clear all 
clc

%% Solve 1D schodinger equation numerically. 
e = 1.602177e-19; % electron charge (C)
hbar = 1.054572e-34; % reduced Planck constant ( J . s )
m_e = 9.109384e-31; % electron mass (kg)

L=200e-9;  %% width of the discretized region. 
N=6000; 
dx=L/(N-1);
mesh_x=[0:dx:L];
mesh_x(1)=[];
mesh_x(length(mesh_x))=[];

t0=hbar^2./(2*m_e*dx^2);

E_harm=0.1; %% unit eV
omega=E_harm*e./hbar;
%V=zeros(1,N-2);

% Harmonic oscillator 
E_harm=0.1; %% unit eV
omega=E_harm*e./hbar;
V=1/2*m_e*omega^2.*(mesh_x-L/2).^2;

%plot(mesh_x*1e9,V./e)

Hamil=zeros(N-2,N-2);
Hamil(1,1)=2*t0;
Hamil(1,2)=-1*t0;
Hamil(N-2,N-3)=-1*t0;
Hamil(N-2,N-2)=2*t0;

for count=2:1:N-3
    
    Hamil(count,count-1)=-1*t0;
    Hamil(count,count)=2*t0;
    Hamil(count,count+1)=-1*t0;

end

Hamil=Hamil+diag(V);
[Vector,E_eig]=eig(Hamil);
Band=real(diag(E_eig))./e;

scatter([0:1:19],Band(1:20),100,'k')
hold on 
x=[0:0.0001:20];
% y=0.1*x-0.05;
y=hbar*omega*(x+1/2)./e;
plot(x,y,'color','k','linewidth',1.5)
hold on 
ylim([0,2])
set(gca,'fontsize',28)
xlabel(['N'],'FontSize',28)
ylabel(['E (eV)'],'FontSize',28)
set(gcf,'Position',[500 300 800 600])
legend({' Eigen energy',' E=$\hbar\omega(n+1/2)$'},'FontSize',26,'Interpreter','latex')
legend('boxoff') 
box on 
set(gcf,'PaperOrientation','landscape')
print(gcf, 'Oscillator_xfitting.pdf', '-dpdf','-r0','-bestfit')


plot(mesh_x*1e9,V./e,'color','k','linewidth',1)
hold on 

for plot_count=1:1:10
plot(mesh_x*1e9,Vector(:,plot_count).^2*3+Band(plot_count),'color','k','linewidth',2)
hold on 

end
xlim([94,106])
ylim([0,1.2])

set(gca,'fontsize',28)
xlabel(['z (nm)'],'FontSize',28)
ylabel(['E (eV)'],'FontSize',28)
set(gcf,'Position',[500 300 800 600])


set(gcf,'PaperOrientation','landscape')
print(gcf, 'Oscillator_wave.pdf', '-dpdf','-r0','-bestfit')

