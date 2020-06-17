clear all 
clc

%% Solve 1D schodinger equation numerically. 
e = 1.602177e-19; % electron charge (C)
hbar = 1.054572e-34; % reduced Planck constant ( J . s )
m_e = 9.109384e-31; % electron mass (kg)

%m_e = 0.067*9.109384e-31; % electron mass (kg) GaAs 

L=3000e-9;  %% width of the discretized region. 
N=6000; 
dx=L/(N-1);
mesh_x=[0:dx:L];
mesh_x(1)=[];
mesh_x(length(mesh_x))=[];
t0=hbar^2./(2*m_e*dx^2);


V=zeros(1,N-2);

% plot(mesh_x*1e9,V./e)

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
% 
scatter([1:1:20],Band(1:20),100,'k')
hold on 
x=[0:0.0001:20];
y=0.0000000415*x.^2;
plot(x,y,'color','k','linewidth',1.5)
hold on 
set(gca,'fontsize',28)
xlabel(['N'],'FontSize',28)
ylabel(['E (eV)'],'FontSize',28)
set(gcf,'Position',[500 300 800 600])
legend({' Eigen energy',' y=x^2 fitting'},'FontSize',26)
legend('boxoff') 
box on 
set(gcf,'PaperOrientation','landscape')
print(gcf, '2000nm_xsquarefitting.pdf', '-dpdf','-r0','-bestfit')


%%
plot(mesh_x*1e9,V./e,'color','k','linewidth',1)
hold on 


for plot_count=1:1:10
plot(mesh_x*1e9,Vector(:,plot_count).^2*0.0002+Band(plot_count),'color','k','linewidth',2)
hold on 

end

 % 3.5 nm well 0.1233 nm
set(gca,'fontsize',28)
xlabel(['z (nm)'],'FontSize',28)
ylabel(['E (eV)'],'FontSize',28)
set(gcf,'Position',[500 300 800 600])


set(gcf,'PaperOrientation','landscape')
print(gcf, '2000nm_freespace.pdf', '-dpdf','-r0','-bestfit')

