clear all 
clc

%% Solve 1D schodinger equation numerically. 
e = 1.602177e-19; % electron charge (C)
hbar = 1.054572e-34; % reduced Planck constant ( J . s )
%m_e = 9.109384e-31; % electron mass (kg)

m_e = 0.067*9.109384e-31; % electron mass (kg) GaAs 

L=20e-9;  %% width of the discretized region. 
N=3000; 
dx=L/(N-1);
mesh_x=[0:dx:L];
mesh_x(1)=[];
mesh_x(length(mesh_x))=[];
t0=hbar^2./(2*m_e*dx^2);

%% Quantum well width 6 nm, depth 0.4 eV.
Width=3.5e-9; % witdh of quantum well 
Well1=L/2-Width/2;
Well2=L/2+Width/2;
[test1,test2]=min(abs((mesh_x-Well1)));
[test3,test4]=min(abs((mesh_x-Well2)));
V=zeros(1,N-2);
V(1:test2)=0.87*0.3*e;
V(test4:length(V))=0.87*0.3*e;
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
% scatter([1:1:5],Band(1:5))
% hold on 

plot(mesh_x*1e9,V./e,'color','k','linewidth',1)
hold on 
plot(mesh_x*1e9,Vector(:,1).^2*40+Band(1),'color','k','linewidth',2)
hold on 
 % 3.5 nm well 0.1233 nm
set(gca,'fontsize',28)
xlabel(['z (nm)'],'FontSize',28)
ylabel(['E (eV)'],'FontSize',28)
set(gcf,'Position',[500 300 800 600])

xlim([0,20])
ylim([-0.05,0.3])

set(gcf,'PaperOrientation','landscape')
print(gcf, '3point5nm quantum well.pdf', '-dpdf','-r0','-bestfit')

