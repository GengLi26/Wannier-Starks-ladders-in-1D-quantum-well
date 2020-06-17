clear all 
clc
% close all

Left_contact=40e-9; % Left contact 30 nm, right 30 nm. 
Well=3.5e-9; % barrier width 3.5 nm. Well width 3.5 nm. 
Well_number=40;  % to get WS ladder ,set well nubmer to be 10 
Stop=40; % the number of eigen states/energies to be displayed
E_contact=-0e5; % Field strength V/m   converted from 30 kV/cm 
N=7000;  %% Number of mesh points. Run a convergence test before trusting the results. 
%% Solve 1D schodinger equation numerically.( Finite Difference Method)
e = 1.602177e-19; % electron charge (C)
hbar = 1.054572e-34; % reduced Planck constant ( J . s )
%m_e = 9.109384e-31; % electron mass (kg)
m_e = 0.067*9.109384e-31; % electron mass (kg) GaAs 

L=Left_contact*2+Well*Well_number*2;  %% width of the discretized region. 
dx=L/(N-1);
mesh_x=[0:dx:L];
mesh_x(1)=[]; % Remove the first mesh points Psi_0
mesh_x(length(mesh_x))=[];% Remove the last mesh points Psi_(N+1)

% Definition
t0=hbar^2./(2*m_e*dx^2);


% Multiple quantum well peridicity 20*7 nm = 140 nm. 

% find the start point of the first well. 
[test1,test2]=min(abs((mesh_x-Left_contact)));
% find the end point of the first well.
[test3,test4]=min(abs((mesh_x-Left_contact-Well)));

% \delta E= (0.87 +_ 0.04)* x eV 
%Ref{HH Liu Quantum well infrared photodetector book}
Height=0.87*0.3*e; % Barrer height of quantum well Al0.3Ga0.7As/GaAs 
V=ones(1,N-2).*Height; % define intial potential
period=(test4-test2)*2;% How many mesh points in a well? 

% define potential 
for V_count=1:1:Well_number
    
    V((test2+V_count*period):(test4+V_count*period))=0;
end

% Apply Electrical field to the potential 
Efield=e.*E_contact.*mesh_x;
V=V+Efield;

% % plot the potential if needed
% plot(mesh_x*1e9,V./e)


%% Generate Hamiltonian for the system
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

Hamil=Hamil+diag(V); % Add potential to the Hamiltonian.
[Vector,E_eig]=eig(Hamil); % Solve eigen value and eigen functions.
Band=real(diag(E_eig))./e; % Convert eigen energies to eV. 

% %%Plot potential and wavefunction 
figure(1)
% potential 
plot(mesh_x*1e9,V./e,'color','k','linewidth',1)
hold on 

% wavefunctions
for plot_count=1:1:Stop
    plot(mesh_x*1e9,Vector(:,plot_count).^2*5+Band(plot_count)...
        ,'color','k','linewidth',2)
hold on 
end

% edit figure
set(gca,'fontsize',28)
xlabel(['z (nm)'],'FontSize',28)
ylabel(['E (eV)'],'FontSize',28)
set(gcf,'Position',[500 300 800 600])

% plot eigen energies
figure(2)
scatter([0:1:Stop-1],Band(1:Stop))
hold on 
%% Edit figure properties
% % xlim([37,95])
% % ylim([-0.5,0.1])
% % 
% % annotation('textbox',...
%     [0.74 0.80 0.25 0.063],...
%     'String',{'50 kV/cm'},...
%     'LineWidth',2,...
%     'LineStyle','none',...
%     'FontSize',26);
% 
% set(gcf,'PaperOrientation','landscape')
% print(gcf, '3.5 nm well _ E=50kV_cm.pdf', '-dpdf','-r0','-bestfit')
% figure(2)
% Band_test1=Band;
% Band_test1(1)=[];
% Band_test1(length(Band))= 0;
% Band_test2=Band_test1-Band;
% plot(Band_test2(1:20))

