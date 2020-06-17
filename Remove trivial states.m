%% Remove trivial states 
% Give bands
% V_ref=Height+Efield;

% Etract the first 'Stop'-th eigen energies
Band_remove_trivial=Band(1:Stop);

%% Remove trivial states if E field is high 
Ref_number=1; % use the first wavefunction as reference. 
Point_ref_1=L-Left_contact+Well; % unit m the start position of right contact
Cut_number_1=round(Point_ref_1./L*length(mesh_x)); % Matrix lable of 'Point_ref_1;
test=Vector(Cut_number_1:length(mesh_x),Ref_number).^2; % Summation of wavefunction in right contact. 

test_ref=sum(test)*3; % Adjustable to properly remove trivial states

count_rem=1; % Intial number for storing label of trivial states. 
for count_ref=1:1:Stop
    %the start position of V_ref=Height+Efield line
    Point_ref=(Band(count_ref).*e-Height)./(e.*E_contact); % unit m
    Cut_number=round(Point_ref./L*length(mesh_x));
    % Summation of wavefunction beyond V_ref=Height+Efield line
    test_count=sum(Vector(Cut_number:length(mesh_x),count_ref).^2);
    
    % if test_count > test_ref, it is a trivial states. 
    if test_count>test_ref
        Vector_remove(count_rem)= count_ref;
        count_rem=count_rem+1;
    end
    
end

%% Remove trivial eigen energies and wavefunctions. 
Band_remove_trivial(Vector_remove)=[];
Vector_remove_trivial=Vector(:,1:Stop);
Vector_remove_trivial(:,Vector_remove)=[];

%% Test the final eigen energies and wavefunctions. 
%% Trivial parts should be removed. 
figure(3)
plot(mesh_x*1e9,V./e,'color','k','linewidth',1)
hold on 

for count_remove_trivial=1:1:length(Band_remove_trivial)
    plot(mesh_x*1e9,Vector_remove_trivial(:,count_remove_trivial).^2*5+Band_remove_trivial(count_remove_trivial)...
        ,'color','k','linewidth',2)
hold on 
end

figure(4)
scatter([0:1:length(Band_remove_trivial)-1],Band_remove_trivial)

     
     
     