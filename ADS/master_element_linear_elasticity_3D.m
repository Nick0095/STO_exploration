function [B N N_diff J H J_inv x y z]= master_element_linear_elasticity_3D(order,ksi,eta,gamma,coordinate)
%% this function calculate the reshape function, reshape devirative , Jocabian in the master element

switch order 
    
    case 'linear'
    N1        = zeros(8,1);
    N1_diff   = zeros(8,3);
    
    N         = zeros(24,3);
    N_diff    = zeros(9,24);
    B = zeros(6,24);
    
    J         = zeros(9,9);
    J_inv     = zeros(9,9);
   
    N1(1)     = 0.125*(1+ksi)*(1+eta)*(1-gamma);
    N1(2)     = 0.125*(1-ksi)*(1+eta)*(1-gamma);
    N1(3)     = 0.125*(1-ksi)*(1-eta)*(1-gamma);
    N1(4)     = 0.125*(1+ksi)*(1-eta)*(1-gamma);
    N1(5)     = 0.125*(1+ksi)*(1+eta)*(1+gamma);
    N1(6)     = 0.125*(1-ksi)*(1+eta)*(1+gamma);
    N1(7)     = 0.125*(1-ksi)*(1-eta)*(1+gamma);
    N1(8)     = 0.125*(1+ksi)*(1-eta)*(1+gamma);
    
    N1_diff(1,:) = [ 0.125*(1+eta)*(1-gamma)  0.125*(1+ksi)*(1-gamma) -0.125*(1+ksi)*(1+eta)];
    N1_diff(2,:) = [-0.125*(1+eta)*(1-gamma)  0.125*(1-ksi)*(1-gamma) -0.125*(1-ksi)*(1+eta)];
    N1_diff(3,:) = [-0.125*(1-eta)*(1-gamma) -0.125*(1-ksi)*(1-gamma) -0.125*(1-ksi)*(1-eta)];
    N1_diff(4,:) = [ 0.125*(1-eta)*(1-gamma) -0.125*(1+ksi)*(1-gamma) -0.125*(1+ksi)*(1-eta)];
    N1_diff(5,:) = [ 0.125*(1+eta)*(1+gamma)  0.125*(1+ksi)*(1+gamma)  0.125*(1+ksi)*(1+eta)];
    N1_diff(6,:) = [-0.125*(1+eta)*(1+gamma)  0.125*(1-ksi)*(1+gamma)  0.125*(1-ksi)*(1+eta)];
    N1_diff(7,:) = [-0.125*(1-eta)*(1+gamma) -0.125*(1-ksi)*(1+gamma)  0.125*(1-ksi)*(1-eta)];
    N1_diff(8,:) = [ 0.125*(1-eta)*(1+gamma) -0.125*(1+ksi)*(1+gamma)  0.125*(1+ksi)*(1-eta)];

    N(1:3:24,1)=N1;N(2:3:24,2)=N1;N(3:3:24,3)=N1;
    
    B(1,1:3:24)=N1_diff(:,1)'; B(2,2:3:24)=N1_diff(:,2)'; B(3,3:3:24)=N1_diff(:,3)'; 
    B(4,2:3:24)=N1_diff(:,1)'; B(4,1:3:24)=N1_diff(:,2)'; B(5,2:3:24)=N1_diff(:,3)'; 
    B(6,3:3:24)=N1_diff(:,1)'; B(5,3:3:24)=N1_diff(:,2)'; B(6,1:3:24)=N1_diff(:,3)'; 
    
    N_diff(1,1:3:24)=N1_diff(:,1)';N_diff(2,1:3:24)=N1_diff(:,2)';N_diff(3,1:3:24)=N1_diff(:,3)'; 
    N_diff(4,2:3:24)=N1_diff(:,1)';N_diff(5,2:3:24)=N1_diff(:,2)';N_diff(6,2:3:24)=N1_diff(:,3)';
    N_diff(7,3:3:24)=N1_diff(:,1)';N_diff(8,3:3:24)=N1_diff(:,2)';N_diff(9,3:3:24)=N1_diff(:,3)';    
    
    H=zeros(6,9);H(1,1)=1;H(2,5)=1;H(3,9)=1;H(4,[2 4])=1;H(5,[6 8])=1;H(6,[3 7])=1;
    J1=N1_diff'*coordinate;
    J(1:3,1:3)=J1;J(4:6,4:6)=J1;J(7:9,7:9)=J1;
    J_inv(1:3,1:3)=inv(J1);J_inv(4:6,4:6)=inv(J1);J_inv(7:9,7:9)=inv(J1);
    
    x=N1'*coordinate(:,1);    
    y=N1'*coordinate(:,2);
    z=N1'*coordinate(:,3);
    
    case 'quadratic'
    error('under construction')

    case 'cubic'
    error('under construction')
end

end
