function [Gaussian_point weight]=Gaussian_quadrature_3D(order)
%% 
switch order
    
    case 'first'
    
    Gaussian_point=zeros(6,3);
    Gaussian_point(1,:)=[ 1 0 0];
    Gaussian_point(1,:)=[-1 0 0];
    Gaussian_point(1,:)=[ 0 1 0];
    Gaussian_point(1,:)=[ 0-1 0];
    Gaussian_point(1,:)=[ 0 0 1];
    Gaussian_point(1,:)=[ 0 0-1];
    weight=[1;1;1;1;1;1;1;1]*4/3; 
    case 'second'
     
     Gaussian_point=zeros(8,3);        
     Gaussian_point(1,:)=[ 1/sqrt(3)  1/sqrt(3)  -1/sqrt(3)];
     Gaussian_point(2,:)=[-1/sqrt(3)  1/sqrt(3)  -1/sqrt(3)];
     Gaussian_point(3,:)=[-1/sqrt(3) -1/sqrt(3)  -1/sqrt(3)];
     Gaussian_point(4,:)=[ 1/sqrt(3) -1/sqrt(3)  -1/sqrt(3)];
     Gaussian_point(5,:)=[ 1/sqrt(3)  1/sqrt(3)   1/sqrt(3)];
     Gaussian_point(6,:)=[-1/sqrt(3)  1/sqrt(3)   1/sqrt(3)];
     Gaussian_point(7,:)=[-1/sqrt(3) -1/sqrt(3)   1/sqrt(3)];
     Gaussian_point(8,:)=[ 1/sqrt(3) -1/sqrt(3)   1/sqrt(3)];
     weight=[1;1;1;1;1;1;1;1];
     
    case 'third'
    error('under construction')
end
