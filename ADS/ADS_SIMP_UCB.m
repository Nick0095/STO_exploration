%%%% A CODE OF SIMP CONBINED WITH UCB FOR ADS BY SUN. H and MA. L %%%%
function [ITER,C,IOU,C_difference,xsubopt,nsubopt,maxstressz] = ADS_SIMP_UCB(volfrac,penal,rmin,rho,xsubopt)
% USER-DEFINED LOOP PARAMETERS
maxloop = 200;    % Maximum number of iterations
tolx = 0.01;      % Termination criterion
displayflag = 1;  % Display structure flag
% USER-DEFINED MATERIAL PROPERTIES
E0 = 68.9e9 * 2;           % Young's modulus of solid material
Emin = 1e-9;      % Young's modulus of void-like material
nu = 0.3;         % Poisson's ratio
coef = 9e6;
Nsubopt=10;
%% Define coordinates of elements and nodes
load('coordinate_node.mat');
X = coordinate_node(1,:);
Y = coordinate_node(2,:);
Z = coordinate_node(3,:);
cenx = min(X);
ceny = min(Y);
cenz = max(Z);
X = X - cenx;
Y = Y - ceny;
Z = -(Z - cenz);
X = round(X/20) + 1;
Y = round(Y/20) + 1;
Z = round(Z/20) + 1;
new_coor = unique([X; Y; Z]','rows','stable');
X = new_coor(:,1);
Y = new_coor(:,2);
Z = new_coor(:,3);
nelx = max(X); nely = max(Y); nelz = max(Z);
X = [X; nelx]; Y = [Y; 1]; Z = [Z; 1];
n = length(X);

nodeX = [X; X; X+1; X+1; X; X; X+1; X+1]';
nodeY = [Y; Y; Y; Y; Y+1; Y+1; Y+1; Y+1]';
nodeZ = [Z; Z+1; Z; Z+1; Z; Z+1; Z; Z+1]';
new_coor = unique([nodeX; nodeY; nodeZ]','rows','stable');
nodeX = new_coor(:,1);
nodeY = new_coor(:,2);
nodeZ = new_coor(:,3);
n2 = length(nodeX);

%% Solid elements
nele = nelx*nely*nelz;
solidele = (Z - 1)*nelx*nely + (X - 1)*nely + Y;
voidele = setdiff(1:nele,solidele);

%% elements and nodes
element = zeros(nele, 9);
for k = 1 : nelz
    for j = 1 : nely
        for i = 1 : nelx
            element((k-1)*nelx*nely+(i-1)*nely+j,:) = [(k-1)*nelx*nely+(i-1)*nely+j ...
                (k-1)*(nelx+1)*(nely+1)+i*(nely+1)+j+1 (k-1)*(nelx+1)*(nely+1)+(i-1)*(nely+1)+j+1 ...
                (k-1)*(nelx+1)*(nely+1)+(i-1)*(nely+1)+j (k-1)*(nelx+1)*(nely+1)+i*(nely+1)+j ...
                k*(nelx+1)*(nely+1)+i*(nely+1)+j+1 k*(nelx+1)*(nely+1)+(i-1)*(nely+1)+j+1 ...
                k*(nelx+1)*(nely+1)+(i-1)*(nely+1)+j k*(nelx+1)*(nely+1)+i*(nely+1)+j];
        end
    end
end
node = zeros((nelx+1)*(nely+1)*(nelz+1), 4);
for k = 1:nelz + 1
    for i = 1:nelx + 1
        for j = 1:nely + 1
            node((k-1)*(nelx+1)*(nely+1) + (i-1)*(nely+1) + j, :) = [(k-1)*(nelx+1)*(nely+1) + (i-1)*(nely+1) + j, i-1, j-1, k-1];
        end
    end
end

%% USER-DEFINED SUPPORT FIXED DOFs
% Plane1: view
theta1 = 72;
d = 400;
d1 = -tand(theta1)*(nodeX-1) + (nodeY-1);
d2 = (((nodeY-1) - (520+20)/20)/(d/20/2*sind(theta1))).^2 + ((nodeZ-1)/(d/20/2)).^2;
d1min = (520 + 55*cosd(theta1))/20 - tand(theta1)*(576/2 - 123.5 - 55*sind(theta1))/20;
d1max = (520 + 85*cosd(theta1))/20 - tand(theta1)*(576/2 - 123.5 - 85*sind(theta1))/20;
index1 = find(d1<d1max & d1>d1min & d2<=1);
% Plane2: arm
alpha = 4;
d1 = -cosd(52)*cosd(alpha)*(nodeX-1) + sind(alpha)*(nodeY-1) + sind(52)*cosd(alpha)*(nodeZ-1);
d2 = (190 + 280*sind(4) + 294/2*cosd(4) + 10)/20; d3 = (190 + 280*sind(4) - 294/2*cosd(4) - 10)/20;
d4 = (((nodeX-1.5) - (576/2 - 280*cosd(alpha)*cosd(52))/20)/(294/20/2*cosd(38))).^2 + (((nodeY-1.5)-(190+280*sind(4))/20)/(294/20/2*cosd(4))).^2;
d1min = (cosd(52)*cosd(alpha)*(576/2 - 280*cosd(alpha)*cosd(52)) + sind(alpha)*(190 + 280*sind(alpha)) + sind(52)*cosd(alpha)*280*cosd(alpha)*sind(52))/20;
index2 = find(abs(d1-d1min)<6.6 & (nodeY-1)<=d2 & (nodeY-1)>=d3 & d4<=0.8);
% Plane3: bottom
index3 = find(nodeY==1);
index = [index1; index2; index3];
fixednid = (nodeZ(index) - 1)*(nelx+1)*(nely+1) + (nodeX(index) - 1)*(nely+1) + nodeY(index);
fixedindex = zeros((nelx+1)*(nely+1)*(nelz+1),1);
fixedindex(fixednid) = 1; 
fixeddof = [3*fixednid(:); 3*fixednid(:)-1; 3*fixednid(:)-2]; % DOFs
% Plane4: x-y plane symmetry
index4 = find(nodeZ==1);
fixednid2 = (nodeZ(index4) - 1)*(nelx+1)*(nely+1) + (nodeX(index4) - 1)*(nely+1) + nodeY(index4);
fixedindex(fixednid2) = 1; 
fixeddof = [fixeddof; fixednid2(:)]; % DOFs

boundary_displace=[fixednid(:) 3*ones(length(fixednid(:)),1) zeros(length(fixednid(:)),1); ...
    fixednid(:) ones(length(fixednid(:)),1) zeros(length(fixednid(:)),1); ...
    fixednid(:) 2*ones(length(fixednid(:)),1) zeros(length(fixednid(:)),1); ...
    fixednid2(:) 3*ones(length(fixednid2(:)),1) zeros(length(fixednid2(:)),1)];
surface_traction=zeros(1,9);

%% Fixed elements
fixedele = [];
for i = 1:n
    if ~isempty(intersect(element(solidele(i),2:end),fixednid))
        fixedele = [fixedele; solidele(i)];
    end
end

%% USER-DEFINED LOAD DOFs
%% find inner nodes
x = zeros(nely, nelx, nelz);
border = zeros(nelx,nely);
for i = 1:n
    x(Y(i),X(i),Z(i)) = 1;
    border(X(i),Y(i)) = max(border(X(i),Y(i)),Z(i));
end

% part1: border
for j = 1:nely
    for i = 1:nelx
        if border(i,j) > 0
            x(j,i,1:border(i,j)) = x(j,i,1:border(i,j)) + 0.2;
        end
    end
end

% part2: hole: plane x-y
x(7,6:8,1:nelz) = x(7,6:8,1:nelz) + 0.2; 
x(8,5:9,1:nelz) = x(8,5:9,1:nelz) + 0.2; 
x(9,4:9,1:nelz) = x(9,4:9,1:nelz) + 0.2; 
x(10,4:10,1:nelz) = x(10,4:10,1:nelz) + 0.2; 
x(11,4:10,1:nelz) = x(11,4:10,1:nelz) + 0.2; 
x(12,4:10,1:nelz) = x(12,4:10,1:nelz) + 0.2; 
x(13,4:10,1:nelz) = x(13,4:10,1:nelz) + 0.2; 
x(14,4:10,1:nelz) = x(14,4:10,1:nelz) + 0.2; 
x(15,4:9,1:nelz) = x(15,4:9,1:nelz) + 0.2; 
x(16,6:8,1:nelz) = x(16,6:8,1:nelz) + 0.2; 

% figure(1); x(x>1)=1; clf; display_3D(x);
% part3: redundancy
x(8:14,1,1:nelz) = floor(x(8:14,1,1:nelz));
x([6 14:16],2,1:nelz) = floor(x([6 14:16],2,1:nelz));
x(16:19,3,1:nelz) = floor(x(16:19,3,1:nelz));
x(18,4,1:nelz) = floor(x(18,4,1:nelz));
% figure(2); x(x>1)=1; clf; display_3D(x);

innerele = find(x > 0 & x < 1);

%% loading
x0 = zeros(nely,nelx,nelz);
x0(innerele) = 0.2;
x0(solidele) = 1;
% x0(fixedele) = 0.2;
x0(16,3,4:5) = 0;
% figure(1); clf; display_3D(x0); % test
xFull = zeros(size(x0) + 2);
xFull(:,:,1) = 0.2;
xFull(2:nely+1, 2:nelx+1, 2:nelz+1) = x0;
diff1 = circshift(xFull,[0,-1,0]) - xFull;
diff2 = circshift(xFull,[0,1,0]) - xFull;
diff3 = circshift(xFull,[-1,0,0]) - xFull;
diff4 = circshift(xFull,[1,0,0]) - xFull;
diff5 = circshift(xFull,[0,0,-1]) - xFull;
diff6 = circshift(xFull,[0,0,1]) - xFull;
loadx = zeros(size(x0));
loady = zeros(size(x0));
loadz = zeros(size(x0));
F0 = zeros(3*(nelx+1)*(nely+1)*(nelz+1),1);
loadedele = [];
for i = 1:nelx
    for j = 1:nely
        for k = 1:nelz
            % detect if fixed
            if ~isempty(intersect((k - 1)*nelx*nely + (i - 1)*nely + j, fixedele))
                continue;
            end              
            flag = 0;
            if diff1(j+1,i+1,k+1) == -1
                flag = 1;
                F0(3*((k-1)*(nelx+1)*(nely+1) + i*(nely+1) + j) - 2) = F0(3*((k-1)*(nelx+1)*(nely+1) + i*(nely+1) + j) - 2) - 0.25;
                F0(3*(k*(nelx+1)*(nely+1) + i*(nely+1) + j) - 2) = F0(3*(k*(nelx+1)*(nely+1) + i*(nely+1) + j) - 2) - 0.25;
                F0(3*((k-1)*(nelx+1)*(nely+1) + i*(nely+1) + j+1) - 2) = F0(3*((k-1)*(nelx+1)*(nely+1) + i*(nely+1) + j+1) - 2) - 0.25;
                F0(3*(k*(nelx+1)*(nely+1) + i*(nely+1) + j+1) - 2) = F0(3*(k*(nelx+1)*(nely+1) + i*(nely+1) + j+1) - 2) - 0.25;
            end
            if diff2(j+1,i+1,k+1) == -1
                flag = 1;
                F0(3*((k-1)*(nelx+1)*(nely+1) + (i-1)*(nely+1) + j) - 2) = F0(3*((k-1)*(nelx+1)*(nely+1) + (i-1)*(nely+1) + j) - 2) + 0.25;
                F0(3*(k*(nelx+1)*(nely+1) + (i-1)*(nely+1) + j) - 2) = F0(3*(k*(nelx+1)*(nely+1) + (i-1)*(nely+1) + j) - 2) + 0.25;
                F0(3*((k-1)*(nelx+1)*(nely+1) + (i-1)*(nely+1) + j+1) - 2) = F0(3*((k-1)*(nelx+1)*(nely+1) + (i-1)*(nely+1) + j+1) - 2) + 0.25;
                F0(3*(k*(nelx+1)*(nely+1) + (i-1)*(nely+1) + j+1) - 2) = F0(3*(k*(nelx+1)*(nely+1) + (i-1)*(nely+1) + j+1) - 2) + 0.25;             
            end
            if diff3(j+1,i+1,k+1) == -1
                flag = 1;
                F0(3*((k-1)*(nelx+1)*(nely+1) + (i-1)*(nely+1) + j+1) - 1) = F0(3*((k-1)*(nelx+1)*(nely+1) + (i-1)*(nely+1) + j+1) - 1) - 0.25;
                F0(3*(k*(nelx+1)*(nely+1) + (i-1)*(nely+1) + j+1) - 1) = F0(3*(k*(nelx+1)*(nely+1) + (i-1)*(nely+1) + j+1) - 1) - 0.25;
                F0(3*((k-1)*(nelx+1)*(nely+1) + i*(nely+1) + j+1) - 1) = F0(3*((k-1)*(nelx+1)*(nely+1) + i*(nely+1) + j+1) - 1) - 0.25;
                F0(3*(k*(nelx+1)*(nely+1) + i*(nely+1) + j+1) - 1) = F0(3*(k*(nelx+1)*(nely+1) + i*(nely+1) + j+1) - 1) - 0.25;
            end
            if diff4(j+1,i+1,k+1) == -1
                flag = 1;
                F0(3*((k-1)*(nelx+1)*(nely+1) + (i-1)*(nely+1) + j) - 1) = F0(3*((k-1)*(nelx+1)*(nely+1) + (i-1)*(nely+1) + j) - 1) + 0.25;
                F0(3*(k*(nelx+1)*(nely+1) + (i-1)*(nely+1) + j) - 1) = F0(3*(k*(nelx+1)*(nely+1) + (i-1)*(nely+1) + j) - 1) + 0.25;
                F0(3*((k-1)*(nelx+1)*(nely+1) + i*(nely+1) + j) - 1) = F0(3*((k-1)*(nelx+1)*(nely+1) + i*(nely+1) + j) - 1) + 0.25;
                F0(3*(k*(nelx+1)*(nely+1) + i*(nely+1) + j) - 1) = F0(3*(k*(nelx+1)*(nely+1) + i*(nely+1) + j) - 1) + 0.25; 
            end 
            if diff5(j+1,i+1,k+1) == -1
                flag = 1;
                F0(3*(k*(nelx+1)*(nely+1) + (i-1)*(nely+1) + j)) = F0(3*(k*(nelx+1)*(nely+1) + (i-1)*(nely+1) + j)) - 0.25;
                F0(3*(k*(nelx+1)*(nely+1) + (i-1)*(nely+1) + j+1)) = F0(3*(k*(nelx+1)*(nely+1) + (i-1)*(nely+1) + j+1)) - 0.25;
                F0(3*(k*(nelx+1)*(nely+1) + i*(nely+1) + j)) = F0(3*(k*(nelx+1)*(nely+1) + i*(nely+1) + j)) - 0.25;
                F0(3*(k*(nelx+1)*(nely+1) + i*(nely+1) + j+1)) = F0(3*(k*(nelx+1)*(nely+1) + i*(nely+1) + j+1)) - 0.25;
            end
            if diff6(j+1,i+1,k+1) == -1
                flag = 1;
                F0(3*((k-1)*(nelx+1)*(nely+1) + (i-1)*(nely+1) + j)) =  F0(3*((k-1)*(nelx+1)*(nely+1) + (i-1)*(nely+1) + j))+ 0.25;
                F0(3*((k-1)*(nelx+1)*(nely+1) + (i-1)*(nely+1) + j+1)) = F0(3*((k-1)*(nelx+1)*(nely+1) + (i-1)*(nely+1) + j+1)) + 0.25;
                F0(3*((k-1)*(nelx+1)*(nely+1) + i*(nely+1) + j)) = F0(3*((k-1)*(nelx+1)*(nely+1) + i*(nely+1) + j)) + 0.25;
                F0(3*((k-1)*(nelx+1)*(nely+1) + i*(nely+1) + j+1)) = F0(3*((k-1)*(nelx+1)*(nely+1) + i*(nely+1) + j+1)) + 0.25;
            end
            if flag == 1
                loadedele = [loadedele; (k-1)*nelx*nely + (i-1)*nely + j, i, j, k];
            end
        end
    end
end
indexx = find(F0(1:3:end)~=0);
indexy = find(F0(2:3:end)~=0);
indexz = find(F0(3:3:end)~=0);
point_force=[indexx ones(length(indexx),1) coef*F0(3*(indexx-1)+1); ...
    indexy 2*ones(length(indexy),1) coef*F0(3*(indexy-1)+2); ...
    indexz 3*ones(length(indexz),1) coef*F0(3*(indexz-1)+3)];
x0 = zeros(nely,nelx,nelz);
x0(solidele) = 1;
% figure(2); clf; display_3D(x0);  % test
x0(loadedele(:,1)) = 0.2;
% figure(3); clf; display_3D(x0);  % test

%% PREPARE FILTER
iH = ones(nele*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for k1 = 1:nelz
    for i1 = 1:nelx
        for j1 = 1:nely
            e1 = (k1-1)*nelx*nely + (i1-1)*nely+j1;
            for k2 = max(k1-(ceil(rmin)-1),1):min(k1+(ceil(rmin)-1),nelz)
                for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
                    for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
                        e2 = (k2-1)*nelx*nely + (i2-1)*nely+j2;
                        k = k+1;
                        iH(k) = e1;
                        jH(k) = e2;
                        sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2+(k1-k2)^2));
                    end
                end
            end
        end
    end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);

Na = ones(nele,1); iter = 0;
C=zeros(Nsubopt,1);ITER=zeros(Nsubopt,1);C_difference=zeros(Nsubopt,1);IOU=zeros(Nsubopt,1);maxstressz=zeros(Nsubopt,1);maxstressele=zeros(Nsubopt,1);
nsubopt=[]; 
for z = 1:Nsubopt
% INITIALIZE ITERATION
x = volfrac * ones(nely,nelx,nelz);
x(voidele) = 0;
xnew = x;
xPhys = x; 
loop = 0; 
change = 1;
% START ITERATION
while change > tolx && loop < maxloop
    loop = loop+1;
    [displacement strain stress von_stress c(loop) dc0] = Mechanical_analysis(element(solidele,:),node,boundary_displace,point_force,F0,surface_traction,xPhys(solidele),loadedele,penal,E0,Emin,nu);    
    dc = zeros(nele, 1);
    dc(solidele) = dc0;
    maxstress(loop) = max(von_stress);
    maxstressz(z) = max(von_stress);
    % OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    dv = ones(nely,nelx,nelz);
    % FILTERING AND MODIFICATION OF SENSITIVITIES
    dc(:) = H*(x(:).*dc(:))./Hs./max(1e-3,x(:));
  % The end of counting
    if z == 1
        qi = loop;
    end 
    % OPTIMALITY CRITERIA UPDATE
    if z > 1  
        iter = iter+1;
        dc(loadedele(:,1)) = dc(loadedele(:,1)) + rho * sqrt(log(iter+1)/2./Na(loadedele(:,1)));
    end     
    l1 = 0; l2 = 1e9; move = 0.2;
    while (l2-l1)/(l1+l2) > 1e-3
        lmid = 0.5*(l2+l1);
        xnew(loadedele(:,1)) = max(0.2, max(x(loadedele(:,1))-move,min(1,min(x(loadedele(:,1))+move,x(loadedele(:,1)).*sqrt(dc(loadedele(:,1))./dv(loadedele(:,1))/lmid)))));
        xPhys = xnew;
        if mean(xPhys(loadedele(:,1))) > volfrac, l1 = lmid; else l2 = lmid; end
    end
    if loop>5;
        change=abs(sum(c(loop-5:loop-3))-sum(c(loop-2:loop)))/sum(c(loop-2:loop));
    end 
    % Update Sa, Na
    Na(loadedele(:,1)) = Na(loadedele(:,1)) + max(0, xnew(loadedele(:,1)) - x(loadedele(:,1)));
    x = xnew;
    % PRINT RESULTS
    fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,c(loop),mean(xPhys(:)),change);
    % PLOT DENSITIES
%     figure(z);
%     if displayflag, clf; display_3D(xPhys); end
end
figure(z);
clf; display_3D(xPhys);

%% evaluation
C(z)=c(loop); ITER(z)=loop;
C_difference(z)=C(z)/C(1)-1;
if size(xsubopt,1) == 0
    x=reshape(x,1,nele);
    xsubopt=[xsubopt;x]; nsubopt=[nsubopt;z];
else
    if C_difference(z) < 10
        x=reshape(x,1,nele);
        for j=1:size(xsubopt,1)
            IOU(z)=max(IOU(z),sum(min(x(loadedele(:,1)),xsubopt(j,loadedele(:,1))))/sum(max(x(loadedele(:,1)),xsubopt(j,loadedele(:,1)))));
            if IOU(z)>0.95
                break;
            end
            if j==size(xsubopt,1)
                xsubopt=[xsubopt;x]; nsubopt=[nsubopt;z];
            end
        end
    end
end
% if size(xsubopt,1)==10
%     break
% end
end
end
% ===================== AUXILIARY FUNCTIONS ===============================
% GENERATE ELEMENT STIFFNESS MATRIX
function [KE] = lk_H8(nu)
A = [32 6 -8 6 -6 4 3 -6 -10 3 -3 -3 -4 -8;
    -48 0 0 -24 24 0 0 0 12 -12 0 12 12 12];
k = 1/72*A'*[1; nu];
% GENERATE SIX SUB-MATRICES AND THEN GET KE MATRIX
K1 = [k(1) k(2) k(2) k(3) k(5) k(5);
    k(2) k(1) k(2) k(4) k(6) k(7);
    k(2) k(2) k(1) k(4) k(7) k(6);
    k(3) k(4) k(4) k(1) k(8) k(8);
    k(5) k(6) k(7) k(8) k(1) k(2);
    k(5) k(7) k(6) k(8) k(2) k(1)];
K2 = [k(9)  k(8)  k(12) k(6)  k(4)  k(7);
    k(8)  k(9)  k(12) k(5)  k(3)  k(5);
    k(10) k(10) k(13) k(7)  k(4)  k(6);
    k(6)  k(5)  k(11) k(9)  k(2)  k(10);
    k(4)  k(3)  k(5)  k(2)  k(9)  k(12)
    k(11) k(4)  k(6)  k(12) k(10) k(13)];
K3 = [k(6)  k(7)  k(4)  k(9)  k(12) k(8);
    k(7)  k(6)  k(4)  k(10) k(13) k(10);
    k(5)  k(5)  k(3)  k(8)  k(12) k(9);
    k(9)  k(10) k(2)  k(6)  k(11) k(5);
    k(12) k(13) k(10) k(11) k(6)  k(4);
    k(2)  k(12) k(9)  k(4)  k(5)  k(3)];
K4 = [k(14) k(11) k(11) k(13) k(10) k(10);
    k(11) k(14) k(11) k(12) k(9)  k(8);
    k(11) k(11) k(14) k(12) k(8)  k(9);
    k(13) k(12) k(12) k(14) k(7)  k(7);
    k(10) k(9)  k(8)  k(7)  k(14) k(11);
    k(10) k(8)  k(9)  k(7)  k(11) k(14)];
K5 = [k(1) k(2)  k(8)  k(3) k(5)  k(4);
    k(2) k(1)  k(8)  k(4) k(6)  k(11);
    k(8) k(8)  k(1)  k(5) k(11) k(6);
    k(3) k(4)  k(5)  k(1) k(8)  k(2);
    k(5) k(6)  k(11) k(8) k(1)  k(8);
    k(4) k(11) k(6)  k(2) k(8)  k(1)];
K6 = [k(14) k(11) k(7)  k(13) k(10) k(12);
    k(11) k(14) k(7)  k(12) k(9)  k(2);
    k(7)  k(7)  k(14) k(10) k(2)  k(9);
    k(13) k(12) k(10) k(14) k(7)  k(11);
    k(10) k(9)  k(2)  k(7)  k(14) k(7);
    k(12) k(2)  k(9)  k(11) k(7)  k(14)];
KE = 1/((nu+1)*(1-2*nu))*...
    [ K1  K2  K3  K4;
    K2'  K5  K6  K3';
    K3' K6  K5' K2';
    K4  K3  K2  K1'];
end
% DISPLAY 3D TOPOLOGY (ISO-VIEW)
function display_3D(rho)
[nely,nelx,nelz] = size(rho);
hx = 1; hy = 1; hz = 1;            % User-defined unit element size
face = [1 2 3 4; 2 6 7 3; 4 3 7 8; 1 5 8 4; 1 2 6 5; 5 6 7 8];
set(gcf,'Name','ISO display','NumberTitle','off');
for k = 1:nelz
    z = (k-1)*hz;
    for i = 1:nelx
        x = (i-1)*hx;
        for j = 1:nely
            y = nely*hy - (j-1)*hy;
            if (rho(j,i,k) > 0.1)  % User-defined display density threshold
                vert = [x y z; x y-hx z; x+hx y-hx z; x+hx y z; x y z+hx;x y-hx z+hx; x+hx y-hx z+hx;x+hx y z+hx];
                vert(:,[2 3]) = vert(:,[3 2]); vert(:,2,:) = -vert(:,2,:);
                patch('Faces',face,'Vertices',vert,'FaceColor',[0.2+0.8*(1-rho(j,i,k)),0.2+0.8*(1-rho(j,i,k)),0.2+0.8*(1-rho(j,i,k))]);
                hold on;
            end
        end
    end
end
axis equal; axis tight; axis off; box on; view([30,30]); pause(1e-6);
end
