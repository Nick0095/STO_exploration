%%%% A CODE OF BESO CONBINED WITH EPSILON GREEDY POLICY FOR 3D CANTILEVER BY SUN. H and MA. L %%%%
function [ITER,C,IOU,C_difference,xsubopt,nsubopt]=epsilon_BESO(nelx,nely,nelz,volfrac,er,penal,rmin,E0,nu,emax,win,sym,xsubopt);
% USER-DEFINED LOOP PARAMETERS
maxloop = 200;    % Maximum number of iterations
tolx = 0.01;      % Termination criterion
displayflag = 1;  % Display structure flag
% USER-DEFINED MATERIAL PROPERTIES
Emin = 1e-9;      % Young's modulus of void-like material
unfixvol=1-volfrac/2;
nele = nelx*nely*nelz;
Nsubopt=10; 
C=zeros(Nsubopt,1);ITER=zeros(Nsubopt,1);C_difference=zeros(Nsubopt,1);IOU=zeros(Nsubopt,1);
nsubopt=[];
% USER-DEFINED LOAD DOFs
il = nelx; jl = nely/2; kl = nelz/2;                         % Coordinates
loadnid = kl*(nelx+1)*(nely+1)+il*(nely+1)+(nely+1-jl); % Node IDs
loaddof = 3*loadnid(:) - 1;                             % DOFs
% USER-DEFINED SUPPORT FIXED DOFs
[jf,kf] = meshgrid(1:nely+1,1:nelz+1);                  % Coordinates
fixednid = (kf-1)*(nely+1)*(nelx+1)+jf;                 % Node IDs
fixeddof = [3*fixednid(:); 3*fixednid(:)-1; 3*fixednid(:)-2]; % DOFs
% PREPARE FINITE ELEMENT ANALYSIS
ndof = 3*(nelx+1)*(nely+1)*(nelz+1);
F = sparse(loaddof,1,-1,ndof,1);
U = zeros(ndof,1);
freedofs = setdiff(1:ndof,fixeddof);
KE = lk_H8(nu);
nodegrd = reshape(1:(nely+1)*(nelx+1),nely+1,nelx+1);
nodeids = reshape(nodegrd(1:end-1,1:end-1),nely*nelx,1);
nodeidz = 0:(nely+1)*(nelx+1):(nelz-1)*(nely+1)*(nelx+1);
nodeids = repmat(nodeids,size(nodeidz))+repmat(nodeidz,size(nodeids));
edofVec = 3*nodeids(:)+1;
edofMat = repmat(edofVec,1,24)+ ...
    repmat([0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1 ...
    3*(nely+1)*(nelx+1)+[0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1]],nele,1);
iK = kron(edofMat,ones(24,1))';
jK = kron(edofMat,ones(1,24))';
% PREPARE FILTER
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

for z = 1 :10
    % INITIALIZE ITERATION
    x = ones(nely,nelx,nelz); loop = 0; change = 1; vol=1;
% START ITERATION
while ~(vol <= volfrac && change < tolx) && loop < maxloop 
    loop = loop+1;
    vol = max(vol*(1-er),volfrac);
    if loop >1; 
        olddc = dc; 
    end
    % FE-ANALYSIS
    sK = KE(:)*(Emin+x(:)'.^penal*(E0-Emin));
    K = sparse(iK(:),jK(:),sK(:)); K = (K+K')/2;
    U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:);
    % OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    c(loop) = 0.;
    ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),[nely,nelx,nelz]);
    c(loop) = 0.5*sum(sum(sum((Emin+x.^penal*(E0-Emin)).*ce)));
    dc = 0.5*penal*(E0-Emin)*x.^(penal-1).*ce;
    % FILTERING AND MODIFICATION OF SENSITIVITIES
    dc(:) = H*(dc(:)./Hs);  
    % STABLIZATION OF EVOLUTIONARY PROCESS
    if loop > 1; 
        dc = (dc+olddc)/2.; 
    end
    % The end of counting
    if z == 1 && (vol-volfrac)>0.01
        qi=loop+5;
    end   
    % RL / BESO DESIGN UPDATE
    if z > 1 && (vol-volfrac)>0.01
        e=emax+(0-emax)*min(1,1/(qi-5)*loop)^3;
        [x]=RL_ADDDEL(nelx,nely,nelz,unfixvol,volnext(loop),dc,e,win,sym);
    else
        [x] = ADDDEL(nele,vol,dc);
        x=reshape(x,nely,nelx,nelz);
    end    
    % PRINT RESULTS
    if loop>10;
        change=abs(sum(c(loop-9:loop-5))-sum(c(loop-4:loop)))/sum(c(loop-4:loop));
    end
    if z == 1 && (vol-volfrac)>0.01
        volnext(loop)=mean(x(:));
    end
    fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,c(loop),mean(x(:)),change);
    % PLOT DENSITIES
    if displayflag, figure(z);clf; display_3D(x,1); end
end
%% evaluation
C(z)=c(loop);ITER(z)=loop;
C_difference(z)=C(z)/C(1)-1;
if size(xsubopt,1) == 0
    x=reshape(x,1,nele);
    xsubopt=[xsubopt;x]; nsubopt=[nsubopt;z];
else
    if C_difference(z) < 1.0
        x=reshape(x,1,nele);
        for j=1:size(xsubopt,1)
            IOU(z)=max(IOU(z),length(find((x+xsubopt(j,:))>1.5))/length(find((x+xsubopt(j,:))>0.5)));
            if IOU(z)>0.9
                break;
            end
            if j==size(xsubopt,1)
                xsubopt=[xsubopt;x]; nsubopt=[nsubopt;z];
            end
        end
    end
end
end
end

%% ===================== AUXILIARY FUNCTIONS ===============================
% OPTIMALITY CRITERIA UPDATE
function [x]=ADDDEL(nele,vol,dc)
    l1 = min(min(min(dc)));    l2 = max(max(max(dc)));
    while ((l2-l1)/l2 > 1.0e-5)
        th = (l1+l2)/2.0;
        x = max(0.001,sign(dc-th));
        if mean(x(:))-vol > 0;
            l1 = th;
        else
            l2 = th;
        end
    end
end

%%%%%%%%%% RL UPDATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x]=RL_ADDDEL(nelx,nely,nelz,unfixvol,volnext,dc,e,win,sym)
% symmetric constraint
if sym==1
    nely=nely/2;
    nelz=nelz/2;
    dc=dc(1:nely,:,1:nelz);
end 

% determine V by both V and dc
nele = nelx*nely*nelz;
dc=reshape(dc,1,nele);
[dc0,index]=sort(dc);

% delete 1-e elements
base=floor((1-e)*(1-volnext)*nele);
x=ones(1,nele);
dc_critical=dc0(base); % RANDOM CHOICE
index_critical=find(dc0==dc_critical);
base2=index_critical(1);
index(index_critical)= index(index_critical(1)-1+randperm(length(index_critical)));
x(index(1:base))=0.001;
x=reshape(x,nely,nelx,nelz);
% delete e elements by window
m=floor(e*(1-volnext)*nele/(2*win+1)^3);
if m>0
    n=unfixvol*nele-base;
    p=randperm(n);
    for i=1:m
        ind=index(base+p(i));
        indz=ceil(ind/nelx/nely);
        indx=ceil(mod(ind,nelx*nely)/nely);
        indy=mod(mod(ind,nelx*nely),nely);
        if indx==0
            indx=nelx;
        end
        if indy==0
            indy=nely;
        end
        for k=1:win
            x(indy,max(1,indx-win):min(indx+win,nelx),indz)=0.001;
            for j=1:winy
                x(max(1,indy-j),max(1,indx-win):min(indx+win,nelx),indz)=0.001;
                x(min(nely,indy+j),max(1,indx-win):min(indx+win,nelx),indz)=0.001;
            end
        end
    end
end
% delete rest elements
n=floor((mean(x(:))-volnext)*nele);
if n>0
    ind2=find(x(index(1:unfixvol*nele))==1);
    p=randperm(length(ind2));
    x(index(ind2(p(1:n))))=0.001;
end
if sym==1
    x=[x;flipud(x)];
    x=cat(3,x,flip(x,3));
end 
end

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
function display_3D(rho,n)
[nely,nelx,nelz] = size(rho);
hx = 1; hy = 1; hz = 1;            % User-defined unit element size
face = [1 2 3 4; 2 6 7 3; 4 3 7 8; 1 5 8 4; 1 2 6 5; 5 6 7 8];
set(gcf,'Name',sprintf('ISO display %i',n),'NumberTitle','off');
for k = 1:nelz
    z = (k-1)*hz;
    for i = 1:nelx
        x = (i-1)*hx;
        for j = 1:nely
            y = nely*hy - (j-1)*hy;
            if (rho(j,i,k) > 0.5)  % User-defined display density threshold
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
