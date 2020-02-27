%%%% A CODE OF SIMP CONBINED WITH UCB FOR 2D CANTILEVER BY SUN. H and MA. L %%%%
function [ITER,C,IOU,C_difference,xsubopt,nsubopt] = SIMP_UCB_cantilever2d(nelx,nely,volfrac,penal,rmin,ft,rho,xsubopt)
Nsubopt=10; Sa = zeros(nely*nelx,1); Na = ones(nely*nelx,1);
maxloop = 200;
C=zeros(Nsubopt,1);ITER=zeros(Nsubopt,1);C_difference=zeros(Nsubopt,1);IOU=zeros(Nsubopt,1);
nsubopt=[];
%% MATERIAL PROPERTIES
E0 = 1;
Emin = 1e-9;
nu = 0.3;
%% PREPARE FINITE ELEMENT ANALYSIS
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
U = zeros(2*(nely+1)*(nelx+1),1);
F = sparse(2*(nelx+1)*(nely+1)-nely,1,-1,2*(nely+1)*(nelx+1),1); %right middle
fixeddofs=[1:2*(nely+1)]; %left
alldofs = [1:2*(nely+1)*(nelx+1)];
freedofs = setdiff(alldofs,fixeddofs);
%% PREPARE FILTER
iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for i1 = 1:nelx
  for j1 = 1:nely
    e1 = (i1-1)*nely+j1;
    for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
      for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
        e2 = (i2-1)*nely+j2;
        k = k+1;
        iH(k) = e1;
        jH(k) = e2;
        sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
      end
    end
  end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);

for z = 1 :Nsubopt
%% INITIALIZE ITERATION
x = repmat(volfrac,nely,nelx);
xPhys = x;
loop = 0;
change = 1;
%% START ITERATION
while change > 0.01 && loop < maxloop
  loop = loop + 1;
  %% FE-ANALYSIS
  sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
  K = sparse(iK,jK,sK); K = (K+K')/2;
  U(freedofs) = K(freedofs,freedofs)\F(freedofs);
  %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
  ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
  c(loop) = 0.5*sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce));
  dc = -0.5*penal*(E0-Emin)*xPhys.^(penal-1).*ce;
  dv = ones(nely,nelx);
    if z > 1  & loop <= qi
        dc(index0) = dc(index0) - rho * sqrt(log(loop+1)/2./Na(index0)');
        dc(index1) = min(dc(:));
    end
  %% FILTERING/MODIFICATION OF SENSITIVITIES
  if ft == 1
    dc(:) = H*(x(:).*dc(:))./Hs./max(1e-3,x(:));
  elseif ft == 2
    dc(:) = H*(dc(:)./Hs);
    dv(:) = H*(dv(:)./Hs);
  end
  % The end of counting
    if z == 1
        qi = loop;
    end 
  %% OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES
  l1 = 0; l2 = 1e9; move = 0.2;
  while (l2-l1)/(l1+l2) > 1e-3
    lmid = 0.5*(l2+l1);
    xnew = max(0,max(x-move,min(1,min(x+move,x.*sqrt(-dc./dv/lmid)))));
    if ft == 1
      xPhys = xnew;
    elseif ft == 2
      xPhys(:) = (H*xnew(:))./Hs;
    end
    if sum(xPhys(:)) > volfrac*nelx*nely, l1 = lmid; else l2 = lmid; end
  end
    if loop>10;
        change=abs(sum(c(loop-9:loop-5))-sum(c(loop-4:loop)))/sum(c(loop-4:loop));
    end 
% Update Na
    if loop <= qi
        Na = Na + max(0, xnew(:) - x(:));
    end
    x = xnew;
  %% PRINT RESULTS
  fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,c(loop), ...
    mean(xPhys(:)),change);
  %% PLOT DENSITIES
  figure(z);
  colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;
end
if z == 1
    [dc0,index] = sort(reshape(-dc,1,nelx*nely));
    index0 = index;
    index1 = [];    
end
%% evaluation
C(z)=c(loop);ITER(z)=loop;
C_difference(z)=C(z)/C(1)-1;
if size(xsubopt,1) == 0
    x=reshape(x,1,nelx*nely);
    xsubopt=[xsubopt;x]; nsubopt=[nsubopt;z];
else
    if C_difference(z) < 1.0
        x=reshape(x,1,nelx*nely);
        for j=1:size(xsubopt,1)
            IOU(z)=max(IOU(z),sum(min(x,xsubopt(j,:)))/sum(max(x,xsubopt(j,:))));
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
