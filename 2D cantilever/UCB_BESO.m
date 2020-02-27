%%%% A CODE OF BESO CONBINED WITH UCB FOR 2D CANTILEVER BY SUN. H and MA. L %%%%
function [ITER,C,IOU,C_difference,xsubopt,nsubopt,csubopt]=UCB_BESO(nelx,nely,volfrac,er,rmin,E,nu,rho,xsubopt,csubopt);
% INITIALIZE
penal = 3.; first = ones(1,nely*nelx); iter = 0; iter0 = 0;
pic=zeros(100,1); Sa = zeros(1,nely*nelx); Na = ones(1,nely*nelx);
Nsubopt=10; 
C=zeros(Nsubopt,1);ITER=zeros(Nsubopt,1);C_difference=zeros(Nsubopt,1);IOU=zeros(Nsubopt,1);c =zeros(Nsubopt,100);
nsubopt=[];
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
    x(1:nely,1:nelx) = 1.; vol0=1.; i = 0; change = 1.;
    vol=zeros(100,1); rank=zeros(100,nelx*nely); 
% START iTH ITERATION
while (change > 0.001)&&(i<100)
    i = i + 1; 
	vol0 = max(vol0*(1-er),volfrac);
    vol(i) = vol0;
    if i > 1; 
        olddc = dc; 
    end
% FE-ANALYSIS
    [U]=FE(nelx,nely,x,penal,E,nu);
% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    [KE] = lk(E,nu);
    c(z,i) = 0.;
    for ely = 1:nely
        for elx = 1:nelx
            n1 = (nely+1)*(elx-1)+ely;
            n2 = (nely+1)* elx +ely;
            Ue = U([2*n1-1;2*n1; 2*n2-1;2*n2; 2*n2+1;2*n2+2;2*n1+1;2*n1+2],1);
            c(z,i) = c(z,i) + 0.5*x(ely,elx)^penal*Ue'*KE*Ue;
            dc(ely,elx) = 0.5*x(ely,elx)^(penal-1)*Ue'*KE*Ue;
        end
    end
% FILTERING OF SENSITIVITIES
    dc(:) = H*(dc(:)./Hs);
% The end of counting
    if z == 1 && (vol(i)-volfrac)>0.01
        qi=i+5;
    end    
% STABLIZATION OF EVOLUTIONARY PROCESS
    if i > 1; 
        dc = (dc+olddc)/2.; 
    end
% Update Sa, Na
    if i > 1       
        if (vol(i)-volfrac)>0.01
            iter = iter + 1;
            index1 = find(x == 1);
            Na(index1) = Na(index1) + 1;
            first(index1) = 0;
        end
    end
% Calculate Rank
    dc0 = reshape(dc,1,nelx*nely);
    [dc1,index_dc]=sort(dc0);
    P = 1:-1/nelx/nely:1/nelx/nely;    
    rank(i,index_dc) = rank(i,index_dc) + P;    
% BESO DESIGN UPDATE
    xold = x;
    if z > 1 && (vol(i)-volfrac)>0.01
        dc0 = dc;
        index0 = find(first>=0);
        dc0(index0) = dc0(index0) + rho * sqrt(log(iter)/2./Na(index0));
        index1 = find(first==-1);
        dc0(index1) = dc0(index1)+ rho * sqrt(log(iter)/2./min(Na(index0)));
        [x] = ADDDEL(nelx,nely,vol(i),dc0,x);  
    else
        [x] = ADDDEL(nelx,nely,vol(i),dc,x);  
    end
% PRINT RESULTS
    if i>10;
        change=abs(sum(c(z,i-9:i-5))-sum(c(z,i-4:i)))/sum(c(z,i-4:i));
    end
    disp([' It.: ' sprintf('%4i',i)  ' Obj.: ' sprintf('%10.4f',c(z,i)) ' Vol.: ' sprintf('%6.3f',sum(sum(x))/(nelx*nely)) ' ch.: ' sprintf('%6.3f',change )])
% % PLOT DENSITIES
    figure(z);
    colormap(gray); imagesc(-x); axis equal; axis tight; axis off; 
end
if z == 1
    [dc0,index] = sort(reshape(dc,1,nelx*nely));
    first(index((1-volfrac/2)*nelx*nely+1:end))=-1;
end
%% evaluation
C(z)=c(z,i);ITER(z)=i;
C_difference(z)=C(z)/C(1)-1;
if size(xsubopt,1) == 0
    x=reshape(x,1,nelx*nely);
    xsubopt=[xsubopt;x]; nsubopt=[nsubopt;z]; csubopt=[csubopt;c(z,:)];
else
    if C_difference(z) < 1.0
        x=reshape(x,1,nelx*nely);
        for j=1:size(xsubopt,1)
            IOU(z)=max(IOU(z),length(find((x+xsubopt(j,:))>1.5))/length(find((x+xsubopt(j,:))>0.5)));
            if IOU(z)>0.9
                break;
            end
            if j==size(xsubopt,1)
                xsubopt=[xsubopt;x]; nsubopt=[nsubopt;z]; csubopt=[csubopt;c(z,:)];
            end
        end
    end
end
% if length(nsubopt)==10
%     break
% end
end
end

%%%%%%%%%% OPTIMALITY CRITERIA UPDATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x]=ADDDEL(nelx,nely,volfra,dc,x)
l1 = min(min(dc)); l2 = max(max(dc));
while ((l2-l1)/l2 > 1.0e-5)
    th = (l1+l2)/2.0;
    x = max(0.001,sign(dc-th));
    if sum(sum(x))-volfra*(nelx*nely) > 0;
        l1 = th;
    else
        l2 = th;
    end
end
end

%%%%%%%%%% MESH-INDEPENDENCY FILTER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dcf]=check(nelx,nely,rmin,x,dc)
dcf=zeros(nely,nelx);
for i = 1:nelx
    for j = 1:nely
        sum=0.0;
        for k = max(i-floor(rmin),1):min(i+floor(rmin),nelx)
            for l = max(j-floor(rmin),1):min(j+floor(rmin),nely)
                fac = rmin-sqrt((i-k)^2+(j-l)^2);
                sum = sum+max(0,fac);
                dcf(j,i) = dcf(j,i) + max(0,fac)*dc(l,k);
            end
        end
        dcf(j,i) = dcf(j,i)/sum;
    end
end
end

%%%%%%%%%% FE-ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [U]=FE(nelx,nely,x,penal,E,nu)
[KE] = lk(E,nu);
K = sparse(2*(nelx+1)*(nely+1), 2*(nelx+1)*(nely+1));
F = sparse(2*(nely+1)*(nelx+1),1); 
U = zeros(2*(nely+1)*(nelx+1),1);
for elx = 1:nelx
    for ely = 1:nely
        n1 = (nely+1)*(elx-1)+ely;
        n2 = (nely+1)* elx +ely;
        edof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];
        K(edof,edof) = K(edof,edof) + x(ely,elx)^penal*KE;
    end
end
% DEFINE LOADS AND SUPPORTS (Cantilever)
fixeddofs=[1:2*(nely+1)]; %left
F = sparse(2*(nelx+1)*(nely+1)-nely,1,-1000/10,2*(nely+1)*(nelx+1),1); %right middle
alldofs = [1:2*(nely+1)*(nelx+1)];
freedofs = setdiff(alldofs,fixeddofs);
% SOLVING
U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);
U(fixeddofs,:)= 0;
end

%%%%%%%%%% ELEMENT STIFFNESS MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [KE]=lk(E,nu)
%% MATERIAL PROPERTIES
k=[1/2-nu/6 1/8+nu/8 -1/4-nu/12 -1/8+3*nu/8 -1/4+nu/12 -1/8-nu/8 nu/6 1/8-3*nu/8];

KE = E/(1-nu^2)*[ k(1) k(2) k(3) k(4) k(5) k(6) k(7) k(8); ...
k(2) k(1) k(8) k(7) k(6) k(5) k(4) k(3); ...
k(3) k(8) k(1) k(6) k(7) k(4) k(5) k(2); ...
k(4) k(7) k(6) k(1) k(8) k(3) k(2) k(5); ...
k(5) k(6) k(7) k(8) k(1) k(2) k(3) k(4); ...
k(6) k(5) k(4) k(3) k(2) k(1) k(8) k(7); ...
k(7) k(4) k(5) k(2) k(3) k(8) k(1) k(6); ...
k(8) k(3) k(2) k(5) k(4) k(7) k(6) k(1)];
end
