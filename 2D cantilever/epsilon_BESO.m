%%%% A CODE OF BESO CONBINED WITH EPSILON GREEDY POLICY FOR 2D CANTILEVER BY SUN. H and MA. L %%%%
function [ITER,C,IOU,C_difference,xsubopt,nsubopt]=epsilon_BESO(nelx,nely,volfrac,er,rmin,E,nu,emax,winx,winy,sym,xsubopt);
% INITIALIZE
penal = 3.; unfixvol=1-volfrac/2;
Nsubopt=10; 
C=zeros(Nsubopt,1);ITER=zeros(Nsubopt,1);C_difference=zeros(Nsubopt,1);IOU=zeros(Nsubopt,1);
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
    vol=zeros(100,1);
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
    c(i) = 0.;
    for ely = 1:nely
        for elx = 1:nelx
            n1 = (nely+1)*(elx-1)+ely;
            n2 = (nely+1)* elx +ely;
            Ue = U([2*n1-1;2*n1; 2*n2-1;2*n2; 2*n2+1;2*n2+2;2*n1+1;2*n1+2],1);
            c(i) = c(i) + 0.5*x(ely,elx)^penal*Ue'*KE*Ue;
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
% RL / BESO DESIGN UPDATE
    if z > 1 && (vol(i)-volfrac)>0.01
        e=emax+(0-emax)*min(1,1/(qi-5)*i)^3;
        [x]=RL_ADDDEL(nelx,nely,unfixvol,volnext(i),dc,e,winx,winy,sym);
    else
        [x] = ADDDEL(nelx,nely,vol(i),dc,x);  
    end
% PRINT RESULTS
    if i>10;
        change=abs(sum(c(i-9:i-5))-sum(c(i-4:i)))/sum(c(i-4:i));
    end  
    if z == 1 && (vol(i)-volfrac)>0.01
        volnext(i)=mean(x(:));
    end
    disp([' It.: ' sprintf('%4i',i)  ' Obj.: ' sprintf('%10.4f',c(i)) ' Vol.: ' sprintf('%6.3f',sum(sum(x))/(nelx*nely)) ' ch.: ' sprintf('%6.3f',change )])
% PLOT DENSITIES
    figure(z);
    colormap(gray); imagesc(-x); axis equal; axis tight; axis off; 
end
%% evaluation
C(z)=c(i);ITER(z)=i;
C_difference(z)=C(z)/C(1)-1;
if size(xsubopt,1) == 0
    x=reshape(x,1,nelx*nely);
    xsubopt=[xsubopt;x]; nsubopt=[nsubopt;z];
else
    if C_difference(z) < 1.0
        x=reshape(x,1,nelx*nely);
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

%%%%%%%%%% RL UPDATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x]=RL_ADDDEL(nelx,nely,unfixvol,volnext,dc,e,winx,winy,sym)
% symmetric constraint
if sym==1
    nely=nely/2;
    dc=dc(1:nely,:);
end 
% delete 1-e elements
[dc0,index]=sort(reshape(dc,1,nelx*nely));
base=floor((1-e)*(1-volnext)*(nelx*nely));
x=ones(1,nely*nelx);
dc_critical=dc0(base); % RANDOM CHOICE
index_critical=find(dc0==dc_critical);
base2=index_critical(1);
index(index_critical)= index(index_critical(1)-1+randperm(length(index_critical)));
x(index(1:base))=0.001;
x=reshape(x,nely,nelx);

% delete e elements by window
m=floor(e*(1-volnext)*(nelx*nely)/(2*winx+1)/(2*winy+1));
if m>0
n=unfixvol*nelx*nely-base;
p=randperm(n);
for i=1:m
    ind=index(base+p(i));
    indx=ceil(ind/nely);
    indy=mod(ind,nely);
    if indy==0
        indy=nely;
    end
    x(indy,max(1,indx-winx):min(indx+winx,nelx))=0.001;
    if winy>0
        for j=1:winy
            x(max(1,indy-j),max(1,indx-winx):min(indx+winx,nelx))=0.001;
            x(min(nely,indy+j),max(1,indx-winx):min(indx+winx,nelx))=0.001;
        end
    end
end
end

% delete rest elements
n=floor((mean(x(:))-volnext)*(nelx*nely));
if n>0
    ind2=find(x(index(1:unfixvol*nelx*nely))==1);
    p=randperm(length(ind2));
    x(index(ind2(p(1:n))))=0.001;
end
if sym==1
    x=[x;flipud(x)];
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
