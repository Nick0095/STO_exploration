%%%% A CODE OF BESO CONBINED WITH TS(Gaussian) FOR 2D CANTILEVER BY SUN. H and MA. L %%%%
function [ITER,C,IOU,C_difference,xsubopt,nsubopt]=TS_BESO(nelx,nely,volfrac,er,rmin,E,nu,rate,xsubopt);
% INITIALIZE
penal = 3.; first = ones(1,nely*nelx); iter = 0;
pic=zeros(100,1); thetas = []; n_arms = nelx*nely;
Nsubopt=10; 
C=zeros(Nsubopt,1);ITER=zeros(Nsubopt,1);C_difference=zeros(Nsubopt,1);IOU=zeros(Nsubopt,1);
nsubopt=[];
% init_lists
Sa = zeros(n_arms,1);     % cumulative reward of arm a
Na = zeros(n_arms,1);     % number of times a has been pulled
% init_prior
mu = 0 * ones(1,n_arms);
sigma = 1 * ones(1,n_arms);

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

for z = 1 :10
    x(1:nely,1:nelx) = 1.; vol0=1.; i = 0; change = 1.;
    vol=zeros(100,1);
% START iTH ITERATION
while (change > 0.01)&&(i<100)
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
% STABLIZATION OF EVOLUTIONARY PROCESS
    if i > 1; 
        dc = (dc+olddc)/2.; 
    end
% calculate T
    T = floor((1-vol(i))*nelx*nely);
% calculate rank & MAB & thetas
    rank=zeros(1,nely*nelx);
    dc=reshape(dc,1,nelx*nely);
    [dc0,index]=sort(dc);
    dc=reshape(dc,nely,nelx);
    P = 1:-1/nelx/nely:1/nelx/nely;   
    rank(index) = rank(index) + P;
    thetas = [thetas;rank];
    if iter > 49
        thetas(1,:) = [];
    end
    if iter == 0
        MAB_mu = thetas;
        MAB_sigma = ones(1,nelx*nely);
    else
        [MAB_mu MAB_sigma] = normfit(thetas);
        MAB_sigma(MAB_sigma==0)=10^(-10);
        MAB_sigma = rate * MAB_sigma;
    end
% Choose action
    % active
    iter = iter + 1;
    if vol(i) > volfrac + 0.05
         [Sa,Na,mu,sigma,arm_sequence] = TS(Sa,Na,mu,sigma,n_arms,index,thetas,T,MAB_mu, MAB_sigma);
    end
%% DESIGN UPDATE
    x(1:nely,1:nelx) = 1.;
% TS DESIGN UPDATE    
    if vol(i) > volfrac + 0.05
        x(arm_sequence) = 0.001;
% BESO DESIGN UPDATE
    else
        [x] = ADDDEL(nelx,nely,vol(i),dc,x);
    end
% PRINT RESULTS
    if i>10;
        change=abs(sum(c(i-9:i-5))-sum(c(i-4:i)))/sum(c(i-4:i));
    end
    disp([' It.: ' sprintf('%4i',i)  ' Obj.: ' sprintf('%10.4f',c(i)) ' Vol.: ' sprintf('%6.3f',sum(sum(x))/(nelx*nely)) ' ch.: ' sprintf('%6.3f',change )])
% % PLOT DENSITIES
    figure(z)
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
x=reshape(x,nely,nelx);
end
end

%% Thompson Sampling
function [Sa,Na,mu,sigma,arm_sequence] = TS(Sa,Na,mu,sigma,n_arms,index,thetas,T,MAB_mu, MAB_sigma)
threshold = 0.99;
reward = zeros(T,1); % (T,)
arm_sequence = zeros(T,1);
number = 1:n_arms;
Sa0 = Sa; Na0 = Na; mu0 = mu; sigma0 =sigma;
% reserve higher class
Sa(index(end-0.25*n_arms+1:end)) = []; Na(index(end-0.25*n_arms+1:end)) = []; mu(index(end-0.25*n_arms+1:end)) = []; sigma(index(end-0.25*n_arms+1:end)) = []; 
thetas(:,index(end-0.25*n_arms+1:end)) = []; MAB_mu(index(end-0.25*n_arms+1:end)) = []; MAB_sigma(index(end-0.25*n_arms+1:end)) = [];
number(index(end-0.25*n_arms+1:end))=[]; 
for t = 1:T
    thetas = normrnd(mu,sigma);
    thetas(isnan(thetas)) = 1;
    [thetas0, index] = max(thetas);
    arm = index;        
    % update_lists
    [Sa, Na, reward, arm_sequence] = update_lists(t, arm, Sa, Na, reward, arm_sequence, MAB_mu, MAB_sigma);
    % Posterior update (update_approx)
    [mu, sigma] = update_posterior(arm, reward(t), sigma, mu, MAB_sigma, MAB_mu);
    % delete arm
    Sa0(number(arm)) = Sa(arm); Na0(number(arm)) = Na(arm); mu0(number(arm)) = mu(arm); sigma0(number(arm)) =sigma(arm);
    arm_sequence(t) = number(arm);
    Sa(arm) = []; Na(arm) = []; mu(arm) = []; sigma(arm) = []; thetas(:,arm) = []; MAB_mu(arm) = []; MAB_sigma(arm) = [];
    number(arm)=[];   
end
Sa = Sa0; Na = Na0; mu = mu0; sigma =sigma0;
end

function [Sa, Na, reward, arm_sequence] = update_lists(t, arm, Sa, Na, reward, arm_sequence, MAB_mu, MAB_var)
Na(arm) = Na(arm) + 1; arm_sequence(t) = arm; 
new_reward = normrnd(MAB_mu(arm),MAB_var(arm));
reward(t) = new_reward; Sa(arm) = Sa(arm) + new_reward;
end

function [mu, sigma] = update_posterior(arm, r, sigma, mu, MAB_sigma, MAB_mu)
mu(arm) = (MAB_sigma(arm)^2 * mu(arm) + r * sigma(arm) ^ 2) / (MAB_sigma(arm)^2 + sigma(arm) ^ 2);
sigma(arm) = sqrt((MAB_sigma(arm) * sigma(arm)) ^ 2 / (MAB_sigma(arm) ^ 2 + sigma(arm) ^ 2));
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
