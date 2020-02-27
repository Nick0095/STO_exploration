%%%% A CODE OF BESO CONBINED WITH IDS(Gaussian) FOR 2D CANTILEVER BY SUN. H and MA. L %%%%
function [ITER,C,IOU,C_difference,xsubopt,nsubopt]=VIDS_BESO(nelx,nely,volfrac,er,rmin,E,nu,rate,xsubopt);
% INITIALIZE
iter = 0; rank=zeros(1000,nelx*nely); penal = 3.; n_arms = 10;
Nsubopt=10; 
C=zeros(Nsubopt+1,1);ITER=zeros(Nsubopt+1,1);C_difference=zeros(Nsubopt+1,1);IOU=zeros(Nsubopt+1,1);
nsubopt=[];
mu = 0 * zeros(1,n_arms); sigma = 1 * ones(1,n_arms);
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

%% VIDS_sample: using MC sampling for Gaussian Bandit Problems with normal prior
[rank0, qi, x, c]=BESO_Gaussian(nelx,nely,volfrac,er,rmin,E,nu);
iter = iter + qi;
M = qi;
data = rank0(1:M,:); rank0 = data;
[MAB_mu0 MAB_sigma0]=normfit(data);
[MAB_mu1,index]=sort(MAB_mu0);

C(1)=c;
if size(xsubopt,1) == 0
    x=reshape(x,1,nelx*nely);
    xsubopt=[xsubopt;x]; nsubopt=[nsubopt;1];
end

%% theta
M = 10000;

for z = 1:10
%% START iTH ITERATION
N = 80*50/n_arms; i = 0; rem = ones(1,n_arms) * N; rem(1:0.2*n_arms) = 0;
x(1:nely,1:nelx) = 1.; change = 1.;  vol0 = 1; vol = zeros(1,100);
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
%% thetas    
    if vol(i) > volfrac + 0.05
        iter = iter + 1;
        dc0 = reshape(dc,1,nelx*nely);
        % rank
        [dc1,index_dc]=sort(dc0);
        P = 1:-1/nelx/nely:1/nelx/nely;    
        rank(iter,index_dc) = rank(iter,index_dc) + P;
        MAB_mu1 = rank(iter,:);
        % MAB_sigma2
        [MAB_mu2 MAB_sigma2]=normfit(rank(1:iter,:));
        for j = 1:n_arms
            MAB_sigma(j) = rate * mean(MAB_sigma2(index(80*50/n_arms*(j-1)+1:80*50/n_arms*j)));
        end
        thetas = zeros(n_arms,M);
        for arm = 1:n_arms
            MAB_mu(arm) = sum(x(index(80*50/n_arms*(arm-1)+1:80*50/n_arms*arm)).*...
                           MAB_mu1(index(80*50/n_arms*(arm-1)+1:80*50/n_arms*arm)))/...
                           sum(x(index(80*50/n_arms*(arm-1)+1:80*50/n_arms*arm)));
            thetas(arm,:) = normrnd(MAB_mu(arm),MAB_sigma(arm),1,M);
        end
    end 
%% DESIGN UPDATE
    % calculate T
    if i == 1
        T = floor((1-vol(i))*nelx*nely); 
    else
        T = floor((vol(i-1)-vol(i))*nelx*nely); 
    end
    if vol(i) > volfrac + 0.05
        if T > 0
            [mu,sigma,arm_sequence,rem] = VID2(N,n_arms,thetas,M,T,MAB_mu,MAB_sigma,rem,mu,sigma);
        end
        x = ones(nely,nelx);
        for j = 0.2*n_arms + 1:n_arms
            if T > 0
                del = length(find(arm_sequence==j));
            end
            % reverse sampling 
            value = normrnd(MAB_mu1(N*(j-1)+1:N*j),MAB_sigma2(N*(j-1)+1:N*j));
            [value0, index2] = sort(value,'descend');
            x(index(N*(j-1) + index2(1:N-rem(j)))) =0.001;
        end
    else
        [x]=ADDDEL(nelx,nely,vol(i),dc,x,i);
    end
% PRINT RESULTS
    if i>10;
        change=abs(sum(c(i-9:i-5))-sum(c(i-4:i)))/sum(c(i-4:i));
    end
    disp([' It.: ' sprintf('%4i',i)  ' Obj.: ' sprintf('%10.4f',c(i)) ' Vol.: ' sprintf('%6.3f',sum(sum(x))/(nelx*nely)) ' ch.: ' sprintf('%6.3f',change )])
% PLOT DENSITIES
    figure(z+1);
    colormap(gray); imagesc(-x); axis equal; axis tight; axis off; 
end
%% evaluation
C(z+1)=c(i);ITER(z+1)=i;
C_difference(z+1)=C(z+1)/C(1)-1;
if size(xsubopt,1) == 0
    x=reshape(x,1,nelx*nely);
    xsubopt=[xsubopt;x]; nsubopt=[nsubopt;z+1];
else
    if C_difference(z+1) < 1.0
        x=reshape(x,1,nelx*nely);
        for j=1:size(xsubopt,1)
            IOU(z+1)=max(IOU(z+1),length(find((x+xsubopt(j,:))>1.5))/length(find((x+xsubopt(j,:))>0.5)));
            if IOU(z+1)>0.9
                break;
            end
            if j==size(xsubopt,1)
                xsubopt=[xsubopt;x]; nsubopt=[nsubopt;z+1];
            end
        end
    end
end
x=reshape(x,nely,nelx);
end
end

function [mu0,sigma0,arm_sequence1,rem] = VID2(N,n_arms,thetas,M0,T,MAB_mu, MAB_sigma,rem,mu,sigma)
rem(1:0.2*n_arms) = 0;
index0 = find(rem==0);
index1 = find(rem>0);
n_arms = length(index1);
mu0 = mu; sigma0 = sigma;
thetas(index0,:) = []; MAB_mu(:,index0) = []; MAB_sigma(:,index0) = [];
threshold = 0.99;
%% init_lists
Sa = zeros(n_arms,1);     % cumulative reward of arm a
Na = zeros(n_arms,1);     % number of times a has been pulled
reward = zeros(T,1); % (T,)
arm_sequence = zeros(T,1);
arm_sequence1 = zeros(T,1);
Maap = zeros(n_arms,n_arms);
p_a = zeros(n_arms,1);

%% main
for t = 1:T
    % computeIDS
    [arm, p_a, optimal_arm] = computeVIDS(n_arms, threshold, Maap, p_a, thetas, M0);
    % update_lists
    [Sa, Na, reward, arm_sequence] = update_lists(t, arm, Sa, Na, reward, arm_sequence, MAB_mu, MAB_sigma);
    arm_sequence1(t) = index1(arm);
    % Posterior update (update_approx)
    [mu, sigma] = update_posterior(arm, reward(t), sigma, mu, MAB_sigma, MAB_mu);
    rem(index1(arm)) = rem(index1(arm)) - 1;
    mu0(index1(arm)) = mu(arm); sigma0(index1(arm)) = sigma(arm);
    if rem(index1(arm)) == 0        
        index1 = find(rem>0);
        n_arms = length(index1);
        Sa(arm) = []; Na(arm) = []; Maap(arm,:) = []; Maap(:,arm) = [];
        p_a(arm) = []; mu(arm) = []; sigma(arm) = []; MAB_mu(:,arm) = []; MAB_sigma(:,arm) = [];
        thetas(arm,:) = [];         
    else
        thetas(arm,:) = normrnd(mu(arm),sigma(arm),1,M0);
    end
end
end

function [mu, sigma] = update_posterior(arm, r, sigma, mu, MAB_sigma, MAB_mu)
mu(arm) = (MAB_sigma(arm)^2 * mu(arm) + r * sigma(arm) ^ 2) / (MAB_sigma(arm)^2 + sigma(arm) ^ 2);
sigma(arm) = sqrt((MAB_sigma(arm) * sigma(arm)) ^ 2 / (MAB_sigma(arm) ^ 2 + sigma(arm) ^ 2));
end

function [arm, p_a, optimal_arm] = computeVIDS(n_arms, threshold, Maap, p_a, thetas, M)
mu = mean(thetas, 2);
[theta0, theta_hat] = max(thetas);
for a = 1:n_arms 
    for ap = 1:n_arms 
        t = thetas(ap, find(theta_hat == a));
        Maap(ap, a) = mean(t);        
        Maap(isnan(Maap)) = 0;
        Maap(isinf(Maap)) = 1e308;
        if ap == a
            p_a(a) = length(t)/M;
        end
    end
end
optimal_arm = 5;
if max(p_a) >= threshold
    [p_a0, index] = max(p_a);
    optimal_arm = index(1);
    arm = optimal_arm;
else
    rho_star = 0;
    for a = 1:n_arms
        rho_star = rho_star + p_a(a) * Maap(a, a);
    end
    delta = rho_star - mu;
    
    v = [];
    for a = 1:n_arms
        v0 = 0;
        for ap = 1:n_arms
            v0 = v0 + p_a(ap) * (Maap(a, ap) - mu(a)) ^ 2;
        end
        v = [v v0];
    end
    % IDSAction
    arm = IDSAction(n_arms, delta, v);
end
end

function arm = IDSAction(n_arms, delta, g)
Q = zeros(n_arms, n_arms);
IR = ones(n_arms, n_arms) * inf;
q = 0:1/999:1;
flag = 0;
for a = 1:n_arms - 1
    for ap = a + 1:n_arms
        if g(a) < 1e-6 || g(ap) < 1e-6
            arm = rd_argmax(-g); flag = 1; break;
        end
        da = delta(a);dap = delta(ap);ga = g(a);gap = g(ap);
        qaap = q(rd_argmax(-(q * da + (1 - q) * dap) .^ 2 ./ (q * ga + (1 - q) * gap + 1e-10)));
        IR(a, ap) = (qaap * (da - dap) + dap) ^ 2 / (qaap * (ga - gap) + gap + 1e-10);
        Q(a, ap) = qaap;
    end
    if flag == 1
        break
    end
end
if flag == 0
    amin = rd_argmax(-reshape(IR, n_arms * n_arms, 1)) - 1;
    ap = floor(amin/n_arms); a = mod(amin, n_arms);
    b = binornd(1,Q(a+1, ap+1));
    arm = floor(b * a + (1 - b) * ap) + 1;
end
end  

function a = rd_argmax(vector)
% Compute random among eligible maximum indices
m = max(vector);
indices = find(vector == m);
a = indices(randperm(length(indices),1));
end

function [Sa, Na, reward, arm_sequence] = update_lists(t, arm, Sa, Na, reward, arm_sequence, MAB_mu, MAB_var)
Na(arm) = Na(arm) + 1; arm_sequence(t) = arm; 
new_reward = normrnd(MAB_mu(arm),MAB_var(arm));
reward(t) = new_reward; Sa(arm) = Sa(arm) + new_reward;
end

%%%%%%%%%% OPTIMALITY CRITERIA UPDATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x]=ADDDEL(nelx,nely,vol,dc,x,i)
l1 = min(min(dc)); l2 = max(max(dc));
while ((l2-l1)/l2 > 1.0e-5)
    th = (l1+l2)/2.0;
    x = max(0.001,sign(dc-th));
    if sum(sum(x))-vol*(nelx*nely) > 0;
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
F = sparse(2*(nelx+1)*(nely+1)-nely,1,-1,2*(nely+1)*(nelx+1),1); %right middle
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
