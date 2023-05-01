%% create a gaussian mixture using gm function
% parameters
mu = [-10;10];
sigma = reshape([1,1],1,1,2);
p = [0.1,0.9];
gm = gmdistribution(mu,sigma,p);
% plot
x = (-50:0.1:50).';
plot(x,pdf(gm,x)),grid

%% create random numbers using randn()
sigma1 = 1;
sigma2 = 10;
mu1 = 8;
mu2 = -8;
mu = [mu1;mu1+mu2];
sigma = reshape([sigma1,sigma1+sigma2],1,1,2);
p = 0.5;
p1 = [1-p,p];
gm = gmdistribution(mu,sigma,p1);
tiledlayout(1,3)
for ii = 1:5000
    X = randn()*sqrt(sigma1) + mu1;
    Y = randn()*sqrt(sigma2) + mu2;
    b = (rand>=1-p);
    Z = X + b*Y;
    l1(ii) = Z;
end
x = (-50:0.1:50).';
nexttile
histogram(l1,'Normalization','pdf')
hold on 
plot(x,pdf(gm,x),"LineWidth",1.5)

sigma1 = 1;
sigma2 = 10;
mu1 = 6;
mu2 = -6;
mu = [mu1;mu1+mu2];
sigma = reshape([sigma1,sigma1+sigma2],1,1,2);
p = 0.5;
p1 = [1-p,p];
gm = gmdistribution(mu,sigma,p1);

for ii = 1:5000
    X = randn()*sqrt(sigma1) + mu1;
    Y = randn()*sqrt(sigma2) + mu2;
    b = (rand>=1-p);
    Z = X + b*Y;
    l1(ii) = Z;
end
nexttile
histogram(l1,'Normalization','pdf')
hold on 
plot(x,pdf(gm,x),"LineWidth",1.5)


sigma1 = 1;
sigma2 = 10;
mu1 = 4;
mu2 = -4;
mu = [mu1;mu1+mu2];
sigma = reshape([sigma1,sigma1+sigma2],1,1,2);
p = 0.5;
p1 = [1-p,p];
gm = gmdistribution(mu,sigma,p1);

for ii = 1:5000
    X = randn()*sqrt(sigma1) + mu1;
    Y = randn()*sqrt(sigma2) + mu2;
    b = (rand>=1-p);
    Z = X + b*Y;
    l1(ii) = Z;
end
nexttile
histogram(l1,'Normalization','pdf')
hold on 
plot(x,pdf(gm,x),"LineWidth",1.5)

%% group 
n = 2;
for ii = 1:1000
    for jj = 1:n
        X = randn()*sqrt(sigma1) + mu1;
        Y = randn()*sqrt(sigma2) + mu2;
        b = (rand>=1-p);
        Z = X + b*Y;
        l2(ii,jj) = Z;
    end
end
mu = mu1 + p*mu2;
sigma = sigma1 + p*sigma2 + p*(1-p)*mu2*mu2;
s = sum(l2,2);
for kk = 1:1000
    u(kk) = (s(kk) - n*mu)/sqrt(n*sigma) ;
end
tiledlayout(3,3)
nexttile
histogram(u,'Normalization','pdf')
n = 3;
for ii = 1:1000
    for jj = 1:n
        X = randn()*sqrt(sigma1) + mu1;
        Y = randn()*sqrt(sigma2) + mu2;
        b = (rand>=1-p);
        Z = X + b*Y;
        l2(ii,jj) = Z;
    end
end
mu = mu1 + p*mu2;
sigma = sigma1 + p*sigma2 + p*(1-p)*mu2*mu2;
s = sum(l2,2);
for kk = 1:1000
    u(kk) = (s(kk) - n*mu)/sqrt(n*sigma) ;
end
nexttile
histogram(u,'Normalization','pdf')
n = 4;
for ii = 1:1000
    for jj = 1:n
        X = randn()*sqrt(sigma1) + mu1;
        Y = randn()*sqrt(sigma2) + mu2;
        b = (rand>=1-p);
        Z = X + b*Y;
        l2(ii,jj) = Z;
    end
end
mu = mu1 + p*mu2;
sigma = sigma1 + p*sigma2 + p*(1-p)*mu2*mu2;
s = sum(l2,2);
for kk = 1:1000
    u(kk) = (s(kk) - n*mu)/sqrt(n*sigma) ;
end
nexttile
histogram(u,'Normalization','pdf')
n = 5;
for ii = 1:1000
    for jj = 1:n
        X = randn()*sqrt(sigma1) + mu1;
        Y = randn()*sqrt(sigma2) + mu2;
        b = (rand>=1-p);
        Z = X + b*Y;
        l2(ii,jj) = Z;
    end
end
mu = mu1 + p*mu2;
sigma = sigma1 + p*sigma2 + p*(1-p)*mu2*mu2;
s = sum(l2,2);
for kk = 1:1000
    u(kk) = (s(kk) - n*mu)/sqrt(n*sigma) ;
end
nexttile
histogram(u,'Normalization','pdf')
n = 10;
for ii = 1:1000
    for jj = 1:n
        X = randn()*sqrt(sigma1) + mu1;
        Y = randn()*sqrt(sigma2) + mu2;
        b = (rand>=1-p);
        Z = X + b*Y;
        l2(ii,jj) = Z;
    end
end
mu = mu1 + p*mu2;
sigma = sigma1 + p*sigma2 + p*(1-p)*mu2*mu2;
s = sum(l2,2);
for kk = 1:1000
    u(kk) = (s(kk) - n*mu)/sqrt(n*sigma) ;
end
nexttile
histogram(u,'Normalization','pdf')
n = 20;
for ii = 1:1000
    for jj = 1:n
        X = randn()*sqrt(sigma1) + mu1;
        Y = randn()*sqrt(sigma2) + mu2;
        b = (rand>=1-p);
        Z = X + b*Y;
        l2(ii,jj) = Z;
    end
end
mu = mu1 + p*mu2;
sigma = sigma1 + p*sigma2 + p*(1-p)*mu2*mu2;
s = sum(l2,2);
for kk = 1:1000
    u(kk) = (s(kk) - n*mu)/sqrt(n*sigma) ;
end
nexttile
histogram(u,'Normalization','pdf')
n = 50;
for ii = 1:1000
    for jj = 1:n
        X = randn()*sqrt(sigma1) + mu1;
        Y = randn()*sqrt(sigma2) + mu2;
        b = (rand>=1-p);
        Z = X + b*Y;
        l2(ii,jj) = Z;
    end
end
mu = mu1 + p*mu2;
sigma = sigma1 + p*sigma2 + p*(1-p)*mu2*mu2;
s = sum(l2,2);
for kk = 1:1000
    u(kk) = (s(kk) - n*mu)/sqrt(n*sigma) ;
end
nexttile
histogram(u,'Normalization','pdf')
n = 100;
for ii = 1:1000
    for jj = 1:n
        X = randn()*sqrt(sigma1) + mu1;
        Y = randn()*sqrt(sigma2) + mu2;
        b = (rand>=1-p);
        Z = X + b*Y;
        l2(ii,jj) = Z;
    end
end
mu = mu1 + p*mu2;
sigma = sigma1 + p*sigma2 + p*(1-p)*mu2*mu2;
s = sum(l2,2);
for kk = 1:1000
    u(kk) = (s(kk) - n*mu)/sqrt(n*sigma) ;
end
nexttile
histogram(u,'Normalization','pdf')

n = 5000;
for ii = 1:1000
    for jj = 1:n
        X = randn()*sqrt(sigma1) + mu1;
        Y = randn()*sqrt(sigma2) + mu2;
        b = (rand>=1-p);
        Z = X + b*Y;
        l2(ii,jj) = Z;
    end
end
mu = mu1 + p*mu2;
sigma = sigma1 + p*sigma2 + p*(1-p)*mu2*mu2;
s = sum(l2,2);
for kk = 1:1000
    u(kk) = (s(kk) - n*mu)/sqrt(n*sigma) ;
end
nexttile
histogram(u,'Normalization','pdf')
%% get figures
clear; clc;
fig(1) = open('Dmu5.fig');
fig(2) = open('Dmu6.fig');

for i = 1:2
    ax(i) = get(fig(i),'CurrentAxes');
end

figure
for iloop = 1:2
    subplot(2,:,iloop)                                                    
    axChildren = get(ax(iloop),'Children');                                     
    copyobj(axChildren, gca);                                              
end


