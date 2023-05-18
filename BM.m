%% Brownian motion
%% generate
clear;
n = 100;
T = 1000;
N = 1000;
t = T/N;
B = zeros(1,N+1);
B(1) = 0;
for jj = 1:n
    for ii = 1:T
        B(ii+1) = B(ii)+sqrt(t)*randn;
    end
    x = linspace(0,T,N+1);
    plot(x,B);
    hold on;
end
xlabel('t')
ylabel('B(t)')
%% equation
clear;
T = 100;
N = 1000;
t = T/N;
v = 1;
x0 = 0;
X = zeros(1,N+1);
X(1) = x0;
n = 10;
sigma = 0;
alpha = 1;
tiledlayout(2,2)
nexttile
for jj = 1:n
    for ii = 1:N
        dB = sqrt(t)*randn;
        X(ii+1) = (1-alpha*t)*X(ii) + alpha*v*t + sigma*dB;
    end
    x = linspace(0,T,N+1);
    plot(x,X)
    hold on
end
xlabel('t')
ylabel('X(t)')


sigma = 1;
nexttile
for jj = 1:n
    for ii = 1:N
        dB = sqrt(t)*randn;
        X(ii+1) = (1-alpha*t)*X(ii) + alpha*v*t + sigma*dB;
    end
    x = linspace(0,T,N+1);
    plot(x,X)
    hold on
end
xlabel('t')
ylabel('X(t)')
nexttile

sigma = 2;
for jj = 1:n
    for ii = 1:N
        dB = sqrt(t)*randn;
        X(ii+1) = (1-alpha*t)*X(ii) + alpha*v*t + sigma*dB;
    end
    x = linspace(0,T,N+1);
    plot(x,X)
    hold on
end
xlabel('t')
ylabel('X(t)')
nexttile
sigma = 3;
for jj = 1:n
    for ii = 1:N
        dB = sqrt(t)*randn;
        X(ii+1) = (1-alpha*t)*X(ii) + alpha*v*t + sigma*dB;
    end
    x = linspace(0,T,N+1);
    plot(x,X)
    hold on
end
xlabel('t')
ylabel('X(t)')



%% monte carlo method (n = 1000)
clear;
T = 100;
N = 100;
t = T/N;
v = 0;
sigma = 1;
alpha = 1;
x0 = 0;
re = 10;
n = 1000;
mu = zeros(1,n);
s = zeros(1,n);
Y = zeros(1,n);
for jj = 1:1000
    for ii = 1:n
        dB = sqrt(t)*randn;
        Y(ii) = (1-alpha*t)*x0 + alpha*v*t + sigma*dB;
    end
    mu(jj) = mean(Y);
    s(jj) = var(Y);
end
h = histogram(s,'Normalization','probability');
hold on
xvals = (h.BinEdges(2:end)+h.BinEdges(1:end-1))/2;
plot(xvals, h.Values, 'r');

%% n = 100
E2 = zeros(1,re);
D2 = zeros(1,re);
n = 100;
Y = zeros(1,n);
for ii = 1:10
    for jj = 1:n
        dB = sqrt(t)*randn;
        Y(jj) = (1-alpha*t)*x0 + alpha*v*t + sigma*dB;    
    end
    E2(ii) = mean(Y);
    D2(ii) = var(Y);
end
%% n = 10
E1 = zeros(1,re);
D1 = zeros(1,re);
n = 10;
Y = zeros(1,n);
for ii = 1:10
    for jj = 1:n
        dB = sqrt(t)*randn;
        Y(jj) = (1-alpha*t)*x0 + alpha*v*t + sigma*dB;    
    end
    E1(ii) = mean(Y);
    D1(ii) = var(Y);
end
%% n = 10000
E4 = zeros(1,re);
D4 = zeros(1,re);
n = 10;
Y = zeros(1,n);
for ii = 1:10
    for jj = 1:n
        dB = sqrt(t)*randn;
        Y(jj) = (1-alpha*t)*x0 + alpha*v*t + sigma*dB;    
    end
    E4(ii) = mean(Y);
    D4(ii) = var(Y);
end
%% plot
% mean
tiledlayout(1,2)
nexttile
for ii = 1:10
    scatter(10,E1(ii),'filled'); hold on
    scatter(100, E2(ii),'filled'); hold on
    scatter(1000,E3(ii),'filled'); hold on
    
    set(gca,'xscale','log');
    plot([10,100,1000],[mean(E1),mean(E2),mean(E3)],'DisplayName','E(x)');  
end
% variance
nexttile
for jj = 1:10
    scatter(10,D1(jj),'filled'); hold on
    scatter(100,D2(jj),'filled'); hold on
    scatter(1000,D3(jj),'filled'); hold on
    set(gca,'xscale','log');
    plot([10,100,1000],[mean(D1),mean(D2),mean(D3)],'DisplayName','D(X)');  
end
%% the second equation
clear;
T = 100;
N = 1000;
t = T/N;

sigma = 1;
sigma1 = 1;
sigma2 = 1;
alpha = 0;
theta = 1;
x0 = 0;
s0 = 0;
X = zeros(1,N+1);
S = zeros(1,N+1);
X(1) = x0;
S(1) = s0;
v = 0;
tiledlayout(3,2)
nexttile
    for ii = 1:N
        dB = sqrt(t)*randn;
        dW = sqrt(t)*randn;
        X(ii+1) = (1-alpha*t)*X(ii)+alpha*v*t+sigma*dB;
        S(ii+1) = (1-theta*t)*S(ii)+theta*t*X(ii)+sigma1*dB+sigma2*dW;
    end
    x = linspace(0,T,N+1);
    plot(x,X); hold on
    plot(x,S);hold on
    
    alpha = 1;
nexttile
    for ii = 1:N
        dB = sqrt(t)*randn;
        dW = sqrt(t)*randn;
        X(ii+1) = (1-alpha*t)*X(ii)+alpha*v*t+sigma*dB;
        S(ii+1) = (1-theta*t)*S(ii)+theta*t*X(ii)+sigma1*dB+sigma2*dW;
    end
    x = linspace(0,T,N+1);
    plot(x,X); hold on
    plot(x,S);hold on
    %plot(x,v+0*x)
    alpha = 2;
    nexttile
    for ii = 1:N
        dB = sqrt(t)*randn;
        dW = sqrt(t)*randn;
        X(ii+1) = (1-alpha*t)*X(ii)+alpha*v*t+sigma*dB;
        S(ii+1) = (1-theta*t)*S(ii)+theta*t*X(ii)+sigma1*dB+sigma2*dW;
    end
    x = linspace(0,T,N+1);
    plot(x,X); hold on
    plot(x,S);hold on
    %plot(x,v+0*x)
    alpha = 4;
    nexttile
    for ii = 1:N
        dB = sqrt(t)*randn;
        dW = sqrt(t)*randn;
        X(ii+1) = (1-alpha*t)*X(ii)+alpha*v*t+sigma*dB;
        S(ii+1) = (1-theta*t)*S(ii)+theta*t*X(ii)+sigma1*dB+sigma2*dW;
    end
    x = linspace(0,T,N+1);
    plot(x,X); hold on
    plot(x,S); hold on
    %plot(x,v+0*x)
alpha = 10;
    nexttile
    for ii = 1:N
        dB = sqrt(t)*randn;
        dW = sqrt(t)*randn;
        X(ii+1) = (1-alpha*t)*X(ii)+alpha*v*t+sigma*dB;
        S(ii+1) = (1-theta*t)*S(ii)+theta*t*X(ii)+sigma1*dB+sigma2*dW;
    end
    x = linspace(0,T,N+1);
    plot(x,X); hold on
    plot(x,S); hold on
    alpha = 20;
    nexttile
    for ii = 1:N
        dB = sqrt(t)*randn;
        dW = sqrt(t)*randn;
        X(ii+1) = (1-alpha*t)*X(ii)+alpha*v*t+sigma*dB;
        S(ii+1) = (1-theta*t)*S(ii)+theta*t*X(ii)+sigma1*dB+sigma2*dW;
    end
    x = linspace(0,T,N+1);
    plot(x,X); hold on
    plot(x,S); hold on
%% monte caro again
clear;
T = 100;
N = 1000;
t = T/N;
sigma = 1;
sigma1 = 1;
sigma2 = 1;
alpha = 1;
theta = 1;
x0 = 0;
s0 = 0;
X = zeros(1,N+1);
S = zeros(1,N+1);
X(1) = x0;
S(1) = s0;
v = 0;
for jj = 1:1000
    for ii = 1:N
        dB = sqrt(t)*randn;
        dW = sqrt(t)*randn;
        X(ii) = (1-alpha*t)*x0+alpha*v*t+sigma*dB;
        S(ii) = (1-theta*t)*s0+theta*t*X(ii)+sigma1*dB+sigma2*dW;
    end
    mu(jj) = mean(S);
    s(jj) = var(S);
end
h = histogram(s,'Normalization','probability');
hold on
xvals = (h.BinEdges(2:end)+h.BinEdges(1:end-1))/2;
plot(xvals, h.Values, 'r');   