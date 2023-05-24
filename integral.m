%% the first question
clear;
a = 0.0000000000000001;
b = 1;
%n = 2;
I = 0;
I1 = 0;
I2 = 0;
f = func(a)+func(b);
for kk = 1:10
    n = 2.^kk;
    h = (b-a)/n;
    h2 = h/2;
% composite trapenoid
    for ii = 1:n-1
        I = I + func(a+ii*h);
    end
    T(kk) = I*h + h2*f;
    ET(kk) = abs(T(kk) + 4/9);
% composite simpson

    for jj = 1:2*n-1 
        if mod(jj,2)==0
            I2 = func(a+jj*0.5*h);
        else
            I1 = func(a+jj*0.5*h);
        end
    end
    S(kk) = (4*I1+2*I2+f)*h2/3;
    ES(kk) = abs(S(kk) + 4/9);
end

%% adaptive
clear;
int=0; 
n=1;
a(1) = 0.0000000000000001;
b(1) = 1;
tol(1)=0.0001; 
app(1)=(b(1)-a(1))*(func(a(1))+func(b(1)))/2;
while n>0 % n is current position at end of the list
    c=(a(n)+b(n))/2; 
    oldapp=app(n);
    app(n)=(c-a(n))*(func(c)+func(a(n)))/2;
    app(n+1)=(b(n)-c)*(func(b(n))+func(c))/2;
    
    if abs(oldapp-(app(n)+app(n+1)))>3*tol(n)
        %b(n+1)=b(n); b(n)=c; % set up new intervals
        a(n+1) = a(n); a(n) = c
        %a(n+1)=c;
        b(n+1)=c
        tol(n)=tol(n)/2; 
        tol(n+1)=tol(n);
        n=n+1;
    else
        int=int+app(n)+app(n+1); % success
        n=n-1; 
        break
    end
end

%% romberg
clear;
a = 0.0000000000000001;
b = 1;
n = 10;
h=(b-a)./(2.^(0:n-1));
r(1,1)=(b-a)*(func(a)+func(b))/2;
E(1,1)= abs(r(1,1)+ 4/9);
for j=2:n
    subtotal = 0;
    for i=1:2^(j-2)
        subtotal = subtotal + func(a+(2*i-1)*h(j));
    end
    r(j,1) = r(j-1,1)/2+h(j)*subtotal;
    E(j,1)= abs(r(j,1)+ 4/9);
    for k=2:j
    r(j,k)=(4^(k-1)*r(j,k-1)-r(j-1,k-1))/(4^(k-1)-1);
    E(j,k)= abs(r(j,k)+ 4/9);
    end
end
%% Gauss
clear;
[x,w] = lgwt(4,-1,1);
w = w.';
for ii = 1:4
    y(ii)=func2(x(ii));
end
integ = dot(w,y);
e = abs(integ+4/9);
%% func itself
clear;
x=linspace(0.0000000000000001,1,1000);
for i=1:1000
    y(i)=func(x(i));
end
plot(x,y)



