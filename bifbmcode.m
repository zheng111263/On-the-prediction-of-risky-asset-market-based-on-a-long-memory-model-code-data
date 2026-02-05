%%%%bifbm monte carlo%%%%
H=0.8;K=0.6; 
T=1;
n=64; 
n=2*n;
u=0.5;sigma=0.2;

BB=[];
UU=[];
CC=[];
HH=[];
for k=1:1000   

t=T/n:T/n:T;
S=zeros(n,n); 
for i=1:n
    for j=1:n
        S(i,j)=(1/(2^K))*((t(i)^(2*H)+t(j)^(2*H))^K-abs(t(i)-t(j))^(2*H*K));
    end
end

A=chol(S); 
A=A';

Z=zeros(n,1); 
for i=1:n
    Z(i)=randn;
end
W1=A*Z;

W=W1;
W=W';


t=T/n:T/n:T;

h=T/n;  

X1=u*t-0.5*sigma^2*t.^(2*H*K)+sigma*W;

y1=0;
for i=2:2:n-2 
y1=y1+(abs(X1(i+2)-X1(i)))^(2);
end

S11=zeros(n,n);

for i=1:n
    for j=1:n
        S11(i,j)=(t(i)^(2*H)+t(j)^(2*H))^K-(abs(t(i)-t(j)))^(2*H*K);
    end
end
S1=zeros(64,64);

for i=1:64
    for j=1:64
        S1(i,j)=S11(2*i,2*j);
    end
end

M1=y1;
 
M2=0;
for i=1:n-1 
M2=M2+(abs(X1(i+1)-X1(i)))^2;
end
b=0.5*(1/log(2)*log(M1/M2)+1);

t=T/64:T/64:T;

h=T/64; 
E=1; 
i=1:64;j=i';
A=i'*ones(1,64); B=A'; C=(max(i,j)-min(i,j)).^(2*b);
X=X1(2:2:128);
option=optimset('Display','iter');
x=fminbnd(@(x) funHKK4(x,X,A,B,C,y1,b,t,h,E),b,0.99,option); 
K1=x;
mu=((t/(0.5*h*y1*((A.^(2*b/K1)+B.^(2*b/K1)).^K1-C)))*X.'+0.5*y1/(2^(1-K1)*h^(2*b-1)*E)*(t/(0.5*h*y1*((A.^(2*b/K1)+B.^(2*b/K1)).^K1-C)))*(t.^(2*b)).')/((t/(0.5*h*y1*((A.^(2*b/K1)+B.^(2*b/K1)).^K1-C)))*t.');
sig=(y1/(2^(1-K1)*h^(2*b-1)*E))^(0.5);
mu;
H11=b/K1;
HH=[HH H11];

BB=[BB sig];
UU=[UU mu];
CC=[CC K1];
end
si=mean(BB)
u1=mean(UU)
k1=mean(CC)
h1=mean(HH)
ms=median(BB)
mmu=median(UU)
mk=median(CC)
mh=median(HH)
ssi=std(BB)
su=std(UU)
sk=std(CC)
sh=std(HH)


%%%method B%%%
H=0.6;K=0.7; 
T=1;
n=64; 
u=0.5;sigma=0.2;

BB=[];
UU=[];
tic
for k=1: 1000 
%     

t=T/n:T/n:T;
S=zeros(n,n); 
for i=1:n
    for j=1:n
        S(i,j)=(1/(2^K))*((t(i)^(2*H)+t(j)^(2*H))^K-abs(t(i)-t(j))^(2*H*K));
    end
end

A=chol(S); 
A=A';

Z=zeros(n,1); 
for i=1:n
    Z(i)=randn;
end
W1=A*Z;

W=W1;
W=W';

t=T/n:T/n:T;

h=T/n;  

X=u*t-0.5*sigma^2*t.^(2*H*K)+sigma*W;

x0=[1;1];
S1=zeros(n,n);

for i=1:n
    for j=1:n
        S1(i,j)=(t(i)^(2*H)+t(j)^(2*H))^K-(abs(t(i)-t(j)))^(2*H*K);
    end
end
x=fsolve(@(x) fun7(x,X,S1,H,K,t,h),x0);
mu=x(1);
sig=x(2)^(0.5);
BB=[BB sig];
UU=[UU mu];
end
toc
si=mean(BB)
u1=mean(UU)
ms=median(BB)
mmu=median(UU)
ssi=std(BB)
su=std(UU)



%%%method C%%%
H=0.6;K=0.7; 
T=1;
n=64; 
u=0.5;sigma=0.2;

BB=[];
UU=[];
CC=[];
HH=[];
tic
for k=1:1000  

t=T/n:T/n:T;
S=zeros(n,n); 
for i=1:n
    for j=1:n
        S(i,j)=(1/(2^K))*((t(i)^(2*H)+t(j)^(2*H))^K-abs(t(i)-t(j))^(2*H*K));
    end
end

A=chol(S); 
A=A';

Z=zeros(n,1); 
for i=1:n
    Z(i)=randn;
end
W1=A*Z;

W=W1;
W=W';


t=T/n:T/n:T;

h=T/n;  

X1=u*t-0.5*sigma^2*t.^(2*H*K)+sigma*W;
for i=2:n
    ZZ(i)=X1(i)-X1(i-1);
end
for i=2:n
    tau(i)=(i/n)^(2*H*K)-((i-1)/n)^(2*H*K);
end
xH=tau*n^(2*H*K);
SIG0=S;
SIG1=inv(SIG0)*(eye(n)-(ones(n,1)*ones(1,n)/SIG0)/((ones(1,n)/SIG0)*ones(n,1)));
hatsig2=2*(ZZ*SIG1)*ZZ'/(sqrt(n^2+xH*SIG1*xH'*ZZ*SIG1*ZZ')+n);
hatmu=(1/((ones(1,n)/SIG0)*ones(n,1)))*((ones(1,n)/SIG0)*ZZ'+(0.5*hatsig2*ones(1,n)/SIG0)*xH');
sig=sqrt(hatsig2);
BB=[BB sig];
UU=[UU hatmu];
end
toc
si=mean(BB)
u1=mean(UU)

ms=median(BB)
mmu=median(UU)

ssi=std(BB)
su=std(UU)



%%%method D%%%
H=0.8;K=0.6; 
T=1;
n=64;
u=0.5;sigma=0.2;

BB=[];
UU=[];
CC=[];
HH=[];
tic
for k=1:1000 

t=T/n:T/n:T;
S=zeros(n,n); 
for i=1:n
    for j=1:n
        S(i,j)=(1/(2^K))*((t(i)^(2*H)+t(j)^(2*H))^K-abs(t(i)-t(j))^(2*H*K));
    end
end

A=chol(S);
A=A';

Z=zeros(n,1); 
for i=1:n
    Z(i)=randn;
end
W1=A*Z;

W=W1;
W=W';


t=T/n:T/n:T;

h=T/n;  

X1=u*t-0.5*sigma^2*t.^(2*H*K)+sigma*W;
for i=2:n
    ZZ(i)=X1(i)-X1(i-1);
end
a=0;b=0;c=0;
for i=1:n
    a=a+(i*h)^2;
end
for i=1:n
    b=(i*h)^(2*H*K+1)+b;
end
for i=1:n
    c=(i*h)^(4*H*K)+c;
end
hat1=0;hat2=0;
for i=1:n
    hat1=hat1+(i*h)^(2*H*K)*X1(i);
end
for i=1:n
    hat2=hat2+(i*h)*X1(i);
end
hatmu=(b*hat1-c*hat2)/(b^2-a*c);
hatsig2=2*(a*hat1-b*hat2)/(b^2-a*c);
if ~(hatsig2>=0)
        continue;
end
 
sig=sqrt(hatsig2);
BB=[BB sig];
UU=[UU hatmu];
end
toc
si=mean(BB)
u1=mean(UU)
ms=median(BB)
mmu=median(UU)
ssi=std(BB)
su=std(UU)
toc




%%%%time-varying sigma%%%%
T = 5;
n = 350; 
t_full = T/n:T/n:T;

file_path1 = fullfile('SPXdata350.xlsx');
data1 = readtable(file_path1,'Sheet','Sheet1');
data_vector1 = table2array(data1(:, 1));
X10 = data_vector1'; 

X1 = X10;
X2 = X1;
M1 = 0;
for i = 3:2:n-3 
    if i+2 > n, break; end 
    M1 = M1 + (abs(X2(i+2) - X2(i)))^2;
end
y1 = M1;
M2 = 0;
for i = 1:n-1 
    M2 = M2 + (abs(X1(i+1) - X1(i)))^2;
end
b = 0.5*(1/log(2)*log(M1/M2) + 1);

HH1_lower = b + 1e-6; 
HH1_upper = 0.99 - 1e-6;
KK1_lower = max(b + 1e-6, b / HH1_upper); 
KK1_upper = min(0.99 - 1e-6, b / HH1_lower); 

n2 = n/2; 
h2 = T/n2;
t = T/n2:T/n2:T; 
X = X1(2:2:n); 

g = @(c0, c1) calculate_g(c0, c1, b, h2, T, y1, HH1_lower, HH1_upper);
f = @(c0, c1) b / g(c0, c1); 

initial_params = [0.04, 6e-4, -8e-7]; 

nonlcon = @(x) constraint_fun(x, f, g, b, HH1_lower, HH1_upper, KK1_lower, KK1_upper);

options = optimoptions('fmincon', ...
    'Display', 'iter', ...
    'MaxFunctionEvaluations', 2e5, ...
    'MaxIterations', 8000, ...
    'OptimalityTolerance', 1e-8, ... 
    'FunctionTolerance', 1e-8, ...
    'ConstraintTolerance', 1e-7, ... 
    'Algorithm', 'sqp', ...
    'MaxSQPIter', 1000, ...
    'SubproblemAlgorithm', 'cg');

lb = [0.01, 2e-4, -1.5e-6];
ub = [0.1, 1.5e-3, 1.5e-6];  

[optimal_params, ~, exitflag, output] = fmincon(...
    @(x) nelike4(x, X, t, f, g, b, HH1_lower, HH1_upper, KK1_lower, KK1_upper), ...
    initial_params, ...
    [], [], [], [], ...
    lb, ub, ...
    nonlcon, ...
    options);

mu_opt = optimal_params(1);
a0_opt = optimal_params(2);
a1_opt = optimal_params(3);

KK1 = g(a0_opt, a1_opt);
HH1 = f(a0_opt, a1_opt);


function g_val = calculate_g(c0, c1, b, h2, T, y1, HH1_lower, HH1_upper)
    integrand = @(x) max((c0 + c1*x).^2, 1e-8);
    integral_val = integral(integrand, 0, T, 'RelTol', 1e-4, 'AbsTol', 1e-8);
    log_h2 = log(h2);
    log_integral = log(integral_val);
    log_y1 = log(y1);
    g_raw = 1 + ((2*b - 1)*log_h2 + log_integral - log_y1) / log(2);
    
    HH1_raw = b / g_raw;
    HH1_clamped = max(min(HH1_raw, HH1_upper), HH1_lower);
    g_val = b / HH1_clamped;
    
    KK1_lower = max(b + 1e-6, b / HH1_upper);
    KK1_upper = min(0.99 - 1e-6, b / HH1_lower);
    g_val = max(min(g_val, KK1_upper), KK1_lower);
end

function neg_like = nelike4(params, X, t, f, g, b, HH1_lower, HH1_upper, KK1_lower, KK1_upper)
    mu = params(1);
    a0 = params(2);
    a1 = params(3);
    
    sigma2 = a0 + a1*t + 1e-8;
    sigma2(sigma2 <= 1e-8) = 1e-8;
    log_like = sum( -0.5*log(sigma2) - 0.5*(X - mu*t).^2 ./ sigma2 );
    
    g_val = g(a0, a1); 
    HH1_val = f(a0, a1);
    penalty = 0;
    
    if HH1_val <= HH1_lower || HH1_val >= HH1_upper
        penalty = 1e13 * (HH1_val - (HH1_lower+HH1_upper)/2)^2;
    end
    
    if g_val <= KK1_lower || g_val >= KK1_upper  
        penalty = penalty + 1e10 * (g_val - (KK1_lower+KK1_upper)/2)^2;
    end
    
    neg_like = -log_like + penalty;
end

function [c, ceq] = constraint_fun(x, f, g, b, HH1_lower, HH1_upper, KK1_lower, KK1_upper)
    a0 = x(2);
    a1 = x(3);
    
    g_val = g(a0, a1); 
    HH1_val = f(a0, a1); 
    c(1) = HH1_val - HH1_upper; 
    c(2) = HH1_lower - HH1_val; 
    
    c(3) = g_val - KK1_upper;  
    c(4) = KK1_lower - g_val;  
    
    ceq(1) = HH1_val * g_val - b;
end`



%%%%path simulation%%%%
%%%bifbmSPX500
T=5;
n=350;
h=T/n;
t=h:h:h*n;
n1=n;
file_path1 = fullfile('SPXdata350.xlsx');
data1=readtable(file_path1,'Sheet','Sheet1');
data_vector1=(table2array(data1(:, 1)));
X2=data_vector1';
X2=X2(1:end);
X0=0.0395088133088193;
h1=h;
t1=h:h:n1*h;
M2=0;
for i=1:n1-1 
M2=M2+(abs(X2(i+1)-X2(i)))^2;
end
M1=0;
for i=3:2:n1-3 
M1=M1+(abs(X2(i+2)-X2(i)))^(2);
end
E=1; 
y1=M2;
b=0.5*(1/log(2)*log(M1/M2)+1);
X=X2;
i=1:n1;j=i';
A=i'*ones(1,n1); B=A'; C=(max(i,j)-min(i,j)).^(2*b);
t=t1;
option=optimset('Display','iter');
x=fminbnd(@(x) funHKK4(x,X,A,B,C,y1,b,t,h,E),b,0.99,option);
K=x;
sigma=(y1/(2^(1-K))/(h1^(2*b-1)))^(0.5);
u=((t1/(0.5*h1*y1*((A.^(2*b/K)+B.^(2*b/K)).^K-C)))*X.'+0.5*y1/(2^(1-K)*h1^(2*b-1)*E)*(t1/(0.5*h1*y1*((A.^(2*b/K)+B.^(2*b/K)).^K-C)))*(t1.^(2*b)).')/((t1/(0.5*h1*y1*((A.^(2*b/K)+B.^(2*b/K)).^K-C)))*t1.');
H=b/K;
c=(M2/((T/n1)^(2*b-1)))/(sigma)^2;
YY=[];n2=n;
t2=h:h:h*n2;%out of sample forecasting:n2=100(SPX),3(DNMR);X0=X2(end);
for k=1:20
S=zeros(n2,n2); 
for i=1:n2
    for j=1:n2
        S(i,j)=(1/(2^K))*((t2(i)^(2*H)+t2(j)^(2*H))^K-abs(t2(i)-t2(j)) ^(2*H*K));
    end
end

A=chol(S); 
A=A';

Z=zeros(n2,1); 
for i=1:n2
    Z(i)=randn;
end
W1=A*Z;

W=W1;
W=W';
Y=X0+u*t2-0.5*sigma^2*t2.^(2*H*K)+sigma*W;
YY=[YY;Y];
end
Y1=mean(YY,1);
%out of sample forecasting:YYY=Y1;
dis=sqrt(sum((Y1-X2).^2));
plot([1:n2], X2, 'b');
hold on; 
plot([1:n2], Y1, 'r--');
hold off; 
%%%bifbmDNMR
T=223;
n=300;
h=T/n;
t=h:h:h*n;
n1=n;
file_path1 = fullfile('DNMRdata300.xlsx');
data1=readtable(file_path1,'Sheet','Sheet1');
data_vector1=(table2array(data1(:, 1)));
X1=(data_vector1(1:n))';
X0=0.0257582327413592;
X2=X1;X3=X1;
h1=h;
t1=h:h:n1*h;
M2=0;
for i=1:n1-1 
M2=M2+(abs(X2(i+1)-X2(i)))^2;
end
M1=0;
for i=3:2:n1-3 
M1=M1+(abs(X2(i+2)-X2(i)))^(2);
end
E=1; 
y1=M2;
b=0.5*(1/log(2)*log(M1/M2)+1);

X=X2;
i=1:n1;j=i';
A=i'*ones(1,n1); B=A'; C=(max(i,j)-min(i,j)).^(2*b);
t=t1;
option=optimset('Display','iter');
x=fminbnd(@(x) funHKK4(x,X,A,B,C,y1,b,t,h,E),b,0.99,option);
K=x;
sigma=(y1/(2^(1-K))/(h1^(2*b-1)))^(0.5);
u=((t1/(0.5*h1*y1*((A.^(2*b/K)+B.^(2*b/K)).^K-C)))*X.'+0.5*y1/(2^(1-K)*h1^(2*b-1)*E)*(t1/(0.5*h1*y1*((A.^(2*b/K)+B.^(2*b/K)).^K-C)))*(t1.^(2*b)).')/((t1/(0.5*h1*y1*((A.^(2*b/K)+B.^(2*b/K)).^K-C)))*t1.');
H=b/K;
c=(M2/((T/n1)^(2*b-1)))/(sigma)^2;
YY=[];n2=n;n1=n;
t2=h:h:h*n1;
for k=1:100
S=zeros(n2,n2); 
for i=1:n2
    for j=1:n2
        S(i,j)=(1/(2^K))*((t2(i)^(2*H)+t2(j)^(2*H))^K-abs(t2(i)-t2(j)) ^(2*H*K));
    end
end

A=chol(S); 
A=A';

Z=zeros(n2,1); 
for i=1:n2
    Z(i)=randn;
end
W1=A*Z;

W=W1;
W=W';
Y=X0+u*t2-0.5*sigma^2*t2.^(2*H*K)+sigma*W;
YY=[YY;Y];
end
Y1=mean(YY,1);

plot([1:n1], X1, 'b');
hold on;
plot([1:n1], Y1, 'r--');
hold off;
dis=sqrt(sum((Y1-X1).^2));

%%%fbmSPX500
T=5;
n=350;
h=T/n;
t=h:h:h*n;
XX=[];
n1=n;
rrr=0;
file_path1 = fullfile('SPXdata350.xlsx');
data1=readtable(file_path1,'Sheet','Sheet1');
data_vector1=(table2array(data1(:, 1)));
X0=0.0395088133088193;
X=data_vector1';
X2=X;
h1=h;
t1=h:h:n1*h;
M2=0;
for i=1:n1-1 
M2=M2+(abs(X2(i+1)-X2(i)))^2;
end
M1=0;
for i=3:2:n1-3 
M1=M1+(abs(X2(i+2)-X2(i)))^(2);
end
E=1; 
y1=M2;
b=0.5*(1/log(2)*log(M1/M2)+1);HH=b;

sig2=y1/(h^(2*HH-1)*E*T);
sig=sqrt(sig2);
X=X2;
i=1:n1;j=i';
A=i'*ones(1,n1); B=A'; C=(max(i,j)-min(i,j)).^(2*b);
t=t1;
K=1;
sigma=(y1/(2^(1-K))/(h1^(2*b-1)))^(0.5);
u=((t1/(0.5*h1*y1*((A.^(2*b/K)+B.^(2*b/K)).^K-C)))*X.'+0.5*y1/(2^(1-K)*h1^(2*b-1)*E)*(t1/(0.5*h1*y1*((A.^(2*b/K)+B.^(2*b/K)).^K-C)))*(t1.^(2*b)).')/((t1/(0.5*h1*y1*((A.^(2*b/K)+B.^(2*b/K)).^K-C)))*t1.');
H=b/K;
c=(M2/((T/n1)^(2*b-1)))/(sigma)^2;
YY=[];
t2=h:h:h*n2;
for k=1:20
S=zeros(n2,n2); 
for i=1:n2
    for j=1:n2
        S(i,j)=(1/(2^K))*((t2(i)^(2*H)+t2(j)^(2*H))^K-abs(t2(i)-t2(j)) ^(2*H*K));
    end
end

A=chol(S);
A=A';

Z=zeros(n2,1); 
for i=1:n2
    Z(i)=randn;
end
W1=A*Z;

W=W1;
W=W';
Y=X0+u*t2-0.5*sigma^2*t2.^(2*H*K)+sigma*W;
YY=[YY;Y];
end
X2=mean(YY,1);
dis=sqrt(sum((X2-X).^2));
plot([1:n],X2,'r--',[1:n],X,'b')


%%%fbmDNMR
T=223;
n=300;
h=T/n;
t=h:h:h*n;n1=n;n2=n;
XX=[];
file_path1 = fullfile('DNMRdata300.xlsx');
data1=readtable(file_path1,'Sheet','Sheet1');
data_vector1=(table2array(data1(:, 1)));
X=data_vector1';
X0=0.0257582327413592;
h1=h;
t1=h:h:n1*h;
M2=0;
for i=1:n1-1 
M2=M2+(abs(X(i+1)-X(i)))^2;
end
M1=0;
for i=2:2:n1-2 
M1=M1+(abs(X(i+2)-X(i)))^(2);
end
E=1;
y1=M2;
HH=0.5*(1/log(2)*log(M1/M2)+1);
sig2=y1/(h^(2*HH-1)*E*T);
sig=sqrt(sig2);
i=1:n;j=i';
A=i'*ones(1,n); B=A'; C=(max(i,j)-min(i,j));
Ga=0.5*h^(2*HH)*(A.^(2*HH)+B.^(2*HH)-C.^(2*HH));
mu=(sig2*(t/Ga)*(t.^(2*HH)).'+2*(X/Ga)*(t.'))/(2*(t/Ga)*(t.'));
L=-0.5/n-0.5*log(sig^(2*n)*det(Ga))-0.5/(sig^2)*(X-mu*t+0.5*sig^2*t.^(2*HH))/Ga*(X-mu*t+0.5*sig^2*t.^(2*HH)).';
AIC=2*3-2*L;
for k=1:100
P=fun29(HH,n,T);
P=P';
fbm=zeros(1,n);
fbm=P(2:n+1);
X1=X0+mu*t-0.5*sig2*t.^(2*HH)+sig*fbm;
XX=[XX;X1];
end
X2=mean(XX,1);
Dis=sqrt(sum(abs(X-X2)));
plot([1:n],X2,'r--',[1:n],X,'b')


%%%bmSPX500
T=5;
n=350;
h=T/n;
t=h:h:h*n;
n1=n;
file_path1 = fullfile('SPXdata350.xlsx');
data1=readtable(file_path1,'Sheet','Sheet1');
data_vector1=(table2array(data1(:, 1)));
X0=0.0395088133088193;
X=data_vector1';
XX=[];
tt=mean(t);
wm=mean(X);
wd=var(X);
sig=sqrt(wd/tt);
bmu=wm/tt+0.5*sig^2;
for k=1:20
dw=sqrt(h)*randn(1,n);
w=cumsum(dw); 
bm=w;
X1=X0+bmu*t-0.5*(sig)^2*t+sig*bm;
XX=[XX;X1];
end
X2=mean(XX);
dis=sqrt(sum((X2-X).^2));
plot([1:n],X2,'r--',[1:n],X,'b')


%%%bmDNMR
T=223;
n=300;
h=T/n;
t=h:h:h*n;
n1=n;
file_path1 = fullfile('DNMRdata300.xlsx');
data1=readtable(file_path1,'Sheet','Sheet1');
data_vector1=(table2array(data1(:, 1)));
X=data_vector1';
h2=2*h;XX=[];
tt=mean(t);
wm=mean(X);
wd=var(X);
sig=sqrt(wd/tt);
bmu=wm/tt+0.5*sig^2;
X0=0.0257582327413592;
t=h:h:h*n1;n=n1;
for k=1:100
dw=sqrt(h)*randn(1,n);
w=cumsum(dw); 
bm=w;
X1=X0+bmu*t-0.5*(sig)^2*t+sig*bm;
XX=[XX;X1];
end
X2=mean(XX);
dis=sqrt(sum((X2-X).^2));
plot([1:n],X2,'r--',[1:n],X,'b')




%%%%option pricing%%%%
%%%bifbm
T=1;n=22;h=T/n;t=h:h:h*n;n1=n;

file_path1 = fullfile('SPXpricedata-option.xlsx');
data1=readtable(file_path1,'Sheet','Sheet1');
data_vector1=(table2array(data1(:, 1)));
close_prices = data_vector1;
log_returns = log(close_prices(2:end)./ close_prices(1:end-1));
X2=(log_returns(1:end-1))';
X0=data_vector1(end);
U1=[];n2=n;n1=n;h1=h;
t1=h:h:n1*h;
M2=0;
for i=1:n1-1 
M2=M2+(abs(X2(i+1)-X2(i)))^2;
end
M1=0;
for i=2:2:n1-2 
M1=M1+(abs(X2(i+2)-X2(i)))^(2);
end
E=1; y1=M2;b=0.5*(1/log(2)*log(M1/M2)+1);
X=X2;
i=1:n1;j=i';
A=i'*ones(1,n1); B=A'; C=(max(i,j)-min(i,j)).^(2*b);
t=t1;
option=optimset('Display','iter');
x=fminbnd(@(x) funHKK4(x,X,A,B,C,y1,b,t,h,E),b,0.99,option);
K=x;
sigma=sqrt(252)*(y1/(2^(1-K))/(h1^(2*b-1)))^(0.5);
H=b/K;
KKK=[1800,1900,1950,2000,2050,2100,2150,2200,2250,2300,2350,2375,2400,2425,2450,2475,2500,2525,2540,2550];
r=0.025;T1=1;

for i=1:20  
C1(i)=vpa(X0*normcdf((log(X0/KKK(i))+r*(T1)+0.5*sigma^2*(T1^(2*b)))/(sigma*sqrt(T1^(2*b))), 0, 1)-KKK(i)*exp(-r*(T1))*normcdf(((log(X0/KKK(i))+r*(T1)+0.5*sigma^2*(T1^(2*b)))/(sigma*sqrt(T1^(2*b))))-sigma*sqrt(T1^(2*b)),0,1),60);
end

T2=4;

for i=1:20  
C2(i)=vpa(X0*normcdf((log(X0/KKK(i))+r*(T2)+0.5*sigma^2*(T2^(2*b)))/(sigma*sqrt(T2^(2*b))), 0, 1)-KKK(i)*exp(-r*(T2))*normcdf(((log(X0/KKK(i))+r*(T2)+0.5*sigma^2*(T2^(2*b)))/(sigma*sqrt(T2^(2*b))))-sigma*sqrt(T2^(2*b)),0,1),60);
end

realprice=[1108.85,1008.9,958.9,908.9,858.9,808.9,758.9,708.9,658.95,608.95,558.7,533.7,508.75,483.75,459,434,409,383.8,369.05,359.05];
realprice2=[1109,1009,959.05,909.1,859.1,809.1,759.1,709.1,659.1,609.1,559.2,534.2,509.2,484.2,459.2,434.2,409.2,384.2,369.1,359.2];

P1=(C1-realprice)./realprice;
P2=(C2-realprice2)./realprice2;

%%%fbm
T=1;n=22;h=T/n;t=h:h:h*n;n1=n;

file_path1 = fullfile('SPXpricedata-option.xlsx');
data1=readtable(file_path1,'Sheet','Sheet1');
data_vector1=(table2array(data1(:, 1)));
close_prices = data_vector1;
log_returns = log(close_prices(2:end)./ close_prices(1:end-1));
X2=(log_returns(1:end-1))';
X0=data_vector1(end);
U1=[];n2=n;n1=n;h1=h;
t1=h:h:n1*h;
M2=0;
for i=1:n1-1 
M2=M2+(abs(X2(i+1)-X2(i)))^2;
end
M1=0;
for i=2:2:n1-2 
M1=M1+(abs(X2(i+2)-X2(i)))^(2);
end
E=1; y1=M2;b=0.5*(1/log(2)*log(M1/M2)+1);
t=t1;K=1;H=b;
sigma=sqrt(252)*(y1/(2^(1-K))/(h1^(2*b-1)))^(0.5);
KKK=[1800,1900,1950,2000,2050,2100,2150,2200,2250,2300,2350,2375,2400,2425,2450,2475,2500,2525,2540,2550];
r=0.025;T1=1;

for i=1:20  
C1(i)=vpa(X0*normcdf((log(X0/KKK(i))+r*(T1)+0.5*sigma^2*(T1^(2*b)))/(sigma*sqrt(T1^(2*b))), 0, 1)-KKK(i)*exp(-r*(T1))*normcdf(((log(X0/KKK(i))+r*(T1)+0.5*sigma^2*(T1^(2*b)))/(sigma*sqrt(T1^(2*b))))-sigma*sqrt(T1^(2*b)),0,1),60);
end

T2=4;

for i=1:20  
C2(i)=vpa(X0*normcdf((log(X0/KKK(i))+r*(T2)+0.5*sigma^2*(T2^(2*b)))/(sigma*sqrt(T2^(2*b))), 0, 1)-KKK(i)*exp(-r*(T2))*normcdf(((log(X0/KKK(i))+r*(T2)+0.5*sigma^2*(T2^(2*b)))/(sigma*sqrt(T2^(2*b))))-sigma*sqrt(T2^(2*b)),0,1),60);
end

realprice=[1108.85,1008.9,958.9,908.9,858.9,808.9,758.9,708.9,658.95,608.95,558.7,533.7,508.75,483.75,459,434,409,383.8,369.05,359.05];
realprice2=[1109,1009,959.05,909.1,859.1,809.1,759.1,709.1,659.1,609.1,559.2,534.2,509.2,484.2,459.2,434.2,409.2,384.2,369.1,359.2];

P1=(C1-realprice)./realprice;
P2=(C2-realprice2)./realprice2;

%%%bm
T=1;n=22;h=T/n;t=h:h:h*n;n1=n;

file_path1 = fullfile('SPXpricedata-option.xlsx');
data1=readtable(file_path1,'Sheet','Sheet1');
data_vector1=(table2array(data1(:, 1)));
close_prices = data_vector1;
log_returns = log(close_prices(2:end)./ close_prices(1:end-1));
X2=(log_returns(1:end-1))';
X0=data_vector1(end);
U1=[];n2=n;n1=n;h1=h;
t1=h:h:n1*h;
M2=0;
for i=1:n1-1 
M2=M2+(abs(X2(i+1)-X2(i)))^2;
end
E=1; y1=M2;b=1/2;K=1;
sigma=sqrt(252)*(y1/(2^(1-K))/(h1^(2*b-1)))^(0.5);
KKK=[1800,1900,1950,2000,2050,2100,2150,2200,2250,2300,2350,2375,2400,2425,2450,2475,2500,2525,2540,2550];
r=0.025;T1=1;

for i=1:20  
C1(i)=vpa(X0*normcdf((log(X0/KKK(i))+r*(T1)+0.5*sigma^2*(T1^(2*b)))/(sigma*sqrt(T1^(2*b))), 0, 1)-KKK(i)*exp(-r*(T1))*normcdf(((log(X0/KKK(i))+r*(T1)+0.5*sigma^2*(T1^(2*b)))/(sigma*sqrt(T1^(2*b))))-sigma*sqrt(T1^(2*b)),0,1),60);
end

T2=4;

for i=1:20  
C2(i)=vpa(X0*normcdf((log(X0/KKK(i))+r*(T2)+0.5*sigma^2*(T2^(2*b)))/(sigma*sqrt(T2^(2*b))), 0, 1)-KKK(i)*exp(-r*(T2))*normcdf(((log(X0/KKK(i))+r*(T2)+0.5*sigma^2*(T2^(2*b)))/(sigma*sqrt(T2^(2*b))))-sigma*sqrt(T2^(2*b)),0,1),60);
end

realprice=[1108.85,1008.9,958.9,908.9,858.9,808.9,758.9,708.9,658.95,608.95,558.7,533.7,508.75,483.75,459,434,409,383.8,369.05,359.05];
realprice2=[1109,1009,959.05,909.1,859.1,809.1,759.1,709.1,659.1,609.1,559.2,534.2,509.2,484.2,459.2,434.2,409.2,384.2,369.1,359.2];

P1=(C1-realprice)./realprice;
P2=(C2-realprice2)./realprice2;