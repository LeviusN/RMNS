clear all, clc;
format long
%% Zadane hodnoty
xi=[0,2.7,6.8,9,9.2,8.6,8,7.8,7.8,7.9,8];
xi=xi.*2.5;
n=3;
N=7;
iteracie=N-n;
%% Vypocet matice H a Y

% Matica H
H=zeros(N-n,n);
for i=1:n
    for j=1:N-n
        H(j,i)=xi(j+1+(i-1))-xi(j);
    end
end

% Vektor y
y=zeros(N-n,1);
for i=1:N-n
    y(i)=-(xi(i+1)-xi(i));
end

% H=[-14.16 1; 
%     -12.6 1; 
%     -7.08 1; 
%     -4.92 1];
% y=[-2.04 0 4.08 7.92];
% y=y';
% n=2
% iteracie=4;
%% Program RMNS

% Inicializacia
P=10e10*eye(n);
Theta=(zeros(1,n))';
Q=0;
% Algoritmus
for i=1:iteracie
    i
    d=P*H(i,:)';
    ro=pinv(1+H(i,:)*d);
    e=y(i)-H(i,:)*Theta;
    Q=Q+ro*e^2;
    Theta=Theta+ro*e*d
    P=(eye(n)-ro*d*H(i,:))*P;
end


%% Otestovanie
theta=pinv(H'*H)*H'*y;
Y_test=H*Theta;
y_test=H*theta;
% for i=1:3
%     d=P*H(i+1,:)';
%     ro=pinv(1+H(i+1,:)*d);
%     e=y(i+1)-H(i+1,:)*Theta;
%     Q=Q+ro*e^2;
%     Theta=Theta+ro*e*d;
%     P=(eye(n)-ro*d*H(i+1,:))*P;
% end
% 
% Theta

%% Program RMNS
format long
% Inicializacia
n=3;
G=10e10*eye(4);
ThetaG=(zeros(1,n))';
Z=[H,y];
ThetaM=zeros(iteracie, n)
QM=zeros(iteracie,1)
GM=zeros(iteracie*4,4)
ss=0;
l=1;
g=zeros(4,4);
% Algoritmus
gama=G(n,n);
g=G(n,1:n-1);
Theta=-g/gama
Q=1/gama^2;
for k=1:iteracie 
    q=1;
    f=G*Z(k,:)';
    s(q)=1;
    s(q+1)=sqrt(s(q)^2+f(1)^2);
    s(q+2)=sqrt(s(q+1)^2+f(2)^2);
    s(q+3)=sqrt(s(q+2)^2+f(3)^2);
    s(q+4)=sqrt(s(q+3)^2+f(4)^2);
%s=[1 1.0291 1.6718 3.9994]
% Matica G pomocou jednotlivich g(ij)
    for i=1:4 % riadok
        for j=1:4 % stlpec
            if j>=i
                part=0; % nic sa neprida
            else
                part=0;
                for m=j:i-1
                    part=part+f(m)*G(m,j);        
                end
            end          
            g(i,j)=(s(i)/s(i+1))*(G(i,j)-f(i)/(s(i)^2)*part);
        end
    end
G=g
gama=G(4,4);
g=G(4,1:4-1);
Theta=-g/gama
Q=1/gama^2

ThetaM(k,:)=Theta;
QM(k)=Q;
% l=l+3*ss;
% u=k*3;
% GM(l:u,:)=G;
ss=1;
end

