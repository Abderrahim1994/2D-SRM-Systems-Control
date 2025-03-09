close all;
clear all;
clc;
L=1;
C=0.1;
R=10;
G=0.1;
M=[1/(1+L), 0, 0, 0;0, 1, 0, 0;0, 0, 1/(1+C), 0;0, 0, 0, 1];
N=[1, -1, 0, 0;0, 0, 1, -1;0, 0, 1, 1/C;1, 1/L, 0, 0];
En=[1, 0, 0, L;0, 0, 0, 0;0, 1, C, 0;0, 0, 0, 0];
An=[0, -R, 0, 0;1, 0, -1, 0;0, 0, -G, 0;0, 1, 0, -1];
Bh=[1,0;1,0];
Bv=[0,1;0,1];
Bn=[Bh,Bh;Bv,Bv];
E=M*En*N
A=M*An*N
B=M*Bn
        
Dx=0.1;
Dt=0.1;
A11=Dx*A(1:2,1:2)+eye(2,2);
A12=Dx*A(1:2,3:4);
A21=Dt*A(3:4,1:2);
A22=Dt*A(3:4,3:4)+eye(2,2);

X=100;
T=100;
%-----------------------------------------------------------
clear xh xv
xh(:,:,1)=zeros(X,T);
xh(:,:,2)=zeros(X,T);
xv(:,:,1)=zeros(X,T);
xv(:,:,2)=zeros(X,T);
%initiales conditions
xh0=[230;10];
xv0=[230;10];

for tt=1:T
xh(1,tt,1) = xh0(1)*sin(tt);
xh(1,tt,2) = xh0(2)*cos(tt);
end

kV = 2.3; % Rate of change of voltage along x (V/km) 
kI = 0.1; % Rate of change of current along x (A/km)

for xx=1:X
xv(xx,1,1) = xv0(1)*sin(-kV*xx);
xv(xx,1,2) = xv0(2)*cos(-kI*xx);
end
%---------------------------------------------------

for i=1:1:X-1
    for j=1:1:T-1
        xh(i+1,j,1) = A11(1,1)*xh(i,j,1)+A11(1,2)*xh(i,j,2)  +   A12(1,1)*xv(i,j,1)+A12(1,2)*xv(i,j,2) - (L*xv(i,j+1,2))/(L+1);
        xh(i,j,2) = -(1/A11(2,2))*(A11(2,1)*xh(i,j,1) + A12(2,1)*xv(i,j,1) + A12(2,2)*xv(i,j,2));
        xv(i,j+1,1) = A21(1,1)*xh(i,j,1) + A21(1,2)*xh(i,j,2) + A22(1,1)*xv(i,j,1) + A22(1,2)*xv(i,j,2) - xh(i+1,j,2)/(C+1);
        xv(i,j,2) = -(1/A22(2,2))*(A22(2,1)*xv(i,j,1)  + A21(2,1)*xh(i,j,1) + A21(2,2)*xh(i,j,2));
end
end
x=[0:2:(X-1)*2];
t=[0:0.01:(T-1)*0.01];
figure(1);
% subplot(2,2,1);
mesh(x,t,xh(1:X,1:T,1)');
title('x_h_1(x, t)');
xlabel('x(km)');
ylabel('t(s)');

% subplot(2,2,2);
figure(2);
mesh(x,t,xh(1:X,1:T,2)');
title('x_h_2(x, t)');
xlabel('x(km)');
ylabel('t(s)');

% subplot(2,2,3);
figure(3);
mesh(x,t,xv(1:X,1:T,1)');
title('x_v_1(x, t)');
xlabel('x(km)');
ylabel('t(s)');

% subplot(2,2,4);
figure(4);
mesh(x,t,xv(1:X,1:T,2)');
title('x_v_2(x, t)');
xlabel('x(km)');
ylabel('t(s)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms s s1 s2
I=[s1*eye(2), zeros(2);zeros(2), s2*eye(2)];
p11=det(En*I-An);
p1 = vpa(p11,6)
[coeffs_p1, ~] = coeffs(p1, [s1, s2]);

% Vérifier si tous les coefficients ont le même signe
if all(coeffs_p1 > 0) || all(coeffs_p1 < 0)
       fprintf('Tous les coefficients du polynôme p2 sont de même signe');
       else
fprintf('The 2D polynomial is unstable') 
end
       
X=[];
x_h1=[];  
x_h2=[];  
x_v1=[];  
x_v2=[]; 
i=sqrt(-1);
M=10; 
n1=1;
for n=-M:0.1:M
  y_h=n;
  y_v=n;
    s_h=i*y_h;
    s_v=i*y_v;
%---------------------------------------------------------------------------------------------------------------
%----------------------------------------------------------

syms s1 s2

 B=[ coeffs_p1(4)  , coeffs_p1(3)  , coeffs_p1(2) ;
               0   ,           0   ,           0  ; 
     coeffs_p1(1)  ,           0   ,           0 ]; 

% --------------------------------------
aa=B*[1 s_v s_v^2]';
A0=aa(1);
A1=aa(2);
A2=aa(3);
P1=[A2 A1 A0];         
R1=roots(P1)
REsh=real(R1);
resh1=REsh(1);
resh2=REsh(2);
IMsh=imag(R1);
imsh1=IMsh(1);
imsh2=IMsh(2);
% -----------------------------------------
bb=[1 s_h s_h^2]*B;
B0=bb(1);
B1=bb(2);
B2=bb(3);
P2=[B2 B1 B0];         
R2=roots(P2)
REsv=real(R2);
resv1=REsv(1);
resv2=REsv(2);
IMsv=imag(R2);
imsv1=IMsv(1);
imsv2=IMsv(2);
% % ----------------------------------------------------
figure(5);
% subplot(121);
plot(resh1,y_v,'ro')
hold on
plot(resh2,y_v,'b*')
% title('Roots of polynom P(s_h(i*Im(s_v)))');
ylabel('Imaginary axis');
xlabel('Real axis');
legend('Re(s_h_1) ','Re(s_h_2) ');
grid;
% axis([-(M+5) (M+5) -8 1])
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(6);
% subplot(122);
plot(resv1,y_h,'ro')
hold on
plot(resv2,y_h,'b*')
% title('Roots of polynom P(s_v(i*Im(s_h)))');
ylabel('Imaginary axis');
xlabel('Real axis');
legend('Re(s_v_1) ','Re(s_v_2) ');
grid
% axis([-(M+5) (M+5) -18 1])
hold on

if resh1>=0 || resh2>=0 || resv1>=0 || resv2>=0 
    n2=1;
    if resh1>=0
        X= [X,10]; 
        x_h1=[x_h1,resh1];
    end

    if resh2>=0
        X= [X,10];
        x_h2=[x_h2,resh2];
    end
    
    if resv1>=0
        X= [X,10];
        x_v1=[x_v1,resv1];
    end
    
    if resv2>=0
        X= [X,10];
        x_v2=[x_v2,resv2];
    end

n2=n2+1;    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
n1=n1+1;
end
x_h1=double(x_h1);
x_h2=double(x_h2);
x_v1=double(x_v1);
x_v2=double(x_v2);
Y1 = ['s_h1max=',num2str(x_h1)];
Y2=['s_h2max=',num2str(x_h2)];
Y3=['s_v1max=',num2str(x_v1)];
Y4= ['s_v2max=',num2str(x_v2)];

disp(Y1)
disp(Y2)
disp(Y3)
disp(Y4)
if isempty(X)
fprintf('The 2D polynomial is stable')
else
fprintf('The 2D polynomial is unstable')  
end

 