close all;
clear all;
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M=[1, 0; 0, 1];
N=[1, 1; 0, -1];
En=[1, 1; 0, 0];
An=[-1, 0; 1, -1];
E=M*En*N
A=M*An*N

%-------------------------------------------------
Y=[0 0;
   0 1];
yy=Y'*Y
%------------------------------------------------
%--------------------------------------
%-------------------------------------------
AAA=E'*Y
% -----------------------------------------------
cpusec_m = clock;
LMIs = set([]);
%-----------
Iq=1;
Xh=sdpvar(1,1,'symmetric'); Xh= eval('Xh');
Xv=sdpvar(1,1,'symmetric'); Xv= eval('Xv');
X=[Xh,zeros(1);
   zeros(1),Xv];
% % %-------------------------------------------------------
LMIs=LMIs + set(X>0);
%--------------------------------------------------
Gh=sdpvar(1,1,'symmetric'); Gh= eval('Gh');
Gv=sdpvar(1,1,'symmetric'); Gv= eval('Gv');
G=[Gh,zeros(1);
   zeros(1),Gv];
% G=sdpvar(2,2,'full'); G= eval('G');
%---------------------------------------------------
%------------------------------------------------------------
H=sdpvar(2,2,'full'); H= eval('H');
%------------------------------------------------------------
J=sdpvar(2,2,'full'); J= eval('J');
LMIs=LMIs + set(J>0);
%----------------------------------------------------
%----------------------------------------------------                            
%------------------------------------------------------                   
TT11 =  A'*H' + H*A  ;
TT12 = (X*E+Y*G)' + A'*J' - H;
TT22 = -J - J';
%-------------------------------------------------------
%--------------------------------------------------------- 
TERM_1=[TT11    TT12 ;
        
        TT12'   TT22]; 
%---------------------------------------------------------------------
LMIs = LMIs + set(TERM_1<0);
%---------------------------------------------------------------------
cpusec_m = etime(clock,cpusec_m);
%----------------------------------------------------------------
T = size(getvariables(LMIs),2);
sol = solvesdp(LMIs,[],sdpsettings('verbose',0,'solver','sedumi'))
cpusec = sol.solvertime;
%-----------------------------------------------------------------
if sol.problem==0 
format shortG
    
X=double(X)
G=double(G)
H=double(H)
J=double(J)   
eJ=eigs(J)
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
Dx=0.1;
Dt=0.1;

A11OL=Dx*A(1,1)+eye(1);
A12OL=Dx*A(1,2);
A21OL=Dt*A(2,1);
A22OL=Dt*A(2,2)+eye(1);

X=100;
T=100;
%-----------------------------------------------------------
% Traces en boucle ouvert

xhOL(:,:,1)=zeros(X,T);
xvOL(:,:,1)=zeros(X,T);

%initiales conditions
xh0=1;
xv0=1;

for tt=1:T
xhOL(1,tt,1) = xh0*sin(tt*Dt);
xhCL(1,tt,1) = xh0*sin(tt*Dt);
end

for xx=1:X
xvOL(xx,1,1) = xv0*(1-exp(-xx*Dx));
xvCL(xx,1,1) = xv0*(1-exp(-xx*Dx));
end


%---------------------------------------------------

for i=1:1:X
    for j=1:1:T
         xhOL(i+1,j,1) = A11OL*xhOL(i,j,1) +  A12OL*xvOL(i,j,1);
         xvOL(i,j,1) = -(1/A22OL)*A21OL*xhOL(i,j,1);
        
end
end
x=[0:Dx:(X-1)*Dx];
t=[0:1:T-1];

figure(1);
% title('Trace etats xh1, xh2, xv1 et xv2 en BO');
 subplot(1,2,1);
mesh(x,t,xhOL(1:X,1:T,1)');
title('$\bar{x}_h(x, t)$', 'Interpreter', 'latex');
%title('x_h(x, t)');
xlabel('x');
ylabel('t');

% figure(2);
subplot(1,2,2);
mesh(x,t,xvOL(1:X,1:T,1)');
title('$\bar{x}_v(x, t)$', 'Interpreter', 'latex');
xlabel('x');
ylabel('t');   
   
end