clear
clc

h=1e-3;
d=0.07;
Q=8;
k1=20;

A21=(-1/1.75)*(2-0.1*Q-4*d*Q);
A22=(-1/1.75)*(0.4);
A23=(-1/1.75)*(-0.2);
A24=(-1/1.75)*(-0.1);
A25=(-1/1.75)*(80);

A41=(-1/1.75)*(-0.5+0.2*Q+d*Q);
A42=(-1/1.75)*(-0.1);
A43=(-1/1.75)*(0.4);
A44=(-1/1.75)*(0.2);
A45=(-1/1.75)*(-20);

alpha1=-3;
alpha2=-3;

epsilon=0.1;

c11=-A21/alpha1-A22;
c12=1;
c13=-A23/alpha2-A24;
c14=0;

c21=-A41/alpha1-A42;
c22=0;
c23=-A43/alpha2-A44;
c24=1;


%%% Robustness of control law to system identification error
% load systemidentification.mat
IdentificationError=0;

if IdentificationError==1

A21=A21*0.99;
A22=A22*1.01;
A23=A23*0.99;
A24=A24*1.01;

A41=A41*0.99;
A42=A42*1.01;
A43=A43*0.99;
A44=A44*1.01;

end
%%%

%x10=0.5;x20=0.5;x30=0.5;x40=0.5;
x10=0.1331;x20=0.1265;x30=0.4508;x40=0.1198;
n=5e5;
x=zeros(n,4);
x(1,:)=[x10,x20,x30,x40];

% Airfoil system
f1=@(x,y,z,m) y;
f2=@(x,y,z,m) (A21/alpha1+A22)*y+(A23/alpha2+A24)*m;
f3=@(x,y,z,m) m;
f4=@(x,y,z,m) (A41/alpha1+A42)*y+(A43/alpha2+A44)*m;

% Sliding mode surface
s1=@(x,y,z,m) c11*x+c12*y+c13*z+c14*m;
s2=@(x,y,z,m) c21*x+c22*y+c23*z+c24*m;

% Noise parameters
 % sigmaB1=0.02;
sigmaB2=0.01;
 % sigmaB3=0.02;
sigmaB4=0.01;

alpha1Levy=1.75;
alpha2Levy=1.75;
alpha3Levy=1.75;
alpha4Levy=1.75;
beta1Levy=0;
beta2Levy=0;
beta3Levy=0;
beta4Levy=0;
sigma1Levy=0.0015;  % Levy noise intensity
sigma2Levy=0.0015;
sigma3Levy=0.0015;
sigma4Levy=0.0015;
u1base=[0];
u2base=[0];

% sign for discontinuous control
for k=2:n

    [k,n]

xnow=x(k-1,1);
ynow=x(k-1,2);
znow=x(k-1,3);
mnow=x(k-1,4);

%%% RK4
K1=f1(xnow,ynow,znow,mnow);
K2=f1(xnow+h*K1/2,ynow+h*K1/2,znow+h*K1/2,mnow+h*K1/2);
K3=f1(xnow+h*K2/2,ynow+h*K2/2,znow+h*K2/2,mnow+h*K2/2);
K4=f1(xnow+h*K3,ynow+h*K3,znow+h*K3,mnow+h*K3);
L1=f2(xnow,ynow,znow,mnow);
L2=f2(xnow+h*L1/2,ynow+h*L1/2,znow+h*L1/2,mnow+h*L1/2);
L3=f2(xnow+h*L2/2,ynow+h*L2/2,znow+h*L2/2,mnow+h*L2/2);
L4=f2(xnow+h*L3,ynow+h*L3,znow+h*L3,mnow+h*L3);
M1=f3(xnow,ynow,znow,mnow);
M2=f3(xnow+h*M1/2,ynow+h*M1/2,znow+h*M1/2,mnow+h*M1/2);
M3=f3(xnow+h*M2/2,ynow+h*M2/2,znow+h*M2/2,mnow+h*M2/2);
M4=f3(xnow+h*M3,ynow+h*M3,znow+h*M3,mnow+h*M3);
N1=f4(xnow,ynow,znow,mnow);
N2=f4(xnow+h*N1/2,ynow+h*N1/2,znow+h*N1/2,mnow+h*N1/2);
N3=f4(xnow+h*N2/2,ynow+h*N2/2,znow+h*N2/2,mnow+h*N2/2);
N4=f4(xnow+h*N3,ynow+h*N3,znow+h*N3,mnow+h*N3);

dx=(K1+2*K2+2*K3+K4)*h/6;
dy=(L1+2*L2+2*L3+L4)*h/6;
dz=(M1+2*M2+2*M3+M4)*h/6;
dm=(N1+2*N2+2*N3+N4)*h/6;

% Noise
Bh=sqrt(h)*randn(4,1);

Mx=stblrnd(alpha1Levy,beta1Levy,h^(1/alpha1Levy),0,1,1);
My=stblrnd(alpha2Levy,beta2Levy,h^(1/alpha2Levy),0,1,1);
Mz=stblrnd(alpha3Levy,beta3Levy,h^(1/alpha3Levy),0,1,1);
Mm=stblrnd(alpha4Levy,beta4Levy,h^(1/alpha4Levy),0,1,1);
Levyh=[sigma1Levy*Mx;sigma2Levy*My;sigma3Levy*Mz;sigma4Levy*Mm];

    %x(k,1)=x(k-1,1)+dx+sigmaB1*Bh(1,1);
    x(k,1)=x(k-1,1)+dx;
    u1base(:,end+1)=-h*epsilon*sign(s2(x(k-1,1),x(k-1,2),x(k-1,3),x(k-1,4)))-A25*(x(k-1,1)^3);
    x(k,2)=x(k-1,2)+dy-h*epsilon*sign(s1(x(k-1,1),x(k-1,2),x(k-1,3),x(k-1,4)))+sigmaB2*Bh(2,1)+Levyh(2,1);
    
    %x(k,3)=x(k-1,3)+dz+sigmaB3*Bh(3,1);
    x(k,3)=x(k-1,3)+dz;
    u2base(:,end+1)=-h*epsilon*sign(s2(x(k-1,1),x(k-1,2),x(k-1,3),x(k-1,4)))-A45*(x(k-1,1)^3);
    x(k,4)=x(k-1,4)+dm-h*epsilon*sign(s2(x(k-1,1),x(k-1,2),x(k-1,3),x(k-1,4)))+sigmaB4*Bh(4,1)+Levyh(4,1);
    
end

load PeriodicTrajectory.mat

subplot(3,2,1)
set(gcf,'color','white');
plot(x(:,1),x(:,3),'b','Linewidth',0.5);
hold on
%plot(x((n-1000):n,1),x((n-1000):n,3),'y','Linewidth',1.5);
%hold on
plot(PeriodicTrajectory(:,1),PeriodicTrajectory(:,3),'r','Linewidth',1.5);
hold on
scatter([0],[0],'g','filled');
hold on
scatter([x10],[x30],'g','filled');
xlabel('α');
ylabel('h');

subplot(3,2,2)
set(gcf,'color','white');
plot3(x(:,1),x(:,2),x(:,3),'b','Linewidth',0.5);
hold on
%plot3(x((n-1000):n,1),x((n-1000):n,2),x((n-1000):n,3),'r','Linewidth',2.5);
%hold on
plot3(PeriodicTrajectory(:,1),PeriodicTrajectory(:,2),PeriodicTrajectory(:,3),'r','Linewidth',2.5);
hold on
scatter3([0],[0],[0],'g','filled');
hold on
scatter3([x10],[x20],[x30],'g','filled');
xlabel('α');
ylabel('dα/dt');
zlabel('h');

Xvar=var(x((n-30000):n,1));
Yvar=var(x((n-30000):n,2));
Zvar=var(x((n-30000):n,3));
Mvar=var(x((n-30000):n,4));

Xvarline=[];
Yvarline=[];
Zvarline=[];
Mvarline=[];
tvarline=[];

NormXZ=sqrt(x(:,1).^2+x(:,3).^2);
%NormXZ=NormXZ(200000:500000,1);
NormXYZM=sqrt(x(:,1).^2+x(:,2).^2+x(:,3).^2+x(:,4).^2);
%NormXYZM=NormXYZM(200000:500000,1);

tline=h*linspace(1,size(NormXZ,1),size(NormXZ,1));
LyapunovXZ=log(NormXZ')./tline;
LyapunovXYZM=log(NormXYZM')./tline;

for i=(size(tline,2)-300000):size(tline,2)
    i
    tvarline(:,end+1)=tline(1,i);
    Xvarline(:,end+1)=var(x((size(tline,2)-301000):i,1));
    Yvarline(:,end+1)=var(x((size(tline,2)-301000):i,2));
    Zvarline(:,end+1)=var(x((size(tline,2)-301000):i,3));
    Mvarline(:,end+1)=var(x((size(tline,2)-301000):i,4));

end



% The maximum Lynapnov exponent for the origin without noise is the maximun
% eigvalue
JacobiO=[0,1,0,0;A21,A22,A23,A24;0,0,0,1;A41,A42,A43,A44];
[~,eigenvalue]=eig(JacobiO);
eigenvalue=real(eigenvalue);
Maxeigenvalue=max(max(eigenvalue));

subplot(3,2,3)
set(gcf,'color','white');
%plot(tline(1,1:100000),LyapunovXZ(1,1:100000),'b');
%hold on;
%plot(tline(1,1:100000),0*tline(1,1:100000),'--r');
%ylim([-5,1]);
plot(tvarline,Xvarline,'b');
hold on
plot(tvarline,0*tvarline,'--r');
xlabel('t');
ylabel('Var(α)');

subplot(3,2,4)
set(gcf,'color','white');
%plot(tline(1,1:100000),LyapunovXYZM(1,1:100000),'b');
%hold on;
%plot(tline(1,1:100000),0*tline(1,1:100000),'--r');
%hold on;
%plot(tline(1,1:100000),Maxeigenvalue*ones(1,100000),'--g');
%ylim([-2.5,2.5]);
plot(tvarline,Zvarline,'b');
hold on
plot(tvarline,0*tvarline,'--r');
xlabel('t');
ylabel('Var(H)');

subplot(3,2,5)
set(gcf,'color','white');
plot(tline,u1base,'b');
ylim([-0.0002,0.0002]);
xlabel('t');
ylabel('u1');

subplot(3,2,6)
set(gcf,'color','white');
plot(tline,u2base,'b');
ylim([-0.0002,0.0002]);
xlabel('t');
ylabel('u2');



