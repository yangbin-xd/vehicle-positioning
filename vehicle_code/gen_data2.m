
% boxplot to compare CMDS and WMDS vs 3sigma_g
clear;clc;close all;tic;
total=10000;   % 100/18.74s in CORE i7 8th Gen

cmds=zeros(total,5);wmds=zeros(total,5);gps=zeros(total,5);
for dg=0.5:0.5:2.5

CMDS=zeros(1,total);WMDS=zeros(1,total);GPS=zeros(1,total);
for ii=1:total

N=30;            % node
eta=2;           % dimension
dr=0.1;          % radar range error
du=1;            % UWB range error
r=100;           % radar range
xi=1;            % threshold
%dg=1;            % GPS error
n=15;           % GPS-available vehicles
iter=1000;       % maximun iteration

X=rand(N,eta)*200-100;                       % node distribution [-100£¬100]
G=X+normrnd(0,dg/3,N,eta);                   % GPS coordinate

d=zeros(N);                                  % Euclid distance matrix
for i=1:N
    for j=1:N
        d(i,j)=sqrt((X(i,1)-X(j,1))^2+(X(i,2)-X(j,2))^2);
    end
end

dn=zeros(N);                                 % EDM + noise
for i=1:N
    for j=1:N
        if d(i,j)>r
            dn(i,j)=d(i,j)+normrnd(0,du/3,1,1);
        else
            dn(i,j)=d(i,j)+normrnd(0,dr/3,1,1);
        end
    end
end

dn=dn-diag(diag(dn));                       % main diagonal set zero
dn=1/2*(dn+dn');                            % average distance measurement

J=eye(N)-1/N*ones(N);                       % centering matrix
B=-1/2*J*(dn.^2)*J;                         % double-centering
[V,D]=eigs(B,eta,'la');                     % EVD
Y=V(:,1:eta)*D(1:eta,1:eta).^(1/2);         % initial relative coordinates

delta=zeros(N,N,iter);
for i=1:N
    for j=1:N
        delta(i,j,1)=sqrt((Y(i,1)-Y(j,1))^2+(Y(i,2)-Y(j,2))^2);
    end
end

W=zeros(N);
for i=1:N
    for j=1:N
        if d(i,j)>r
            W(i,j)=1/(du/3)^2;
        else
            W(i,j)=1/(dr/3)^2;
        end
    end
end

phi=zeros(1,iter);
for i=1:N
    for j=1:N
        phi(1)=phi(1)+W(i,j)*(delta(i,j)-dn(i,j))^2;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SMACOF

A=zeros(N);                                 % a matrix help to display tr
I=eye(N);                                   % identity matrix
E=zeros(N);
for i=1:N
    for j=i:N
        A=(I(:,i)-I(:,j))*(I(:,i)-I(:,j))';
        E=E+W(i,j)*A;
    end
end

Z=zeros(N,eta,iter);
Z(:,:,1)=Y;

for k=1:iter
    
s=zeros(N);
for i=1:N
    for j=1:N
        s(i,j)=dn(i,j)/(sqrt((Z(i,1,k)-Z(j,1,k))^2+(Z(i,2,k)-Z(j,2,k))^2));
        s(i,i)=0;
    end
end

F=zeros(N);
for i=1:N
    for j=i:N
        A=(I(:,i)-I(:,j))*(I(:,i)-I(:,j))';
        F=F+W(i,j)*s(i,j)*A;
    end
end
        
Z(:,:,k+1)=pinv(E)*F*Z(:,:,k);

for i=1:N
    for j=1:N
        delta(i,j,k+1)=sqrt((Z(i,1,k+1)-Z(j,1,k+1))^2+(Z(i,2,k+1)-Z(j,2,k+1))^2);
    end
end

for i=1:N
    for j=1:N
        phi(k+1)=phi(k+1)+W(i,j)*(delta(i,j,k+1)-dn(i,j))^2;
    end
end
if (phi(k)-phi(k+1))<xi
    break;
end
end
Z=Z(:,:,k+1);           % optimized relative coordinates

Y1=change(G,Y,n);
Z1=change(G,Z,n);

ye=0;ze=0;ge=0;
for i=1:N
    ye=ye+sqrt((X(i,1)-Y1(i,1))^2+(X(i,2)-Y1(i,2))^2);
    ze=ze+sqrt((X(i,1)-Z1(i,1))^2+(X(i,2)-Z1(i,2))^2);
    ge=ge+sqrt((X(i,1)-G(i,1))^2+(X(i,2)-G(i,2))^2);
end

CMDS(ii)=ye/N;
WMDS(ii)=ze/N;
GPS(ii)=ge/N;

end
CMDS=CMDS';WMDS=WMDS';GPS=GPS';
cmds(:,6-dg*2)=CMDS;wmds(:,6-dg*2)=WMDS;gps(:,6-dg*2)=GPS;
%boxplot([GPS,CMDS,WMDS],'Labels',{'GPS','CMDS','WMDS'},'Whisker',2);
end
save data2;
toc;
