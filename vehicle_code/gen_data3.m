
% boxplot to compare CMDS and WMDS vs k
clear;clc;close all;tic;
total=10000;   % 100/3.29s in CORE i7 8th Gen

cmds=zeros(total,5);wmds=zeros(total,5);gps=zeros(total,5);
CMDS=zeros(1,total);GPS=zeros(1,total);
WMDS1=zeros(1,total);WMDS2=zeros(1,total);WMDS3=zeros(1,total);WMDS4=zeros(1,total);WMDS5=zeros(1,total);
for ii=1:total

N=30;            % node
eta=2;           % dimension
dr=0.1;          % radar range error
du=1;            % UWB range error
r=100;           % radar range
xi=1;            % threshold
dg=1;            % GPS error
n=15;            % GPS-available vehicles
iter=10;      % maximun iteration
X=rand(30,2)*200-100;                       % node distribution [-100£¬100]
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

phi=zeros(1,100);
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

end
%Z=Z(:,:,k+1);           % optimized relative coordinates

Y1=change(G,Y,n);ZZ1=zeros(N,eta);ZZ2=zeros(N,eta);ZZ3=zeros(N,eta);ZZ4=zeros(N,eta);ZZ5=zeros(N,eta);

ZZ1=Z(:,:,1+1);ZZ2=Z(:,:,2+1);ZZ3=Z(:,:,3+1);ZZ4=Z(:,:,4+1);ZZ5=Z(:,:,5+1);
Z1=change(G,ZZ1,n);Z2=change(G,ZZ2,n);Z3=change(G,ZZ3,n);Z4=change(G,ZZ4,n);Z5=change(G,ZZ5,n);

ye=0;ge=0;ze1=0;ze2=0;ze3=0;ze4=0;ze5=0;
for i=1:N
    ye=ye+sqrt((X(i,1)-Y1(i,1))^2+(X(i,2)-Y1(i,2))^2);
    ze1=ze1+sqrt((X(i,1)-Z1(i,1))^2+(X(i,2)-Z1(i,2))^2);
    ze2=ze2+sqrt((X(i,1)-Z2(i,1))^2+(X(i,2)-Z2(i,2))^2);
    ze3=ze3+sqrt((X(i,1)-Z3(i,1))^2+(X(i,2)-Z3(i,2))^2);
    ze4=ze4+sqrt((X(i,1)-Z4(i,1))^2+(X(i,2)-Z4(i,2))^2);
    ze5=ze5+sqrt((X(i,1)-Z5(i,1))^2+(X(i,2)-Z5(i,2))^2);
    ge=ge+sqrt((X(i,1)-G(i,1))^2+(X(i,2)-G(i,2))^2);
end

CMDS(ii)=ye/N;
WMDS1(ii)=ze1/N;WMDS2(ii)=ze2/N;WMDS3(ii)=ze3/N;WMDS4(ii)=ze4/N;WMDS5(ii)=ze5/N;
GPS(ii)=ge/N;

end
CMDS=CMDS';GPS=GPS';WMDS1=WMDS1';WMDS2=WMDS2';WMDS3=WMDS3';WMDS4=WMDS4';WMDS5=WMDS5';
cmds(:,1:5)=[CMDS,CMDS,CMDS,CMDS,CMDS];wmds(:,1:5)=[WMDS1,WMDS2,WMDS3,WMDS4,WMDS5];gps(:,1:5)=[GPS,GPS,GPS,GPS,GPS];
%boxplot([GPS,CMDS,WMDS],'Labels',{'GPS','CMDS','WMDS'},'Whisker',2);
save data3;
toc;
