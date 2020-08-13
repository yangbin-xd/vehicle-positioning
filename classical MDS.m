
% classical MDS positionig error
clear;clc;close all;

N=10;             % node
eta=2;            % dimension
de=1;             % range error

X=rand(N,eta)*100-50;                       % node distribution [-50£¬50]
noise=normrnd(0,de/3,N,N);                  % range error noise

d=zeros(N);                                 % Euclid distance matrix
for i=1:N
    for j=1:N
        d(i,j)=sqrt((X(i,1)-X(j,1))^2+(X(i,2)-X(j,2))^2);
    end
end

dn=d+noise;                                 % distance measurement with noise
dn=dn-diag(diag(dn));                       % main diagonal set zero
dn=1/2*(dn+dn');                            % average distance measurement

J=eye(N)-1/N*ones(N);                       % centering matrix
B=-1/2*J*(d.^2)*J;                          % double-centering
[V,D]=eigs(B,eta,'la');                     % EVD
Y=V(:,1:eta)*D(1:eta,1:eta).^(1/2);         % relative coordinates

% coordinate transformation from Y(relative coordinates) to Z(estimated coordinates)
Z=zeros(N,eta);

Xa=X-(sum(X)'/N*ones(1,N))';                % off-mean
Ya=Y-(sum(Y)'/N*ones(1,N))';                % off-mean
P=Xa'*Ya;
[U,S,V1]=svd(P);                            % SVD
R=U*V1';                                    % rotation matrix
Y1=R*Y';Y1=Y1';
t=sum(X)/N-sum(Y1)/N;                       % translation vector
for i=1:N
    Z(i,:)=(R*Y(i,:)'+t')';
end                                         % estimated coordinates

scatter(X(:,1),X(:,2),'bo');hold on; 
scatter(Z(:,1),Z(:,2),'r+');
num=[];                                     % ID set
for num1=1:N
    num=[num;num1];
end
num=num2cell(num);
text(X(:,1)+1.5,X(:,2)+0.2,num);
legend('actual','estlmate');
xlabel('x(m)');ylabel('y(m)');
axis([-60,60,-60,60]);
set(gca,'XTick',-60:20:60);
set(gca,'YTick',-60:20:60);

