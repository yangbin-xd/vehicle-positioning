
% coordinate transformation
function [Z] = change(X,Y,p)
% ������Y������X���任���õ��任�������Z
% ���� X:N*eta,Y:N*eta ǰ3���������ڱ任
% ��� Z:N*eta
[m,n]=size(X);
X1=X(1:p,:);Y1=Y(1:p,:);
X_mean=X1-(sum(X1)'/p*ones(1,p))';X_mean=X_mean';
Y_mean=Y1-(sum(Y1)'/p*ones(1,p))';Y_mean=Y_mean';
P=X_mean*Y_mean';
[U,S,V]=svd(P);                    % SVD
R=U*V';
Y2=R*Y1';Y2=Y2';
t=sum(X1)/p-sum(Y2)/p;
Z=zeros(2,m);Y=Y';
for i=1:m
    Z(:,i)=R*Y(:,i)+t';
end
Z=Z';                              % �任����������Z
end

