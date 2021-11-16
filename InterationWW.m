function [Q] = InterationQ(X,L,W,p,gamma,d,invSt)
%L: laplacian matrix
%X: data matrix(dim*num)
%gamma: coefficient of L21
%dim: dimension of X
%m: projection dimension of W

INTER_W = 100;
for j = 1:INTER_W
WW =[];
for i =1:length(X)
    WW = [WW W{i}];
end
tempQ = 0.5*p * (sqrt(sum(WW.^2,2)+eps)).^(p-2);
Q = diag(tempQ);

for i = 1:length(X)
    [W{i}] = InterationW(L{i},X{i},gamma,d,Q,invSt{i});
end

t = 0;
for i = 1:length(X)
    t = t + trace(W{i}'*X{i}*L{j}*X{i}'*W{i});
end
w1(j) = t;%w1(j) + trace(W{i}'*X{i}*L{j}*X{i}'*W{i}); % log Tr(WXLXW)
w2(j) = gamma*(sum((sum(WW.^2,2))^(p/2)))^(1/p);% gama*||W||_21
    %w2(i) = gamma*sum(sqrt(sum(W.^2,2)));% gama*||W||_21
WResult(j) = w1(j)+w2(j);

if j > 1 && abs(WResult(j-1)-WResult(j)) < 0.000001
    break;
end;   
end;

%INTER_W = 100;
%Q = eye(dim);
%xlx= X*L*X';
%xlx = (xlx+xlx')/2;
%p=1; % L_2p
%for i = 1:INTER_W
    %tempXLXQ = invSt*(xlx+gamma*Q);
    %[vec,val] = eig(tempXLXQ);
    %[~,di] = sort(diag(val));
    %W = vec(:,di(1:m));
    %W = W*diag(1./sqrt(diag(W'*W)));
    %tempQ = 0.5*p * (sqrt(sum(W.^2,2)+eps)).^(p-2);
    %Q = diag(tempQ);

    %w1(i) = trace(W'*X*L*X'*W); % log Tr(WXLXW)
    %w2(i) = gamma*sum(sqrt(sum(W.^2,2)));% gama*||W||_21
    %WResult(i) = w1(i)+w2(i);

    %if i > 1 && abs(WResult(i-1)-WResult(i)) < 0.000001
    %    break;
    %end;
    
%end;
% WResult = WResult';
end