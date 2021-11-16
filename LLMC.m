function [clusternum,y] = LLMC(X,Q,invSt,gamma,d,c,k,p)
% Input
% X: dim*num data matrix
% gamma: coefficient of L21
% d: projection dim of W(dim*d)
% c: number of clusters
% k: nearest neighobrs

%Output
%id: sorted features by ||w_i||_2
%dim =j;
%dim = ReducedDim;%size(X{1},1);
%Q = eye(dim);
for i = 1:length(X)
num{i} = size(X{i},2);
iter_continue(i) = 1;
%dim{i} = size(X{i},1);
end
for j = 1:length(X)
distX{j} = L2_distance_1(X{j},X{j});
[distX1, idx{j}] = sort(distX{j},2);
A = zeros(num{j});
rr = zeros(num{j},1);
for i = 1:num{j}
    di = distX1(i,2:k+2);
    rr(i) = 0.5*(k*di(k+1)-sum(di(1:k)));
    id = idx{j}(i,2:k+2);
    A(i,id) = (di(k+1)-di)/(k*di(k+1)-sum(di(1:k))+eps);
end;
r = mean(rr);
lambda{j} = r;

A0 = (A+A')/2;
D0 = diag(sum(A0));
L0 = D0 - A0;
L{j} = L0;

[F{j}, ~, ~]=eig1(L0, c{j}, 0);
[W{j}] = InterationW(L0,X{j},gamma,d,Q,invSt{j});
%size(W{j})
end

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
    t = t + trace(W{i}'*X{i}*L{i}*X{i}'*W{i});
end
w1(j) = t;%w1(j) + trace(W{i}'*X{i}*L{j}*X{i}'*W{i}); % log Tr(WXLXW)
w2(j) = gamma*(sum((sum(WW.^2,2)).^(p/2)))^(1/p);% gama*||W||_21
    %w2(i) = gamma*sum(sqrt(sum(W.^2,2)));% gama*||W||_21
WResult(j) = w1(j)+w2(j);

if j > 1 && abs(WResult(j-1)-WResult(j)) < 0.000001
    break;
end;   
end;
NITER = 30;
for iter = 1:NITER
    for j = 1:length(X)
        if iter_continue(j)==1
    distf = L2_distance_1(F{j}',F{j}');
    distx = L2_distance_1(W{j}'*X{j},W{j}'*X{j});
    if iter>5
        [~, idx{j}] = sort(distx,2);
    end;
    A = zeros(num{j});
    for i=1:num{j}
        idxa0 = idx{j}(i,2:k+1);
        dfi = distf(i,idxa0);
        dxi = distx(i,idxa0);
        ad = -(dxi+lambda{j}*dfi)/(2*r);
        A(i,idxa0) = EProjSimplex_new(ad);
    end;
    
    A = (A+A')/2;
    At{j} = A;
    D = diag(sum(A));
    L{j} = D-A;
    
    [W{j}] = InterationW(L{j},X{j},gamma,d,Q,invSt{j});%InterationW(L,X{i},gamma,dim,d);
    
    F_old = F{j};
    [F{j}, ~, ev]=eig1(L{j}, c{j}, 0);
    %evs(:,iter+1) = ev;
    
    fn1 = sum(ev(1:c{j}));
    fn2 = sum(ev(1:c{j}+1));
    if fn1 > 0.000000001
        lambda{j} = 2*lambda{j};
    elseif fn2 < 0.00000000001
        lambda{j} = lambda{j}/2;  F{j} = F_old;
    elseif iter>1
       iter_continue(j) = 0;
    end;
        else
            [W{j}] = InterationW(L{j},X{j},gamma,d,Q,invSt{j});
        end;
    end
    
%update Q
WW =[];
for i =1:length(X)
    WW = [WW W{i}];
end
tempQ = 0.5*p * (sqrt(sum(WW.^2,2)+eps)).^(p-2);
Q = diag(tempQ);

end

for i = 1:length(X)
[clusternum{i}, y{i}]=graphconncomp(sparse(At{i})); 
y{i} = y{i}';
if clusternum{i} ~= c{i}
    sprintf('Can not find the correct cluster number: %d', c{i})
end
end

%sqW = (W.^2);
%sumW = sum(sqW,2);
%[~,id] = sort(sumW,'descend');


