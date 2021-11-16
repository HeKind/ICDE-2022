clear all
load leaves

    X{1} = double(fea{1}');
    Y{1} = gnd;

    X{2} = double(fea{2}');
    Y{2} = gnd;


for i = 1:length(X)
    c{i} = length(unique(Y{i}));
    cc(i) = length(unique(Y{i}));
end



% Use PCA to make the dimensions of each task equal
options.ReducedDim = 0;
for i = 1:length(X)
X0 = X{i}';
mX0 = mean(X0);
X1 = X0 - ones(size(X{i},2),1)*mX0;
scal = 1./sqrt(sum(X1.*X1)+eps);
scalMat = sparse(diag(scal));
X{i} = X1*scalMat;
[eigvector{i}, eigvalue] = PCA(X{i}, options);
s = sum(eigvalue,1);
e = eigvalue(1);
j=1;
while(e/s <=0.9)
    j= j+1;
    e = e + eigvalue(j);
end
rdimlist(i) = j;
end

dim = max(rdimlist);
for i = 1:length(X)
    dimlist(i) = size(eigvector{i},2);
end
%dim = min(dimlist);
ReducedDim = dim;
for i =1:length(X)
%size(eigvector)
newfea = X{i}*eigvector{i}(:,1:dim);
newfea = NormalizeFea(newfea);
X{i} = newfea';
St = X{i}*X{i}';
invSt{i} = inv(St);
end

%dim =j;
%dim = ReducedDim;%size(X{1},1);
Q = eye(dim);
p = 1;
gamma = 1;
d = max(cc);
% k is the number of neighors.
k = 3
[clusternum,y] = LLMC(X,Q,invSt,gamma,d,c,k,p);

result_label = y{1};
disp('Clustering results of Task 1.')
results = ClusteringMeasure(Y{1}, result_label)

result_label = y{2};
disp('Clustering results of Task 2.')
results = ClusteringMeasure(Y{2}, result_label)


