% compute eigs
clear all

fid = fopen('../data/data2.txt','r');
formatSpec = '%d %d %f\n';
sizeA = [3 Inf];
A = fscanf(fid,formatSpec,sizeA);
fclose(fid);

i = A(1,:);
j = A(2,:);
val = A(3,:);

A2 = sparse(i,j,val);
numEigen = 50;

tic,
[eigvec,eigval] = eigs(A2, numEigen);
toc
