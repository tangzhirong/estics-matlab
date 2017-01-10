S = load('out.txt');
B = load('onezero-Matrix.txt');
r = 30;
lambda = 1e-3;
MaxIter = 50;

[L,R,v] = estics(S,B,r,lambda,MaxIter);
X = L*R';

save result.txt X -ascii