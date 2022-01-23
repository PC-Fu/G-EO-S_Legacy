clear
clc

au=load('umatrix00.dat');
au=spconvert(au);

bu=load('umatrix01.dat');
bu=spconvert(bu);

cu=load('umatrix10.dat');
cu=spconvert(cu);

du=load('umatrix11.dat');
du=spconvert(du);

mu = [ au bu; cu du];

b0 = load('urhs0.dat');
b1 = load('urhs1.dat');

bu = [b0'; b1'];

xu = mu \ bu;


%scale = diag(1./diag(mu));
scale = sum(abs(mu),2);
scale = diag(1./scale);
mss = scale*mu;

cond_mss = condest(mss)



a=load('matrix00.dat');
a=spconvert(a);

b=load('matrix01.dat');
b=spconvert(b);

c=load('matrix10.dat');
c=spconvert(c);

d=load('matrix11.dat');
d=spconvert(d);

ms = [ a b; c d];

b0 = load('rhs0.dat');
b1 = load('rhs1.dat');

bs = [b0'; b1'];

xs = ms \ bs;

norm_diff = norm(xu-xs)

cond_mu = condest(mu)
cond_ms = condest(ms)

conditioning_improvement = condest(mu)/condest(ms)

 y0 = load('sol0.dat');
 y1 = load('sol1.dat');
 
 y = [y0'; y1'];
 
 difference = norm(xs-y)
 
 geos_resid = norm(mu*y-bu)
 matlab_scaled_resid = norm(mu*xs-bu)
 matlab_unscaled_resid = norm(mu*xu-bu)
 
% schur=load('schur_est.dat');
% schur=spconvert(schur);
% schur=full(schur)
% 
% ainv = 1./diag(a);
% ainv = diag(ainv);
% ainv = sparse(ainv);
% %ainv = inv(a)
% 
% full(d-b'*ainv*b)
% 
% perturbation_only = full(b'*ainv*b)
% 
% g = xs(1:length(a));
% f = xs(length(a)+1:end);
% 
% norm(g)
% norm(f)
% norm(xu)

