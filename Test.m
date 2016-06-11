clear; clc;
close all;

[ x, rows, cols ] = generateSyn( 50, 200, 0.25, 0.75 );
n = 10000;
W = sprandn(rows*cols, n, 0.01);
W = bsxfun(@rdivide, W, sqrt(sum(W.^2, 1)));
y = W'*x;
y = y + randn(size(y))*0.05;

para.row = rows;
para.col = cols;
para.bias = - mean(y);

para.tol = 1e-9;
para.maxIter = 5000;
para.regFunc = @ LSPreg;
para.dcGrad = @ dcgradLSP;

clear rows cols;

%% ---------------------------------------------------------------

lambda = 1;
mu = 0.01;
theta = 1;
gamma = 1;
% ------------------------------------------------------------------------
method = 1;
[ xp{method}, out{method} ] = SCP( y, W, lambda, theta, mu, gamma, para );

% ------------------------------------------------------------------------
method = 2;
[ xp{method}, out{method} ] = GIST(y, W, lambda, theta, mu, gamma, para);

% ------------------------------------------------------------------------
method = 3;
[ xp{method}, out{method} ] = nmAPG(y, W, lambda, theta, mu, gamma, para);

% ------------------------------------------------------------------------
method = 4;
para.proxFunc = @ proxLSP;
[ xp{method}, out{method} ] = GDPAN (y, W, lambda, theta, mu, gamma, para);
para = rmfield(para, 'proxFunc');

% ------------------------------------------------------------------------
method = 5;
[ xp{method}, out{method} ] = N2C(y, W, lambda, theta, mu, gamma, para);

%% ---------------------------------------------------------------
objMin = min(cat(1, out{1}.obj, out{2}.obj, out{3}.obj, out{4}.obj, out{5}.obj));

figure;
semilogy(out{1}.Time, out{1}.obj - objMin);
hold on
semilogy(out{2}.Time, out{2}.obj - objMin);
semilogy(out{3}.Time, out{3}.obj - objMin);
semilogy(out{4}.Time, out{4}.obj - objMin);
semilogy(out{5}.Time, out{5}.obj - objMin);

legend('SCP','GIST','nmAPG','GD-PAN','N2C');
xlabel('cputime (seconds)');
ylabel('objective value - best');
