function [ s ] = initStepSize( D )

r = randn(size(D, 2), 1);
s = top1svd( D, r, 5, 1e-6);

s = s^2;

end

