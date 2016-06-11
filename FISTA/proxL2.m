function [ X ] = proxL2(Z, lambda)

a = sqrt(sum(Z.^2, 2)) + 1e-10;
X = bsxfun(@rdivide, Z, a);

a = max(a - lambda, 0);
X = bsxfun(@times, X, a);

end