function [ X ] = proxAvg( Y, lambda, theta,  mu, gamma, para )

X1 = para.proxFunc(Y(:), lambda, theta);
X1 = reshape(X1, size(Y,1), size(Y,2));

X2 = sqrt(sum(Y.^2, 2)) + 1e-10;
Y  = bsxfun(@rdivide, Y, X2);
X2 = para.proxFunc(X2, mu, gamma);
X2 = bsxfun(@times, Y, X2);

X = (lambda*X1 + mu*X2)/(lambda + mu);

end
