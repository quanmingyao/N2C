function [ g ] = dcgradLSP( x, lambda, theta, mu, gamma, para )

% for l1 part
s = sign(x);
z1 = abs(x);
z1 = (lambda/theta)*(1 - 1./(1 + z1 ./ theta));
z1 = z1.*s;

% for l2 part
x = reshape(x, para.row, para.col);
s = sqrt(sum(x.^2, 2)) + 1e-10;
z2 = (mu/gamma)*(1 - 1./(1 + s ./ gamma));
s = bsxfun(@rdivide, x, s);
z2 = bsxfun(@times, z2, s);
z2 = z2(:);

g = z1 + z2;

end

