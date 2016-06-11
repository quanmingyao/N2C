function [ X ] = proxL1(Z, lambda)

X = sign(Z) .* max(abs(Z) - lambda, 0);

end