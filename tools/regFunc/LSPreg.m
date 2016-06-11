function [ z ] = LSPreg( x, lambda, theta )

z = lambda*log(1 + abs(x)/theta);

end

