function [z] = proxLSP(z, mu, theta)

z = proximalRegC(z, length(z), mu, theta, 2);

end