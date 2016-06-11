function [ obj ] = getObject( y, W, x, lambda, theta, mu, gamma, para )

obj = y - (W'*x + para.bias);
obj = (1/2)*sum(obj(:).^2);

regl1 = para.regFunc(x(:), lambda, theta);
regl1 = sum(regl1);

x = reshape(x, para.row, para.col);
regl2 = sqrt(sum(x.^2, 2));
regl2 = para.regFunc(regl2, mu, gamma);
regl2 = sum(regl2);

obj = obj + regl1 + regl2;

end

