function [x, obj] = proxGIST(y, lambda, theta,  mu, gamma, para, x0)

sz = 1.5;
x = x0;

maxIter = 100;
obj = zeros(maxIter, 1);

obj(1) = proxGISTobj(y, x, lambda, theta,  mu, gamma, para);
for i = 1:maxIter
    p = para.dcGrad(x, lambda, theta, mu, gamma, para);
    z = (x - y) - p;
    z = x - z/sz;
    
    z = proxL1(z, lambda/(sz*theta));
    z = reshape(z, para.row, para.col);
    x = proxL2(z, mu/(sz*gamma));
    x = x(:);

    obji = proxGISTobj(y, x, lambda, theta,  mu, gamma, para);
    obj(i + 1) = obji;
    
    delta = (obj(i) - obji)/obji;
%     assert(delta > 0);
    if(i >= 3 && delta < 1e-12)
        break;
    end
end

obj = obj(1:i);

end

%% --------------------------------------------------------------
function [ obj ] = proxGISTobj(y, x, lambda, theta,  mu, gamma, para)

obj = y - x;
obj = (1/2)*sum(obj(:).^2);

regl1 = para.regFunc(x(:), lambda, theta);
regl1 = sum(regl1);

x = reshape(x, para.row, para.col);
regl2 = sqrt(sum(x.^2, 2));
regl2 = para.regFunc(regl2, mu, gamma);
regl2 = sum(regl2);

obj = obj + regl1 + regl2;

end