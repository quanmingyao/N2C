function [ x, out ] = GIST( y, W, lambda, theta,  mu, gamma, para )

if (isfield(para, 'maxIter'))
    maxIter = para.maxIter;
else
    maxIter = 1000;
end

if (~isfield(para, 'tol'))
    para.tol = 1e-4;
end

sz = initStepSize( W );
Wy = W*y;

x = zeros(size(W,1), 1);
obj = zeros(maxIter, 1);
Time = zeros(maxIter, 1);

obj(1) = getObject(y, W, x, lambda, theta,  mu, gamma, para);
tt = tic;
for i = 1:maxIter    
    z = W*(W'*x + para.bias) - Wy;
    z = x - z/sz;
    [x, objProx] = proxGIST(z, lambda/sz, theta,  mu/sz, gamma, para, x);
    
    obji = getObject(y, W, x, lambda, theta,  mu, gamma, para);
    obj(i + 1) = obji;
    
    delta = abs(obj(i) - obji)/obj(i);
    fprintf('iter:%d, obj:(%.2d,%.2d), prox:(%d) \n', ...
        i, obji, delta, length(objProx));
    
    Time(i) = toc(tt);
    if(delta < para.tol)
        break;
    end
end

out.obj = obj(1:i);
out.Time = Time(1:i);

end



