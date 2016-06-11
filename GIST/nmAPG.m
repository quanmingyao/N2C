function [ x1, out ] = nmAPG( y, W, lambda, theta,  mu, gamma, para )

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

x0 = zeros(size(W,1), 1);
x1 = x0;
obj = zeros(maxIter, 1);
Time = zeros(maxIter, 1);

t0 = 1;
t1 = 1;

obj(1) = getObject(y, W, x1, lambda, theta,  mu, gamma, para);
tt = tic;
for i = 1:maxIter
    wht = (t0 - 1)/t1;
    
    xi = x1 + wht*(x1 - x0);
    
    zi = W*(W'*xi + para.bias) - Wy;
    zi = xi - zi/sz;
    [xi, objProx] = proxGIST(zi, lambda/sz, theta,  mu/sz, gamma, para, x1);
    
    obji = getObject(y, W, xi, lambda, theta,  mu, gamma, para);
    if(obji > obj(i)*1.05)
        zi = W*(W'*x1 + para.bias) - Wy;
        zi = x1 - zi/sz;
        
        [xi, objProx] = proxGIST(zi, lambda/sz, theta,  mu/sz, gamma, para, x1);
    end
    
    ti = t1;
    t1 = (1 + sqrt(1 + 4*t0^2))/2;
    t0 = ti;
    
    obj(i + 1) = obji; x0 = x1; x1 = xi;
    
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