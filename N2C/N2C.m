function [ x1, out ] = N2C( y, W, lambda, theta,  mu, gamma, para )

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

c = 1;
x0 = zeros(size(W,1), 1);
x1 = x0;
obj = zeros(maxIter, 1);
Time = zeros(maxIter, 1);

obj(1) = getObject(y, W, x1, lambda, theta,  mu, gamma, para);
tt = tic;
for i = 1:maxIter
    wht = (c  - 1)/(c + 2);
    
    xi = x1 + wht*(x1 - x0);
    
    p = para.dcGrad(xi, lambda, theta, mu, gamma, para);
    zi = W*(W'*xi + para.bias) - Wy - p;
    zi = xi - zi/sz;
    zi = proxL1(zi, lambda/(sz*theta));
    zi = reshape(zi, para.row, para.col);
    xi = proxL2(zi, mu/(sz*gamma));
    xi = xi(:);
    
    obji = getObject(y, W, xi, lambda, theta,  mu, gamma, para);
    
    if(i < 5)
        nmobj = inf;
    else
        nmobj = max(obj(i - 4:i));
    end
    if(obji < nmobj)
        c = c + 1;
    else
        c = 1;
    end
    
    obj(i + 1) = obji; x0 = x1; x1 = xi;
    
    delta = abs(obji - obj(i))/obj(i);
    fprintf('iter:%d, obj:(%.2d,%.2d) \n', i, obji, delta);
    
    Time(i) = toc(tt);
    if(delta < para.tol)
        break;
    end
end

out.obj = obj(1:i);
out.Time = Time(1:i);

end

