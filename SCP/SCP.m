function [ x, out ] = SCP( y, W, lambda, theta,  mu, gamma, para )

if (isfield(para, 'maxIter'))
    maxIter = para.maxIter;
else
    maxIter = 1000;
end

if (~isfield(para, 'tol'))
    para.tol = 1e-4;
end

sz = initStepSize( W );
sz = sz*1.15;

Wy = W*y;

x = zeros(size(W,1), 1);
obj = zeros(maxIter, 1);
Time = zeros(maxIter, 1);

obj(1) = getObject(y, W, x, lambda, theta,  mu, gamma, para);
tt = tic;
for i = 1:maxIter
    p = para.dcGrad(x, lambda, theta, mu, gamma, para);
    
    z = W*(W'*x + para.bias) - Wy - p;
    z = x - z/sz;
    z = proxL1(z, lambda/(sz*theta));
    z = reshape(z, para.row, para.col);
    x = proxL2(z, mu/(sz*gamma));
    x = reshape(x, para.row*para.col, 1);
    
    obji = getObject(y, W, x, lambda, theta,  mu, gamma, para);
    obj(i + 1) = obji;
    
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

