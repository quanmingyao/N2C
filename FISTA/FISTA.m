function [ x1, out ] = FISTA( y, W, lambda, mu, para )

if (isfield(para, 'maxIter'))
    maxIter = para.maxIter;
else
    maxIter = 1000;
end

if (~isfield(para, 'tol'))
    para.tol = 1e-4;
end


sz = initStepSize( W );
sz = 1.5*sz;
Wy = W*y;

c = 1;
x0 = zeros(size(W,1), 1);
x1 = x0;
obj = zeros(maxIter, 1);
Time = zeros(maxIter, 1);

obj(1) = fistaObjective(y, W, x1, lambda, mu, para);
tt = tic;
for i = 1:maxIter
    wht = (c - 1)/(c + 2);
    
    xi = x1 + wht*(x1 - x0);
    zi = W*(W'*xi + para.bias) - Wy;
    zi = xi - zi/sz;
    zi = proxL1(zi, lambda/sz);
    zi = reshape(zi, para.row, para.col);
    xi = proxL2(zi, mu/sz);
    xi = xi(:);
    
    obji = fistaObjective(y, W, xi, lambda, mu, para);
    if(obji < obj(i))
        c = c + 1;
    else
        c = 1;
    end
    
    obj(i + 1) = obji; x0 = x1; x1 = xi;
    
    Time(i) = toc(tt);
    delta = abs(obji - obj(i))/obj(i);
    if(delta < para.tol)
        break;
    end
    
    fprintf('iter:%d, obj:(%.2d,%.2d) \n', i, obji, delta);
end

out.obj = obj(1:i);
out.Time = Time(1:i);

end

%% --------------------------------------------------------------
function [ obj ] = fistaObjective(y, W, x, lambda, mu, para)

obj = y - (W'*x + para.bias);
obj = (1/2)*sum(obj(:).^2);

x = reshape(x, para.row, para.col);
obj = obj + lambda*sum(abs(x(:)));
obj = obj + mu*sum(sqrt(sum(x.^2, 2)));

end

