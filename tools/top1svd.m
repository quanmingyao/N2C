function [s, q] = top1svd( A, r, maxIter, tol)

if(~exist('tol', 'var'))
    tol = 1e-5;
end

y = A*r;
[q, ~] = qr(y, 0);
err = zeros(maxIter, 1);
for i = 1:maxIter
    y = A'*q;
    y = A *y;
    [iq, ~] = qr(y, 0);
    
    err(i) = norm(iq(:,1) - q(:,1), 1);
    q = iq;
    
    if(err(i) < tol)
        break;
    end
end

s = A'*q;
s = sqrt(sum(s.^2));

end

