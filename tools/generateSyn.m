function [ x, rows, cols ] = generateSyn( rows, cols, rnz, cnz )

X = randn(rows, cols)*10;
X = bsxfun(@minus, X, mean(X, 2));

nzcol = (rand(size(X)) > 1 - cnz);
X = X .* nzcol;

nzrow = randperm(rows);
nzrow = nzrow(1:floor((1 - rnz)*rows));
X(nzrow, :) = 0;

x = reshape(X, rows*cols, 1);

end

