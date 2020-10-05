function [i, j] = distanceVectorIndToSub(m, k)
% Copied from https://stackoverflow.com/a/41269711 by edelburg

% m = number of observations or points

% r(j) is the no. of elements in cols 1..j, belonging to the upper triangular part
r = cumsum(1:m-1);
% p(j) is the no. elements in cols 1..j, belonging to the lower triangular part
p = cumsum(m-1:-1:1);
% The linear index of value D(k)
q = find(p >= k, 1);
% The subscript indices of value D(k)
[i, j] = ind2sub([m m], k + r(q));

end