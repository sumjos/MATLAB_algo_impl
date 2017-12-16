%This function solves Ax=B for x, where A is a tridiagonal matrix
%represneted by vectors a, b and c, abd B is vector f
function [X] = solveTriadiagonalByGE(a,b,c,f)
for k = 1:length(a)-1
    m = vpa(b(k+1)/a(k));
    a(k+1) = a(k+1) - m*c(k);
    f(k+1) = f(k+1) - m*f(k);
end
X = vpa(zeros(size(a)));
X(end) = vpa(f(end)/a(end));
for r = length(a)-1:-1:1
    X(r) = vpa((f(r) - c(r)*X(r+1))/a(r));
end
end