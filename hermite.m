%-------------------------------------------------------------------------------
% File:         hermite.m
% Authors:      Igor Janjic
% Description:  Implements Hermite interpolation.
%               Computes coefficients of the Hermite polynomial using divided
%               differences.
% Inputs:
%   x           Nodes.
%   fx          Function values in nodes.
%   dfx         function values in nodes.
% Output:
%   a           Coefficients of Hermite polynomial.
%   z           Values of Hermite polynomial.
%-------------------------------------------------------------------------------

function [a, z] = hermite(x, fx, dfx)

% Determine n and reserve space for Q, z and a.
n = length(x) - 1;
Q = zeros(2*(n + 1) + 1);
z = zeros(2*(n + 1) + 1, 1);
a = zeros(2*(n + 1) + 1, 1);
for i = 0:n
    z(2*i + 1) = x(i + 1);
    z(2*i + 2) = x(i + 1);
end

% First two columns of Q and z.
for i = 0:n
    Q(2*i + 1, 1) = fx(i + 1);
    Q(2*i + 2, 1) = fx(i + 1);
    Q(2*i + 2, 2) = dfx(i + 1);
end
for i = 1:n
    Q(2*i + 1, 2) = (Q(2*i + 1, 1) - Q(2*i, 1))/(z(2*i + 1) - z(2*i));
end

% Other columns of Q.
for i = 2:(2*n + 1)
    for j = 2:i
        Q(i + 1, j + 1) = (Q(i + 1, j) - Q(i, j))/(z(i + 1) - z(i - j + 1));
    end
end

% Coefficients a of polynomial.
for i = 0:(2*n + 1)
    a(i + 1) = Q(i + 1, i + 1);
end
end
