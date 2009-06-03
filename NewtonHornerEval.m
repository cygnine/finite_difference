% MATLAB File : NewtonHornerEval.m
%
% * Creation Date : 2009-06-03
%
% * Last Modified : Wed 03 Jun 2009 01:50:22 PM EDT
%
% * Created By : Akil Narayan
%
% * Purpose : Uses Horner's method of evaluation to compute the interpolant at
%   the points z given the modal coefficients c corresponding to a Newton
%   polynomial basis at nodes x.
%   length(c) = n
%   length(x) = n-1 (any additional nodal locations are ignored)
%

function[y] = NewtonHornerEval(z,c,x)

n = length(c);
y = c(n)*ones(size(z));
for q = (n-1):(-1):1
  y = y.*(z-x(q)) + c(q);
end
