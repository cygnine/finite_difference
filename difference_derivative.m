function[d] = difference_derivative(x,y,k,varargin)
% [D] = DIFFERENCE_DERIVATIVE(X,Y,K,{R=0,PERIODIC=false})
%
%   Computes the K'th order finite difference derivative approximation
%   on the unstructured nodal inputs (X,Y). The output satisfies size(D) =
%   size(X).
%
%   The two optional inputs R and PERIODIC are fed right into difference_stencil
%   and have the same meaning as in that function. 

global handles;
newton = handles.speclab.NewtonPolynomials;
fd = handles.finite_difference;

% Create stencil
n = length(x);
[stencil] = fd.difference_stencil(n,k,varargin{:});

% Use stencil to compute interpolants
dd = newton.divided_difference(x(stencil).',y(stencil.'));

% Differentiate and evaluate the interpolants
d = newton.newton_derivative_evaluate(x(stencil).',dd).';
