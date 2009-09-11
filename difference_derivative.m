function[d] = difference_derivative(x,y,k,varargin)
% difference_derivative -- finite-difference derivative
%
% [d] = difference_derivative(x,y,k,{r=0,periodic=false})
%
%   Computes the K'th order finite difference derivative approximation
%   on the unstructured nodal inputs (x,y). The output satisfies size(d) =
%   size(x).
%
%   The two optional inputs r and periodic are fed right into difference_stencil
%   and have the same meaning as in that function. 

global handles;
newton = handles.speclab.newton_polynomials;
fd = handles.finite_difference;

% Create stencil
n = length(x);
[stencil] = fd.difference_stencil(n,k,varargin{:});

% Use stencil to compute interpolants
dd = newton.divided_difference(x(stencil).',y(stencil.'));

% Differentiate and evaluate the interpolants
d = newton.newton_derivative_evaluate(x(stencil).',dd).';
