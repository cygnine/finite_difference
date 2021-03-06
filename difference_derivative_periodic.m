function[d] = difference_derivative_periodic(x,y,k,interval,varargin)
% difference_derivative_periodic -- periodic finite-difference derivative
%
% [d] = difference_derivative_periodic(x,y,k,interval,{r=0})
%
%     Computes the k'th order finite difference derivative approximation on the
%     unstructured nodal inputs (x,y). Because periodicity is assumed here, we
%     need the input interval, a 2-vector specifying the interval of
%     approximation. The output satisfies size(d) = size(x). 
%
%     The optional input r is fed right into difference_stencil and has the same
%     meaning as in that function. 

global packages;
newton = packages.speclab.newton_polynomials;
fd = packages.finite_difference;

xmin = interval(1); xmax = interval(2);

% Create stencil
n = length(x);
if isempty(varargin)
  [stencil,stencil_periodicity] = ...
                  fd.difference_stencil(n,k,'periodic',true);
else
  [stencil,stencil_periodicity] = ...
                  fd.difference_stencil(n,k,'r',varargin{1},'periodic',true);
end

% Compute x values
XInput = x(stencil);
inds = stencil_periodicity==1;
% For indices that wrap down to 1:
XInput(inds) = xmax + (XInput(inds) - xmin);

inds = stencil_periodicity==-1;
% For indices that wrap up to n:
XInput(inds) = xmin - (xmax - XInput(inds));

% Use stencil to compute interpolants
dd = newton.divided_difference(XInput.',y(stencil.'));

% Differentiate and evaluate the interpolants
d = newton.newton_derivative_evaluate(XInput.',dd).';
