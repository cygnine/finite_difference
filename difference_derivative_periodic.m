function[d] = difference_derivative_periodic(x,y,k,interval,varargin)
% [D] = DIFFERENCE_DERIVATIVE_PERIODIC(X,Y,K,INTERVAL,{R=0})
%
%     Computes the K'th order finite difference derivative approximation on the
%     unstructured nodal inputs (X,Y). Because periodicity is assumed here, we
%     need the input interval, a 2-vector specifying the interval of
%     approximation. The output satisfies size(d) = size(x). 
%
%     The optional input R is fed right into difference_stencil and has the same
%     meaning as in that function. 

global handles;
newton = handles.speclab.newton_polynomials;
fd = handles.finite_difference;

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
