function[mat] = derivative_matrix_periodic(x,k,interval,varargin);
% [MAT] = DERIVATIVE_MATRIX_PERIODIC(X,K,INTERVAL,{R=0})
%
%     Creates a sparse finite-difference matrix of order k on the 1D mesh
%     defined by the nodal locations x with periodic continuation over the
%     interval specified by interval. The optional input R is the shift and
%     serves the same purpose as in difference_stencil, where it is explained. 

global handles;
newton = handles.speclab.NewtonPolynomials;
fd = handles.FiniteDifference;

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

% Allocation
mat = spalloc(n,n,n*(k+1));

% This is the really slow, trivial way to do things
for q = 1:n
  y = zeros([n,1]);
  y(q) = 1;

  % Use stencil to compute interpolants
  dd = newton.divided_difference(XInput.',y(stencil.'));

  % Differentiate and evaluate the interpolants
  mat(:,q) = sparse(newton.newton_derivative_evaluate(XInput.',dd).');
end
