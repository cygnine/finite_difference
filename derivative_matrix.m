function[mat] = derivative_matrix(x,k,varargin);
% [MAT] = DERIVATIVE_MATRIX(X,K,{R=0,PERIODIC=False});
%
%     Creates a sparse finite-difference matrix of order k on the 1D
%     mesh defined by the nodal locations x. The optional inputs R and PERIODIC
%     are the shift and periodicity flags, respectively, and serve the same
%     purpose as in difference_stencil, where they are explained. 

global handles;
newton = handles.speclab.NewtonPolynomials;
fd = handles.finite_difference;

% Create stencil
n = length(x);
stencil = fd.difference_stencil(n,k,varargin{:});

mat = spalloc(n,k+1,n*(k+1));

% This is the really slow, trivial way to do things
for q = 1:n
  y = zeros([n,1]);
  y(q) = 1;

  % Use stencil to compute interpolants
  dd = newton.divided_difference(x(stencil).',y(stencil.'));

  % Differentiate and evaluate the interpolants
  mat(:,q) = sparse(newton.newton_derivative_evaluate(x(stencil).',dd).');
end

