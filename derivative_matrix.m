function[mat] = derivative_matrix(x,k,varargin);
% derivative_matrix -- finite-difference derivative matrix
%
% [mat] = derivative_matrix(x,k,{r=0,periodic=false});
%
%     Creates a sparse finite-difference matrix of order k on the 1D
%     mesh defined by the nodal locations x. The optional inputs r and periodic
%     are the shift and periodicity flags, respectively, and serve the same
%     purpose as in difference_stencil, where they are explained. 

persistent difference_stencil divided_difference newton_derivative_evaluate
if isempty(difference_stencil)
  from finite_difference import difference_stencil
  from speclab.newton_polynomials import divided_difference newton_derivative_evaluate
end

% Create stencil
n = length(x);
stencil = difference_stencil(n,k,varargin{:});

mat = spalloc(n,k+1,n*(k+1));

% This is the really slow, trivial way to do things
for q = 1:n
  y = zeros([n,1]);
  y(q) = 1;

  % Use stencil to compute interpolants
  dd = divided_difference(x(stencil).',y(stencil.'));

  % Differentiate and evaluate the interpolants
  mat(:,q) = sparse(newton_derivative_evaluate(x(stencil).',dd).');
end
