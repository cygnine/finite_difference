% Example script for the finite_difference module

% This module creates finite-difference derivative matrices of any order on
% unstructured meshes. Effectively, this is the same as using
% poly_differentiation in the piecewise_polynomial module.

clear; close all;
global handles;
fd = handles.finite_difference;

% Function to differentiate
f = @(x) exp(sin(x));
df = @(x) cos(x).*exp(sin(x));

% Data points
x = linspace(-3,3,200).'; 
y = f(x);  % note that y must be a column matrix if we want to do matvecs

% Form the finite-difference matrix with quadratic approximation (same as
% centralDmat).
dmat2 = fd.derivative_matrix(x,2);

dmat7 = fd.derivative_matrix(x,7);

figure; 
subplot(1,2,1); spy(dmat2);
title('Sparsity pattern for dmat2')
subplot(1,2,2); spy(dmat7);
title('Sparsity pattern for dmat7')

figure;
plot(x,dmat2*y, 'k-.', x, dmat7*y, 'r-.', x, df(x), 'b.');
xlabel('x');
ylabel('Derivative approximations');
legend('Quadratic', '7th order', 'Actual derivative');

% You can check that dmat2 and centralDmat are the same matrix. Also, you can
% check visually with figure 2 that dmat7 is much more accurate, but this is not
% necessarily the case if the data is unstructured (non-equispaced) or the
% function does not look smooth.

% Note that this package also has a function difference_derivative that
% effectively computes dmat*y, i.e. it does not return the matrix, just the
% action of it. This is *much* faster than forming the matrix. In addition,
% although this is merely a special case of
% piecewise_polynomial.poly_differentiation (calling z the same as x),
% difference_derivative is also much faster than this function.
