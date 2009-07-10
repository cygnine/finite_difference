% Debugging file for difference_derivative

clear
cd ..

% For making more random value to start
RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));

%%%%%%%%%%%%%% Try derivative evaluation of known function (nonperiodic)
if true
  f = @(x) sin(x);
  df = @(x) cos(x);

  N = 100;
  xmin = 0;
  xmax = 2*pi;
  x = sort((xmax-xmin)*rand([N,1]) + xmin);
  fx = f(x);

  k = 4;

  dfx = difference_derivative(x,fx,k);

  fprintf('Derivative error is %3.3e\n', norm(dfx-df(x)));

  % While we're here, let's test making the sparse differentiation matrix, too
  spmat = derivative_matrix(x,k);

  dfx2 = spmat*fx;

  fprintf('Derivative error is %3.3e\n', norm(dfx2-df(x)));
end

%%%%%%%%%%%%%% Try derivative evaluation of known function (periodic)
if true
  f = @(x) sin(x);
  df = @(x) cos(x);

  N = 100;
  xmin = 0;
  xmax = 2*pi;

  x = sort((xmax-xmin)*rand([N,1]) + xmin);
  fx = f(x);

  k = 4;

  dfx_periodic = difference_derivative_periodic(x,fx,k,[xmin,xmax]);

  fprintf('Derivative error is %3.3e\n', norm(dfx_periodic-df(x)));

  % While we're here, let's test making the sparse differentiation matrix, too
  spmat_periodic = derivative_matrix_periodic(x,k,[xmin,xmax]);

  dfx2_periodic = spmat_periodic*fx;

  fprintf('Derivative error is %3.3e\n', norm(dfx2_periodic-df(x)));
end

cd debug
