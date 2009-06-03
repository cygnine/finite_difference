% MATLAB File : DividedDifferenceDebug.m
%
% * Creation Date : 2009-06-03
%
% * Last Modified : Wed 03 Jun 2009 02:16:38 PM EDT
%
% * Created By : Akil Narayan
%
% * Purpose : Debugging DivdedDifference
%

clear
cd ..

dbstop if error

RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));

%%%%%%%%%%%%%%%% ORDERED INTERPOLATION
if false
  % Test set:
  x = [1,2,3,4,5];
  y = [1,-2,4,-1,0];

  % Calculate modal coefficients
  dd = DividedDifference(x,y);

  % Plot data points, interpolant
  z = linspace(min(x),max(x),100);
  zy = NewtonHornerEval(z,dd,x);

  plot(z,zy,'b',x,y,'r.');
end

%%%%%%%%%%%%%%%% UNORDERED INTERPOLATION
if false
  % Test set:
  x = [1,3,5,2,4];
  y = [1,4,0,-2,-1];

  x = rand([5,1]);
  y = rand([5,1]);

  % Calculate modal coefficients
  dd = DividedDifference(x,y);

  % Plot data points, interpolant
  z = linspace(min(x),max(x),100);
  zy = NewtonHornerEval(z,dd,x);

  plot(z,zy,'b',x,y,'r.');
end

%%%%%%%%%%%%%%%% UNORDERED INTERPOLATION, VECTORIZED
if true
  % Test set:
  x = rand([5,5]);
  y = rand([5,5]);

  % Calculate modal coefficients
  dd = DividedDifference(x,y);

  % Plot data points, interpolant
  
  figure; hold on;
  for q = 1:size(x,2)
    z = linspace(min(x(:,q)),max(x(:,q)),100);
    zy = NewtonHornerEval(z,dd(:,q),x(:,q));
    plot(z,zy,'b');
  end

  plot(x,y,'r.');
end

cd debug
