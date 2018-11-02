function [Q, X, fcnt] = cumquad(funfcn,a,b,tol,trace,varargin)
%CUMQUAD   Numerically evaluate cumulative integral, adaptive Simpson quadrature (Adam Fudged).
%   [Q,X] = CUMQUAD(FUN,A,B) tries to approximate the cumulative integral of scalar-valued
%   function FUN from A to B to within an error of 1.e-6 using recursive
%   adaptive Simpson quadrature. FUN is a function handle. The function
%   Y=FUN(X) should accept a vector argument X and return a vector result
%   Y, the integrand evaluated at each element of X. The result is Q(X)
%   where Q(X) is the integral of FUN from a to X.
%
%   [Q,X] = CUMQUAD(FUN,A,B,TOL) uses an absolute error tolerance of TOL 
%   instead of the default, which is 1.e-6.  Larger values of TOL
%   result in fewer function evaluations and faster computation,
%   but less accurate results.  The QUAD function in MATLAB 5.3 used
%   a less reliable algorithm and a default tolerance of 1.e-3.
%
%   [Q,X] = CUMQUAD(FUN,A,B,TOL,TRACE) with non-zero TRACE shows the values
%   of [fcnt a b-a Q] during the recursion. Use [] as a placeholder to
%   obtain the default value of TOL.
%
%   [Q,X,FCNT] = CUMQUAD(...) returns the number of function evaluations.
%
%   Use array operators .*, ./ and .^ in the definition of FUN
%   so that it can be evaluated with a vector argument.
%
%   Notes:
%   Function CUMQUADL (TBA) may be more efficient with high accuracies and smooth
%   integrands.
%
%   Example:
%      [Q,X] = cumquad(@myfun,0,2);
%   where myfun.m is the M-file function:
%      %-------------------%
%      function y = myfun(x)
%      y = 1./(x.^3-2*x-5);
%      %-------------------%
%
%   or, use a parameter for the constant:
%      [Q,X] = cumquad(@(x)myfun2(x,5),0,2);
%   where myfun2 is the M-file function:
%      %----------------------%
%      function y = myfun2(x,c)
%      y = 1./(x.^3-2*x-c);
%      %----------------------%
%
%   Class support for inputs A, B, and the output of FUN: 
%      float: double, single
%
%   See also QUAD, QUADV, QUADL, DBLQUAD, TRIPLEQUAD, TRAPZ, FUNCTION_HANDLE.

%   Based on "adaptsim" by Walter Gander.  
%   Ref: W. Gander and W. Gautschi, "Adaptive Quadrature Revisited", 1998.
%   http://www.inf.ethz.ch/personal/gander
%   Copyright 1984-2004 The MathWorks, Inc. 
%   $Revision: 5.26.4.5 $  $Date: 2004/12/06 16:35:07 $
%   Fudged by Adam Wyatt to allow cumulative integrals at marginal extra
%   cost but without extra computation.

f = fcnchk(funfcn);
if nargin < 4 || isempty(tol), tol = 1.e-6; end;
if nargin < 5 || isempty(trace), trace = 0; end;

% Initialize with three unequal subintervals.
h = 0.13579*(b-a);
x = [a a+h a+2*h (a+b)/2 b-2*h b-h b];
y = f(x, varargin{:});
fcnt = 7;

% Fudge endpoints to avoid infinities.
if ~isfinite(y(1))
   y(1) = f(a+eps(superiorfloat(a,b))*(b-a),varargin{:});
   fcnt = fcnt+1;
end
if ~isfinite(y(7))
   y(7) = f(b-eps(superiorfloat(a,b))*(b-a),varargin{:});
   fcnt = fcnt+1;
end

% Call the recursive core integrator.
hmin = eps(b-a)/1024;
[Q1,X1,fcnt,warn(1)] = ...
   quadstep(f,x(1),x(3),y(1),y(2),y(3),tol,trace,fcnt,hmin,varargin{:});
[Q2,X2,fcnt,warn(2)] = ...
   quadstep(f,x(3),x(5),y(3),y(4),y(5),tol,trace,fcnt,hmin,varargin{:});
[Q3,X3,fcnt,warn(3)] = ...
   quadstep(f,x(5),x(7),y(5),y(6),y(7),tol,trace,fcnt,hmin,varargin{:});
Q = cumsum([Q1, Q2, Q3]);
X = [X1, X2, X3];
warn = max(warn);

switch warn
   case 1
      warning('MATLAB:quad:MinStepSize', ...
          'Minimum step size reached; singularity possible.')
   case 2
      warning('MATLAB:quad:MaxFcnCount', ...
          'Maximum function count exceeded; singularity likely.')
   case 3
      warning('MATLAB:quad:ImproperFcnValue', ...
          'Infinite or Not-a-Number function value encountered.')
   otherwise
      % No warning.
end

% ------------------------------------------------------------------------

function [Q,X,fcnt,warn] = quadstep (f,a,b,fa,fc,fb,tol,trace,fcnt,hmin,varargin)
%QUADSTEP  Recursive core routine for function QUAD.

maxfcnt = 10000;

% Evaluate integrand twice in interior of subinterval [a,b].
h = b - a;
c = (a + b)/2;
if abs(h) < hmin || c == a || c == b
   % Minimum step size reached; singularity possible.
   Q = h*fc;
   X = c;
   warn = 1;
   return
end
x = [(a + c)/2 (c + b)/2];
y = f(x, varargin{:}); 
fcnt = fcnt + 2;
if fcnt > maxfcnt
   % Maximum function count exceeded; singularity likely.
   Q = h*fc;
   X = c;
   warn = 2;
   return
end
fd = y(1);
fe = y(2);

% Three point Simpson's rule.
Q1 = (h/6)*(fa + 4*fc + fb);

% Five point double Simpson's rule.
%Q2 = (h/12)*(fa + 4*fd + 2*fc + 4*fe + fb);
Q =(h/12)*[(fa + 4*fd + fc), (fc + 4*fe + fb)];
Qav = sum(Q);

% One step of Romberg extrapolation.
Q2 = Qav + (Qav - Q1)/15;

X = [c, b];

if ~isfinite(Q)
   % Infinite or Not-a-Number function value encountered.
   warn = 3;
   return
end
if trace
   disp(sprintf('%8.0f %16.10f %18.8e %16.10f',fcnt,a,h,Q))
end

% Check accuracy of integral over this subinterval.
if abs(Q2 - Qav) <= tol
   warn = 0;
   return

% Subdivide into two subintervals.
else
   [Qac,Xac,fcnt,warnac] = quadstep(f,a,c,fa,fd,fc,tol,trace,fcnt,hmin,varargin{:});
   [Qcb,Xcb,fcnt,warncb] = quadstep(f,c,b,fc,fe,fb,tol,trace,fcnt,hmin,varargin{:});
   Q = [Qac, Qcb];
   X = [Xac, Xcb];
   warn = max(warnac,warncb);
end
