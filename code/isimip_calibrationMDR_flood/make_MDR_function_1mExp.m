function mdr = make_MDR_function_1mExp(x,scale,shape)
% climada isimip create MDR function of x (hazard intensity) as a function
%   of scale and shape.
% MODULE:
%   isimip
% NAME:
%   make_MDR_function_1mExp
% PURPOSE:
%   create MDR function of x (hazard intensity) as a function
%   of scale and shape.
%   next call: isimip...
% CALLING SEQUENCE:
%   mdr_func=make_MDR_function_1mExp(scale,shape)
% EXAMPLE:
%   mdr_func=make_MDR_function_1mExp(0.5,0.5)
% INPUTS:
%   scale: scale parameter of function scale*(1-exp(-shape*x))
%   shape: scale parameter of function scale*(1-exp(-shape*x))
% OPTIONAL INPUT PARAMETERS:
%   none.
% OUTPUTS:
%   mdr: function of x which returns scale*(1-exp(-shape*x))
%
% MODIFICATION HISTORY:
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 20180911, initial
%-

% For reference, here the relevant example from https://ch.mathworks.com/help/matlab/matlab_prog/nested-functions.html
%make_MDR_function_1mExp returns a function of x that computes
%scale*(1-exp(-shape*x))
%   This function allows to obtain a function to evaluate
%   scale*(1-exp(-shape*x)) as a function of x only.
% function p = makeParabola(a,b,c)
% p = @parabola;
%    function y = parabola(x)
%    y = a*x.^2 + b*x + c;
%    end
% end
% The makeParabola function returns a handle to the parabola function that includes values for coefficients a, b, and c.
% At the command line, call the makeParabola function with coefficient values of 1.3, .2, and 30. Use the returned function handle p to evaluate the polynomial at a particular point:
% p = makeParabola(1.3,.2,30);
% X = 25;
% Y = p(X)
% mdr = @myfunc;
%     function y=myfunc(x)


mdr=scale*(1-exp(-shape*x));
%     end
end

