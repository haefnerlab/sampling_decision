function y = Cut_Gaussian(fct,mu,sigma,x)

% function y = Cut_Gaussian(fct,mu,sigma,x)
% for fct='random', x=size(samples)
% is x is omitted, it is assumed to be 1

if nargin<4, x=1; end

if mu/sigma<-35
  %disp(['mu/sigma=' num2str(mu/sigma) '<-30: approx by exponential']);
  y=exprnd(-sigma^2/mu,x);

else
  
  compC=normcdf(0,-mu,sigma); % =1-comp
  
  switch fct
    case 'pdf'
      y=normpdf(x,mu,sigma)/compC;
      y(x<0)=0;
    case 'cdf'
      comp =normcdf(0, mu,sigma);
      y=(normcdf(x,mu,sigma)-comp)/compC;
      y(x<0)=0;
    case 'inv'
      %y=norminv(compC*x+comp,mu,sigma);
      y=-norminv(compC*x,mu,sigma)+2*mu;
      y(x<0 | x>1)=NaN;
    case 'random'
      y=Cut_Gaussian('inv',mu,sigma,rand(x));
    otherwise
      error(fct);
  end
end