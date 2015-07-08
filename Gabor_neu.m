function y = Gabor_neu(para,x,fct,ana)

% same as Gabor_inv but using sigma and freq, not omega and sigma^2
% function y = Gabor_inv(para,x)
% fct is optional, default value='rect'
% fct='rect': return rectified value
% fct='sqrt': return rectified square-root
% fct='orig': return unrectified Gabor
% fct='orig-sqrt': return signed sqrt of Gabor

persistent sigma;

if nargin<3, fct='rect';  end
if nargin<4, ana='plain'; end

if nargin==1 % set sigma's
  sigma=para;
  y=-1;
  %disp(['Gabor_2_AC: sigma set to: ' num2str(sigma)]);
else

  if isempty(sigma)
    sigma=1;
  end

  switch ana
    case 'plain'
      y=para(1)+para(2)*cos(2*pi*para(3)*(x-para(5))-para(4)).*exp(-(x-para(5)).^2/2/para(6).^2);
      y=y./sigma;
      switch lower(fct)
        case 'rect'
          y(y<0)=0; % simple rectification
        case 'sqrt'
          y(y<0)=0;
          y=sqrt(y);
        case 'orig'
        case 'orig-sqrt'
          y=Pot_Signed(y,0.5);
        otherwise
          warning(['returning original Gabor instead of ' fct]);
      end
      
    case 'derivative'
      if isempty(sigma) || sigma~=1
        sigma=1;
        disp('Gabor_neu: sigma set to 1 for derivative!');
      end
      %y=-para(2)*sin(2*pi*para(3)*(x-para(5))-para(4))*2*pi.*exp(-(x-para(5)).^2/2/para(6).^2)+...
      %   para(2)*cos(2*pi*para(3)*(x-para(5))-para(4)).*exp(-(x-para(5)).^2/2/para(6).^2).*(x-para(5))/para(6);
      % above seems wrong! (independent of spatial frequency)
      y=-para(2)*exp(-(x-para(5)).^2/2/para(6)^2)/para(6)^2.*...
        (2*pi*para(3)*para(6)^2*sin(2*pi*para(3)*(x-para(5))-para(4))+(x-para(5)).*cos(2*pi*para(3)*(x-para(5))-para(4)));
    otherwise
      error(ana);
  end
end
