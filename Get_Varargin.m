function varargout = Get_Varargin(varargin)

% function varargout = Get_Varargin(varargin)
% usage:
% first call: [a b c etc]=Get_Varargin(varargin);
% subsequent calls: [d e f etc]=Get_Varargin;
%
% by Ralf M. Haefner, 22/5/2012

persistent next inputs;

if isempty(next), next=1; end

if nargin>0
  if iscell(varargin{1})
    inputs=varargin{1};
    if nargin<2, next=1; else next=varargin{2}; end
  else
    switch varargin{1}
      case 'GVNext', varargout{1}=next; return;
      otherwise
        error('invalid input');
    end
  end
end 

varargout=cell(1,nargout);
nin =length(inputs);
if next+nargout-1<=nin
  for i=1:nargout
    varargout{i}=inputs{i-1+next};
  end
else
  %noutin=min([nout nargout-next+1]);
  noutin=nin-next+1;
  for i=1:noutin
    varargout{i}=inputs{i-1+next};
  end
  for i=max([1 noutin+1]):nargout
    varargout{i}=[];
  end
end
next=next+nargout;
  
end
