function Plot_Misc(varargin)

ax=get(gca);

fct=varargin{1};
switch upper(fct)
    case 'ID'
        if nargin<2, style='k:'; else style=varargin{2}; end
        if nargin<3, opt='single'; else opt=varargin{3}; end
        x(1)=max([ax.XLim(1) ax.YLim(1)]);
        x(2)=min([ax.XLim(2) ax.YLim(2)]);
        plot(x,x,style); bg;
        switch opt
            case 'single'
            case 'double', Per=2*pi;
                plot(x+Per,x    ,style); bg;
                plot(x    ,x+Per,style); bg;
                plot(x+Per,x+Per,style); bg;
            otherwise, error(['invalid opt: ' opt]);
        end
    case 'AXES'
        if nargin<2, style='k:'; else style=varargin{2}; end
        Plot_Misc('X-AXIS',style);
        Plot_Misc('Y-AXIS',style);
    case 'X-AXIS'
        if nargin<2, style='k:'; else style=varargin{2}; end
        plot([ax.XLim(1) ax.XLim(2)],[0 0],style); bg;
    case 'Y-AXIS'
        if nargin<2, style='k:'; else style=varargin{2}; end
        plot([0 0],[ax.YLim(1) ax.YLim(2)],style); bg;
    case {'V','VERTICAL'}
        x=varargin{2};
        if nargin<3, style='k:'; else style=varargin{3}; end
        for i=1:length(x)
            if nargin<4
                plot([x(i) x(i)],[ax.YLim(1) ax.YLim(2)],style);
            else
                plot([x(i) x(i)],[ax.YLim(1) ax.YLim(2)],style,'color',varargin{4});
            end
            bg;
        end
    case {'H','HORIZONTAL'}
        y=varargin{2};
        if nargin<3, style='k:'; else style=varargin{3}; end
        for i=1:length(y)
            plot([ax.XLim(1) ax.XLim(2)],[y(i) y(i)],style); bg;
        end
    case {'P','POLY'}
        p=varargin{2};
        if nargin<3, style='k:'; else style=varargin{3}; end
        if length(p)==2, x=[ax.XLim(1) ax.XLim(2)];
        else             x=linspace(ax.XLim(1),ax.XLim(2),100);
        end
        plot(x,polyval(p,x),style); bg;
    case {'ERR','ERROR'}
        x =varargin{2};  y=varargin{3};
        dx=varargin{4}; dy=varargin{5};
        if nargin<6, style='g-'; else style=varargin{6}; end
        if length(dx)>1, dx=median(dx); end
        if length(dy)>1, dy=median(dy); end
        plot([x-dx/2 x+dx/2],[y y],style);
        plot([x x],[y-dy/2 y+dy/2],style);
    otherwise
        warning(['invalid fct ' fct]);
end

end

function bg
set(gca,'children',circshift(get(gca,'children'),-1));
end

