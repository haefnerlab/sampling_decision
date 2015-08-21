function state=Subplot(nn,i,ss, nyy)

% function Subplot(number of panels,
%                  index of current panel,
%                  desired ratio of vertical/horizontal panels)
% last parameter optional, =1 per default

persistent n;
persistent s;
persistent nx;
persistent ny;
persistent idx; % current panel

fig=gcf;

if nargin<1,     idx(fig.Number)=mod(idx(fig.Number),n(fig.Number))+1;
elseif nargin<2
    if isnumeric(nn), idx(fig.Number)=nn;
    else
        switch nn
            case 'h', idx(fig.Number)=mod(idx(fig.Number),n(fig.Number))+1;
            case 'v', idx(fig.Number)=mod(idx(fig.Number),n(fig.Number))+nx(fig.Number);
            case 'total', state=n(fig.Number);
            case 'current', state=idx(fig.Number);
            otherwise, error(nn);
        end
        switch nn
            case {'h','v'}
            otherwise, return;
        end
    end
else
    idx(fig.Number)=i;
    if length(s)<fig, s(fig.Number)=0; end
    if length(n)<fig, n(fig.Number)=0; end
    if nargin<3
        if isempty(s) || s(fig.Number)==0 || n(fig.Number)~=nn, ss=2/3;
        else           ss=s(fig.Number);
        end
        s(fig.Number)=ss;
    end
    if nargin<4
        if isempty(n) || isempty(n(fig.Number)) || n(fig.Number)~=nn || s(fig.Number)~=ss
            n(fig.Number)=nn;
            s(fig.Number)=ss;
            nx(fig.Number)=round(sqrt(n(fig.Number)/s(fig.Number))); if nx(fig.Number)<1, nx(fig.Number)=1; end
            ny(fig.Number)=ceil(n(fig.Number)/nx(fig.Number));
            idx(fig.Number)=1;
            changed=1;
            while changed, changed=0;
                if nx(fig.Number)>ny(fig.Number) && (nx(fig.Number)-1)*ny(fig.Number)>=nn, nx(fig.Number)=nx(fig.Number)-1; changed=1; end
                if ny(fig.Number)>nx(fig.Number) && (ny(fig.Number)-1)*nx(fig.Number)>=nn, ny(fig.Number)=ny(fig.Number)-1; changed=1; end
            end
        end
    else
        nx(fig.Number)=ss;
        ny(fig.Number)=nyy;
        n(fig.Number)=nn;
        if nx(fig.Number)*ny(fig.Number)<nn
            warning('nx(fig.Number)*ny(fig.Number)<n(fig.Number)! continuing for the time being' );
        end
    end
end

if ~isempty(idx(fig.Number)) & idx(fig.Number)>0 & idx(fig.Number)<=n(fig.Number), subplot(ny(fig.Number),nx(fig.Number),idx(fig.Number));
elseif idx(fig.Number)~=0
    warning('Subplot: can''t do what you want!');
end

hold on;

state.n=n(fig.Number);
%state.s=s(fig.Number);
state.nx=nx(fig.Number);
state.ny=ny(fig.Number);
state.idx=idx(fig.Number);
