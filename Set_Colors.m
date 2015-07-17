function col = Set_Colors(n,repeat)

% function col = Set_Colors(n)

if nargin<1, n=7; end
if nargin<2, repeat=1; end

i=1;
if n<=7
  col={'b','r','g','c','m','y','k'};
elseif n<=25
  for r=1:-0.5:0
    for g=1:-0.5:0
      for b=1:-0.5:0
        col{i}=[r g b];
        i=i+1;
      end
    end
  end
  col=col(3:end); % first two are white and almost white
elseif n<=57
  for r=1:-1/3:0
    for g=1:-1/3:0
      for b=1:-1/3:0
        col{i}=[r g b];
        i=i+1;
      end
    end
  end
  col=col([5:16 20:end]); % first five are white and almost white
else
  step=1/(n+10)^(1/3);
  for r=1:-step:0
    for g=1:-step:0
      for b=1:-step:0
        col{i}=[r g b];
        i=i+1;
      end
    end
  end
end

for i=2:repeat
  col={col{:} col{:}};
end