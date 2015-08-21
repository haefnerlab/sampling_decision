function ProgressReport(time,index,loopsize)

% function ProgressReport(time,index,loopsize)
%
% remains silent if matlabpool size >0

global PARALLEL;

persistent lasttime;

%if matlabpool('size')==0 % useless since clients report 0, too!
if isempty(PARALLEL), PARALLEL=0; end

if PARALLEL==0
    
    if isempty(lasttime)
        try toc, catch err, tic; end % call tic if it hasn't been called before
        lasttime=toc;
    end
    
    if toc-lasttime>time
        p=index(end);
        for i=1:length(index)-1
            p=p+(index(i)-1)*prod(loopsize(i+1:end));
        end
        p=p/prod(loopsize);
        fprintf([num2str(round(p*100)) '%% ']);
        lasttime=toc;
    elseif toc-lasttime<0 % a tic was run in the meantime somewhere
        lasttime=toc;
    end
end
