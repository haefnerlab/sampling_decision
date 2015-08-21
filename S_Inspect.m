function out = S_Inspect(fct,varargin)

% function out = S_Inspect(fct,varargin)
% 'Exp-X-Simple': input E.X, Histograms for each Xs across repetitions
% 'Exp-LOGS' : detailed inspection for LOGS experiments, text output
%     panel 1:
% 'Exp-X'    :



col='brgcmyk'; col=[col col col]; col=[col col col]; col=[col col col];

switch fct
    
    case 'Exp-All' % everything, e.g. for debugging
        [E thresh]=Get_Varargin(varargin);
        if isempty(thresh), thresh=0.9; end
        if thresh>1, warning('threshhold must be <=1!'); end
        
        [nrep nL nG ntime]=size(E.G);
        nX=size(E.X,2);
        iL=1;
        start=5;
        n0S=E.InputImage.n_zero_signal; % HACK! After 2013/7/11 remove -1!
        Get_Figure('SI:Exp-1');
        Subplot(24,0,6,4);
        col=Set_Colors(4);    style={'-','-','--','--'};
        Subplot; title('X'); xlabel('time');
        imagesc(squeeze(mean(E.X(:,:,start:end),1))); colorbar;
        Subplot; title('G'); xlabel('time');
        imagesc(squeeze(mean(E.G(:,iL,:,start:end),1))); colorbar;
        if length(size(E.Signal))>2
            Subplot; title('Signal'); xlabel('\phi');
            plot(squeeze(mean(E.Signal(:,:,E.Sampling.access(start:10)),1)),'b-');
            plot(squeeze(mean(E.Signal(:,:,E.Sampling.access(11:ntime)),1)),'r-');
            Subplot; title({'mean signal (conditioned) before signal onset','blue=Task1, blue=Task2'});
            s=permute(E.Signal(:,:,1:n0S),[2 1 3]);s=s(:,:); % combine last 2 dimensions
            sm=mean(s,2); plot(sm,'g-');
            i=1; % only vertical vs horizontal stimuli possible at the moment
            idx=find(s((i-1)*nX/4+1,:)-s((i+1)*nX/4+1,:)>0);
            plot(mean(s(:,idx),2),[col{i} '-']);
            plot(mean(s(:,idx),2)-sm,[col{i} '--']);
            Subplot; title('mean signal (conditioned) after signal onset');
            s=permute(E.Signal(:,:,n0S+1:end),[2 1 3]);s=s(:,:); % combine last 2 dimensions
            sm=mean(s,2); plot(sm,'g-');
            i=1; % only vertical vs horizontal stimuli possible at the moment
            idx=find(s((i-1)*nX/4+1,:)-s((i+1)*nX/4+1,:)>0);
            plot(mean(s(:,idx),2),[col{i} '-']);
            plot(mean(s(:,idx),2)-sm,[col{i} '--']);
            Subplot; title('signals before signal onset');
            s=permute(E.Signal(:,:,1:n0S),[2 1 3]);s=s(:,:); % combine last 2 dimensions
            for i=1:4 % 4 possible signals
                %Subplot; title('signals');
                if numel(unique(s((i-1)*nX/4+1,:)))>1
                    Histogram(s((i-1)*nX/4+1,1:n0S),[],[],[col{i} style{i}]);
                end
            end
            Subplot; title('signals after signal onset');
            s=permute(E.Signal(:,:,n0S+1:end),[2 1 3]);s=s(:,:); % combine last 2 dimensions
            for i=1:4 % 4 possible signals
                %Subplot; title('signals');
                if numel(unique(s((i-1)*nX/4+1,:)))>1
                    Histogram(s((i-1)*nX/4+1,:),[],[],[col{i} style{i}]);
                end
            end
        end
        Subplot; title('mean X'); xlabel('\phi');
        plot(squeeze(mean(mean(E.X(1,:,start:n0S),1),3)),'-');
        plot(squeeze(mean(mean(E.X(1,:,n0S+1:ntime),1),3)),'r-');
        Subplot; title('mean G'); xlabel('\phi');
        plot(squeeze(mean(mean(E.G(:,iL,:,start:n0S  ),1),4)),'-');
        plot(squeeze(mean(mean(E.G(:,iL,:,n0S+1:ntime),1),4)),'r-');
        Subplot; title('Task'); xlabel('time');
        plot(squeeze(mean(E.T(:,1,:))),'r*');
        plot(squeeze(mean(E.T(:,2:end,:))),'-');
        Plot_Misc('v',n0S+0.5);
        Subplot; title('O'); xlabel('time');
        plot(squeeze(mean(E.O(:,1,:))),'r*');
        plot(squeeze(mean(E.O(:,2:end,:))),'-');
        Plot_Misc('v',n0S+0.5);
        Subplot; title('S'); xlabel('time');
        plot(mean(E.S,1),'.');
        Plot_Misc('v',n0S+0.5);
        Subplot; title('mean(X)'); xlabel('time');
        plot(squeeze(mean(mean(E.X,1))),'b-');
        plot(squeeze(mean(mean(E.X(:,1:nX/2,:),1))),'r-');
        plot(squeeze(mean(mean(E.X(:,nX/2+1:nX,:),1))),'r--');
        plot(squeeze(mean(mean(E.X(:,[1:nX/4 3/4*nX+1:nX],:),1))),'g-');
        plot(squeeze(mean(mean(E.X(:,nX/4+1:3/4*nX,:),1))),'g--');
        Plot_Misc('v',n0S+0.5);
        Subplot; title('posterior(O)'); xlabel('time');
        plot(squeeze(E.O(1:min([100 size(E.O,1)]),2,:))');
        Subplot; title('logodds(O)'); xlabel('time');
        plot(squeeze(diff(log(E.O(1:min([100 size(E.O,1)]),2:3,:)),[],2))');
        Subplot; title('RT'); xlabel('time from signal onset');
        rt=zeros(1,nrep);
        for i=1:nrep
            pmax=max(E.O(i,2:3,n0S+1:ntime),[],2);
            aux=1+find(pmax<thresh,1,'last');
            if isempty(aux), rt(i)=n0S;
            else             rt(i)=aux+n0S;
            end
        end
        %bar(1:(ntime+1),hist(rt,0.5:ntime+0.5),'barwidth',1);
        plot(1:ntime+1,cumsum(histc(rt,1:ntime+1))/nrep,'b-'); % data
        gauss=abs(cumsum(randn(nrep,ntime-n0S),2));
        gauss_thresh=prctile(gauss(:,end),100*sum(rt==ntime+1)/nrep);
        gauss_rt=zeros(1,nrep);
        for i=1:nrep
            aux=1+find(gauss(i,:)<gauss_thresh,1,'last');
            if isempty(aux), gauss_rt(i)=n0S; else gauss_rt(i)=aux+n0S; end
        end
        plot(1:ntime+1,cumsum(histc(gauss_rt,1:ntime+1))/nrep,'r-'); % normal random walk
        legend('data','Gauss');
        Subplot; title('p(O==1)(t)');
        p_edges=linspace(0,1.0001,20);
        post=histc(squeeze(E.O(:,2,:)),p_edges);
        for i=1:ntime, post(:,i)=post(:,i)/max(post(:,i)); end
        imagesc(1:ntime,p_edges(1:end-1),post(1:end-1,:)); colorbar;
        xlabel('time'); ylabel('rel frequency per time');
        Subplot; title('logodds(O)(t)');
        p_edges=linspace(-3,3,21); % #odd to include zero in bin
        post=histc(squeeze(diff(log(E.O(:,2:3,:)),[],2)),p_edges);
        for i=1:ntime, post(:,i)=post(:,i)/max(post(:,i)); end
        imagesc(1:ntime,p_edges(1:end-1),post(1:end-1,:)); colorbar;
        xlabel('time'); ylabel('rel frequency per time');
        Subplot; title('var of logodds');
        plot(1:ntime,var(squeeze(diff(log(E.O(:,2:3,:)),[],2))),'b-');
        plot(n0S+1:ntime,var(E.Projection.odds_inc*cumsum(randn(nrep,ntime-n0S),2)),'r-');
        plot(n0S+1:ntime,E.Projection.odds_inc*(1:ntime-n0S),'g-');
        legend('data','Gauss');
        Subplot; title('max(p(O))(t)');
        p_edges=linspace(0.5,1.0001,20);
        post=histc(squeeze(max(E.O(:,2:3,:),[],2)),p_edges);
        for i=1:ntime, post(:,i)=post(:,i)/max(post(:,i)); end
        imagesc(1:ntime,p_edges(1:end-1),post(1:end-1,:)); colorbar;
        xlabel('time'); ylabel('rel frequency per time');
        Subplot; title('cumsum(max(p(O)))(t)');
        post=histc(squeeze(max(E.O(:,2:3,:),[],2)),p_edges);
        imagesc(1:ntime,p_edges(end-1:-1:1),cumsum(post(end-1:-1:1,:),1)); colorbar;
        xlabel('time'); ylabel('rel frequency per time');
        Subplot; title('decision');
        bar([sum(E.O(:,2,ntime)>0.5) ...
            sum(E.O(:,3,ntime)>0.5)]);
        Subplot; title('PK'); xlabel('time');
        O_pref=1; ixp=1; ixa=1+nX/2;
        idx_pref=(E.O(:,2,end)>0.5);
        idx_anti=(E.O(:,3,end)>0.5);
        prefpref=mean(E.Signal(idx_pref,ixp,:));
        prefanti=mean(E.Signal(idx_pref,ixa,:));
        antipref=mean(E.Signal(idx_anti,ixp,:));
        antianti=mean(E.Signal(idx_anti,ixa,:));
        kernel=prefpref-prefanti-antipref+antianti;
        plot(squeeze(kernel),'b-');
        Plot_Misc('v',n0S+0.5);
        Plot_Misc('v',n0S+5+0.5);
        Subplot; title('PK(RT>end-5)'); xlabel('time');
        idx_pref=(E.O(:,2,end)>0.5) & (rt'>length(pmax)-5);
        idx_anti=(E.O(:,3,end)>0.5) & (rt'>length(pmax)-5);
        prefpref=mean(E.Signal(idx_pref,ixp,:));
        prefanti=mean(E.Signal(idx_pref,ixa,:));
        antipref=mean(E.Signal(idx_anti,ixp,:));
        antianti=mean(E.Signal(idx_anti,ixa,:));
        kernel=prefpref-prefanti-antipref+antianti;
        plot(squeeze(kernel),'b-');
        Plot_Misc('v',n0S+0.5);
        Plot_Misc('v',n0S+5+0.5);
        
        for i=1:24, Subplot(i); axis tight; end
        
    case 'corr-temp' % everything, e.g. for debugging
        E=Get_Varargin(varargin);
        
        
    case 'Exp-X-Simple' % Histograms for each Xs across repetitions
        X=varargin{1};
        [nrep nX nsamples]=size(X);
        Get_Figure(['SI:' fct]);
        X=reshape(shiftdim(X,1),[nX,nsamples*nrep]);
        col=Set_Colors(nX);
        for i=1:nX
            dx=5*diff(prctile(X(i,:),[16 84]))/sqrt(size(X,2));
            Histogram(X(i,:),dx/5,[dx 0 max(X(:))],col{i});
            Plot_Misc('v',prctile(X(i,:),50),'--',col{i});
            Plot_Misc('v',mean(X(i,:))      ,'-' ,col{i});
        end
        %xlim([0 max(prctile(X',95))]);
        
    case 'Exp-LOGS'
        Exp=varargin{1};
        disp('X:');
        disp(mean(Exp.X,3));
        disp('G:');
        disp(squeeze(mean(mean(Exp.G,4))));
        disp('O:');
        disp(mean(Exp.O,2));
        for i=1:size(Exp.O,1), disp([sum(Exp.O(i,:)==1) sum(Exp.O(i,:)==2)]); end
        disp('L:');
        disp(mean(Exp.L,2));
        for i=1:size(Exp.L,1),
            disp([sum(Exp.L(i,:)==1) sum(Exp.L(i,:)==2) sum(Exp.L(i,:)==3)]);
        end
        disp('S:');
        disp(mean(Exp.S,2));
        for i=1:size(Exp.S,1), disp(prctile(Exp.S(i,:),[5 16 50 84 95])); end
        
    case 'Exp-OG'
        Exp=varargin{1};
        if isfield(Exp,'end')
            X=Exp.X(1:Exp.end,:,:); G=Exp.G(1:Exp.end,:,:); O=Exp.O(1:Exp.end,:);
        else
            X=Exp.X; G=Exp.G; O=Exp.O;
        end
        [nr nX ns]=size(X);
        nG=size(G,2);
        Get_Figure(['SI:' fct]);
        Subplot(4,1,2,2); title('corr(X)');
        cX=zeros(nr,nX,nX);
        for i=1:nr
            cX(i,:,:)=corr(squeeze(X(i,:,:))',squeeze(X(i,:,:))');
        end
        cX=squeeze(mean(cX)); for j=1:nX, cX(j,j)=0; end
        phi_X=Exp.P.phi_x;
        imagesc(phi_X,phi_X,cX); axis image; colorbar;
        xlabel('\Delta\phi'); ylabel('\Delta\phi');
        Subplot; title('corr(G)');
        cG=zeros(nr,nG,nG);
        for i=1:nr
            cG(i,:,:)=corr(squeeze(G(i,:,:))',squeeze(G(i,:,:))');
        end
        cG=squeeze(mean(cG)); for j=1:nG, cG(j,j)=0; end
        phi_G=[0:nG-1]/nG*pi;
        imagesc(phi_G,phi_G,cG); axis image; colorbar;
        xlabel('\phi'); ylabel('\phi');
        Subplot; title('Toeplitz profiles');
        csame=(cX(1:nX/2,1:nX/2)+cX(nX/2+1:nX,nX/2+1:nX))/2;
        phi=Exp.P.phi_x(1:nX/2);
        [profile para]=Fit_Toeplitz(phi,csame,'poly',1);
        plot(phi,profile,'b-');
        plot(phi,polyval(para,phi),'b--');
        cdiff=(cX(1:nX/2,1+nX/2:nX)+cX(nX/2+1:nX,1:nX/2)')/2;
        phi=fliplr(Exp.P.phi_x(1:nX/2)); %-mod(:,pi)
        [profile para]=Fit_Toeplitz(phi,cdiff,'poly',1);
        plot(phi,profile,'r-');
        plot(phi,polyval(para,phi),'r--');
        xlabel('\Delta\phi')
        
        
    case 'Exp-X'
        X=varargin{1};
        [nr nX ns]=size(X);
        Get_Figure(['SI:' fct]);
        Subplot(4,1,2,2); title('combined');
        Xt=reshape(shiftdim(X,1),[nX,ns*nr]);
        for i=1:4
            dx=10*diff(prctile(Xt(i,:),[16 84]))/sqrt(ns*nr);
            Histogram(Xt(i,:),dx/5,[dx 0 max(Xt(i,:))],col(i));
            Plot_Misc('v',prctile(Xt(i,:),50),[col(i) '--']);
            Plot_Misc('v',mean(Xt(i,:))      ,[col(i)  '-']);
        end
        xlim([0 max(prctile(Xt',95))]);
        Subplot(2); title('mean');
        Xt=squeeze(mean(X,3));
        for i=1:4
            dx=10*diff(prctile(Xt(:,i),[16 84]))/(ns*nr)^0.4;
            Histogram(Xt(:,i),dx/5,[dx 0 max(Xt(:,i))],col(i));
            Plot_Misc('v',prctile(Xt(:,i),50),[col(i) '--']);
            Plot_Misc('v',mean(Xt(:,i))      ,[col(i)  '-']);
        end
        xlim([0 max(prctile(Xt',95))]);
        Subplot(8,5,2,4);
        pair=[[1 3];[2 4];[1 2];[3 4]]; sp=[5 5 7 7 5 7]; co=[1 2 1 2 3 3];
        for i=1:4, Subplot(sp(i));
            aux=squeeze(mean(X(:,pair(i,1),:)-X(:,pair(i,2),:),3));
            dx=5/sqrt(nr);
            Histogram(aux,dx/5,dx,col(co(i)));
            disp(['aux: ' num2str([pair(i,:) mean(aux>0) std(aux>0) prctile(aux,[0 16 50 84 100])])]);
        end
        if nX>4
            aux=zeros(1,nr);
            for i=5:nX, Subplot(sp(i));
                for j=1:nr, aux(j)=sum(X(j,i,:)==1)/ns; end
                x=5/sqrt(nr);
                Histogram(aux-0.5,dx/5,dx,col(co(i)));
                disp(['aux: ' num2str([i mean(aux>0.5) std(aux>0.5) prctile(aux,[0 16 50 84 100])])]);
            end
        end
        
        Subplot(16,11,4,4);
        pair=[[1 2];[1 3];[1 4]; [3 4]]; sp=[11 12 15 16];
        for i=1:4, Subplot(sp(i));
            x=mean(X(:,pair(i,1),:),3)-mean(mean(X(:,pair(i,1),:),1),3);
            y=mean(X(:,pair(i,2),:),3)-mean(mean(X(:,pair(i,2),:),1),3);
            scatter(x,y,4,'b','.'); axis tight;
            [c p]=corrcoef(x,y);
            title([Num2StrV(pair(i,:),'f0') '|' Num2StrV(100*[c(1,2) p(1,2)],'f0') '%']);
        end
        
    case 'Attention-X'
        olist=varargin{1};
        tau=S_Retrieve('tau',olist); ntau=length(tau);
        pL=S_Retrieve('pL',olist);   npL=length(pL);
        c1=S_Retrieve('c',olist,1);  nc1=length(c1);
        c2=S_Retrieve('c',olist,4);  nc2=length(c2);
        fig=Get_Figure(['SI' fct]);
        xp=zeros(1000,ntau,nc1,npL);
        xa=zeros(1000,ntau,nc1,npL);
        for it=1:ntau
            for ic=1:nc1
                for ip=1:npL
                    aux=S_Retrieve('1',olist,tau(it),c1(ic),0,pL(ip));
                    xp(:,it,ic,ip)=mean(aux.X(:,1,:),3);
                    xa(:,it,ic,ip)=mean(aux.X(:,2,:),3);
                    vp(:,it,ic,ip)=var(aux.X(:,1,:),[],3)./xp(:,it,ic,ip);
                    va(:,it,ic,ip)=var(aux.X(:,2,:),[],3)./xa(:,it,ic,ip);
                end
            end
        end
        Subplot(8*ntau,0,ntau,8);
        r=zeros(size(xp));
        % pref
        for it=1:ntau, Subplot; axis([0 1 0 2.3]); title(['\tau=' num2str(tau(it))]);
            for ic=1:nc1
                for ip=1:npL, r(:,it,ic,ip)=xp(:,it,ic,ip)./xp(:,it,ic,npL); end
                errorbar(pL,mean(r(:,it,ic,:)),std(r(:,it,ic,:))/sqrt(1000),[col(ic) '.-']);
            end
        end
        for it=1:ntau, Subplot;
            imagesc(squeeze(mean(r(:,it,:,:))),[0 2.2]); axis tight;
        end
        % anti-pref
        for it=1:ntau, Subplot; axis([0 1 0.6 1.7]);
            for ic=1:nc1
                for ip=1:npL, r(:,it,ic,ip)=xa(:,it,ic,ip)./xa(:,it,ic,npL); end
                errorbar(pL,mean(r(:,it,ic,:)),std(r(:,it,ic,:))/sqrt(1000),[col(ic) '.-']);
            end
        end
        for it=1:ntau, Subplot;
            imagesc(squeeze(mean(r(:,it,:,:))),[0 2.2]); axis tight;
        end
        % pref
        for it=1:ntau, Subplot; %axis([0 1 0 2.3]);
            for ic=1:nc1
                for ip=1:npL, r(:,it,ic,ip)=vp(:,it,ic,ip)./vp(:,it,ic,npL); end
                errorbar(pL,mean(r(:,it,ic,:)),std(r(:,it,ic,:))/sqrt(1000),[col(ic) '.-']);
            end
        end
        for it=1:ntau, Subplot;
            imagesc(squeeze(mean(r(:,it,:,:))),[0 2.2]); axis tight;
        end
        % anti-pref
        for it=1:ntau, Subplot; %axis([0 1 0.6 1.7]);
            for ic=1:nc1
                for ip=1:npL, r(:,it,ic,ip)=va(:,it,ic,ip)./va(:,it,ic,npL); end
                errorbar(pL,mean(r(:,it,ic,:)),std(r(:,it,ic,:))/sqrt(1000),[col(ic) '.-']);
            end
        end
        for it=1:ntau, Subplot;
            imagesc(squeeze(mean(r(:,it,:,:)))); axis tight;
        end
        
    case 'Psychometry'
        olist=varargin{1};
        tau=S_Retrieve('tau',olist); ntau=length(tau);
        pL=S_Retrieve('pL',olist);
        pL=pL(pL>=0.5); npL=length(pL);
        c1=S_Retrieve('c',olist,1);  nc1=length(c1);
        c2=S_Retrieve('c',olist,4);  nc2=length(c2);
        fig=Get_Figure(['PNIPS11' fct]);
        %Subplot(3,1,3,1);
        pair=[[1 2];[3 4];[1 3];[2 4];[1 4];[2 3]];
        cor=zeros(ntau,nc1,1,npL,size(pair,1)); sig=zeros(size(cor));
        col='b..rgrrr';
        Subplot(4,0,2,2);
        for it=[1 4] %1:length(tau),
            Subplot; title(tau(it)); axis([0 5 0.5 1]);
            for ip=[1 npL]% npL] %1:length(pL)
                %Subplot; title(num2str([tau(it) pL(ip)]));
                if 1 % Attended condition
                    for ic1=1:length(c1), ic2=1;
                        o=S_Retrieve('get-one',olist,tau(it),c1(ic1),c2(ic2),pL(ip));
                        X=o.X;
                        yx(ic1)=sum(sign(squeeze(mean(X(:,1,:)-X(:,2,:)+...
                            X(:,3,:)-X(:,4,:),3)))==1)/size(X,1);
                        yo(ic1)=(sum(sign(mean(X(:,6,:)-1.5,3))==-1) + ...
                            sum(sign(mean(X(:,6,:)-1.5,3))==0)/2)/size(X,1);
                    end
                    %plot(c1,yx,[col(it) '.--']);
                    %plot(c1,yo,[col(it) '.-']);
                    errorbar(c1,yx,yx.*(1-yx)/sqrt(size(X,1)),[col(ip) '.--']);
                    errorbar(c1,yo,yx.*(1-yx)/sqrt(size(X,1)),[col(ip) '.-']);
                end
                if 1 % unattended Condition
                    for ic2=1:length(c1), ic1=1;
                        o=S_Retrieve('get-one',olist,tau(it),c1(ic1),c2(ic2),pL(ip));
                        X=o.X;
                        yx(ic2)=sum(sign(squeeze(mean(X(:,1,:)-X(:,2,:)+...
                            X(:,4,:)-X(:,3,:),3)))==1)/size(X,1);
                        yo(ic2)=(sum(sign(mean(X(:,6,:)-1.5,3))==1) + ...
                            sum(sign(mean(X(:,6,:)-1.5,3))==0)/2)/size(X,1);
                    end
                end
            end
        end
        
    case 'template'
        olist=varargin{1};
        tau=S_Retrieve('tau',olist); ntau=length(tau);
        pL=S_Retrieve('pL',olist);   npL=length(pL);
        c1=S_Retrieve('c',olist,1);  nc1=length(c1);
        c2=S_Retrieve('c',olist,4);  nc2=length(c2);
        fig=Get_Figure(['PNIPS11' fct]);
        %Subplot(3,1,3,1);
        pair=[[1 2];[3 4];[1 3];[2 4];[1 4];[2 3]];
        cor=zeros(ntau,nc1,1,npL,size(pair,1)); sig=zeros(size(cor));
        for it=1:length(tau), for ip=1:length(pL)
                for ic1=1:length(c1), ic2=1; %for ic2=1:length(c2)
                    o=S_Retrieve('get-one',olist,tau(it),c1(ic1),c2(ic2),pL(ip));
                    X=o.X;
                    for i=1:size(pair,1)
                        x=mean(X(:,pair(i,1),:),3)-mean(mean(X(:,pair(i,1),:),1),3);
                        y=mean(X(:,pair(i,2),:),3)-mean(mean(X(:,pair(i,2),:),1),3);
                        [c p]=corrcoef(x,y);
                        cor(it,ic1,ic2,ip,i)=c(1,2);
                        sig(it,ic1,ic2,ip,i)=p(1,2);
                    end
                end
            end, end
        ip=npL;
        cor(end,:,1,npL,1)=zeros(1,nc1); % must be zero for tau=1!
        for ic=1:3 %nc1-1 % exclude highest contrast (not significant)
            Subplot(1); plot(tau,cor(:,ic,1,npL,1),[col(ic) '.-']);
            Subplot(3); plot(tau,sig(:,ic,1,npL,1),[col(ic) '.-']);
        end
        for i=[1 3], Subplot(i); axis tight; end
        
    case 'Psychometry-Old'
        o=varargin{1}; % takes 4D cell array
        [n_tau n_c1 n_c2 n_pL]=size(o);
        for i=1:n_tau, tau(i)=o{i,1,1,1}.P.tau; end, disp(['tau: ' num2str(tau)]);
        for j1=1:n_c1, c1(j1)=o{1,j1,1,1}.I.c(1); end, disp(['c1: ' num2str(c1)]);
        for j2=1:n_c2, c2(j2)=o{1,1,j2,1}.I.c(4); end, disp(['c2: ' num2str(c2)]);
        for k=1:n_pL,  pL(k) =o{1,1,1,k}.S.pL(1); end, disp(['pL: ' num2str(pL)]);
        Get_Figure('SIP');
        Subplot(n_tau,1);
        for i=1:n_tau, Subplot(i); title(['tau: ' num2str(tau(i))]);
            for k=[1 3]% n_pL]
                if 1 % Attended condition
                    for j=1:n_c1, X=o{i,j,1,k}.X;
                        %yx{j}=sum(sign(squeeze(mean(X(:,1,:)-X(:,2,:),3)))==1)/size(X,1);
                        yx{j}=sum(sign(squeeze(mean(X(:,1,:)-X(:,2,:)+...
                            X(:,3,:)-X(:,4,:),3)))==1)/size(X,1);
                        yo{j}=(sum(sign(mean(X(:,6,:)-1.5,3))==-1) + ...
                            sum(sign(mean(X(:,6,:)-1.5,3))==0)/2)/size(X,1);
                        yxm(j)=mean(yx{j});
                        yom(j)=mean(yo{j});
                    end
                    Plot_Points_Error(c1,yx,95,col(k));
                    plot(c1,yxm,[col(k) 'o--']);
                    Plot_Points_Error(c1,yo,95,col(k));
                    plot(c1,yom,[col(k) 'o-']);
                end
                if 1 % Unattended condition
                    for j=1:n_c2, X=o{i,1,j,k}.X;
                        %yx{j}=sum(sign(squeeze(mean(X(:,1,:)-X(:,2,:),3)))==1)/size(X,1);
                        yx{j}=sum(sign(squeeze(mean(X(:,1,:)-X(:,2,:)+...
                            X(:,4,:)-X(:,3,:),3)))==1)/size(X,1);
                        yo{j}=(sum(sign(mean(X(:,6,:)-1.5,3))==1) + ...
                            sum(sign(mean(X(:,6,:)-1.5,3))==0)/2)/size(X,1);
                        yxm(j)=mean(yx{j});
                        yom(j)=mean(yo{j});
                    end
                    Plot_Points_Error(c1,yx,95,col(k));
                    plot(c1,yxm,[col(k) 'x--']);
                    Plot_Points_Error(c1,yo,95,col(k));
                    plot(c1,yom,[col(k) 'x-']);
                end
            end
        end
        
        
    otherwise
        warning(['option ' fct ' not implemented!']);
end