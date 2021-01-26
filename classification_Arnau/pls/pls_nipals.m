function [T,W,P,U,B,Q,ssq,reg]=pls_nipals(X,Y,ncomp)
% 
% function [T,W,P,U,B,Q,ssq,reg]=pls_nipals(X,Y,ncomp)
% 
% The function calculates the PLS using the NIPALS algorithm
% 
% INPUT ARGUMENTS
% ===============
% X--> Xblock, where each row corresponds to a sensor response. Size n x m
% Y--> Yblock, where each column corresponds to a magintude of interest
% (e.g: concentration) and each row is related to a sample. Size n x p
% ncomp--> Number of components to calculate the PLS model. If ncomp is not
% specified, ncomp is taken as the minimum between n and m, i.e min([n,m])
% 
% OUTPUT ARGUMENTS
% ================
% 
% T --> X-scores matrix. Size n x ncomp
% W --> X inner weights matrix. Size m x ncomp
% P --> X-loadings matrix. Size m x ncomp
% U --> Y-scores matriz. Size n x ncomp
% B --> Y inner regression coefficients. Size ncomp x 1
% Q --> Y-loadings matrix. Size p x ncomp
% ssq --> Residuals sum of squares (%). First column: X-block, 2nd column: Y-block
% reg --> Regression coefficients for "ncomp" LV. Size m x p
% 
import classification.pls.*;

if nargin<2
    error('Not enough input arguments');
end

[nx,m]=size(X);
[ny,py]=size(Y);
if nx~=ny
    error('The number of samples (rows) must be equal in X-block and Y-block');    
end

if nargin<3
    ncomp=min([nx,m]);
end

%Tolerance: convergence criterion
    tol=1e-10;

%Matrices initialization
T=zeros(nx,ncomp);  %X scores
W=zeros(m,ncomp);   %Inner X weights
P=zeros(m,ncomp);   %X loadings
U=zeros(ny,ncomp);  %Y scores
B=zeros(ncomp,1);   %Inner Y weights
Q=zeros(py,ncomp);  %Y loadings
ssq     = zeros(ncomp,2); %Sum of squares
ssqx    = (sum(X(:).^2)'); %Sum of squares X-block
ssqy    = (sum(Y(:).^2)'); %Sum of squares Y-block
   
for k=1:ncomp

%------------------------------------------------------------------------
% NIPALS for 1 LV starts here
%------------------------------------------------------------------------

% disp(['Calculating PLS NIPALS component ',num2str(k),'/',num2str(ncomp)]);

    %Flag for convergence
        final=0;
    %Score vector initialization    
        told=X(:,1); %Creo que la clave de que funcione esta aquí. 
        %told se debe inicializar a una columna de X, en lugar de todo a
        %ceros

    %1)The column of Y with the greatest variance is chosen
        if py>1   
            u0=Y(:,find(var(Y)==max(var(Y)),1));        
        else
            u0=Y(:,1);        
        end
        u=u0;
    
    while final==0
        
        %2)In the X-block
            w = (u'*X)';
            
        %3) X weights normalization
            w = (w'/norm(w'))';
            w=w(:);
        %4)
            t = X*w;
            t = t(:);
            
        if py == 1          
            q = 1;
            break; %Leave from while
        end

        %5)In the Y-block
            q = (t'*Y)';
        %6)y weights normalization
            q = (q'/norm(q'))';
            q = q(:);

        %7)
            u = Y*q;
            u = u(:);

        %8) Convergence criterion
        if sqrt(sum((told-t).^2))/sqrt(sum(t.^2)) < tol
            final=1; %convergence achieved
        end
        told=t;

    end
    U(:,k)=u;
    Q(:,k)=q;

    %9) Calculate the X loadings
        p = (t'*X/(t'*t))';
        p = p(:);
        
    p_norm=norm(p);        
    %10)X Loadings normalization
        p = p/p_norm;
        P(:,k)=p;
    %11)X-scores
        t = t*p_norm;
        T(:,k)=t;
    %12)X inner weights

        w = w*p_norm;
        W(:,k)=w;
%-----------------------------------------------------------------
%End of NIPALS for 1 LV
%-----------------------------------------------------------------

    %13)Find the regression coefficient b
        b=u'*t/(t'*t);
        B(k)=b;
   
    %Calculate residuals
        resX=X-t*p';
        resY=Y-b*t*q';
    %Deflate X and Y blocks
        X=resX;
        Y=resY;
    %Sum of squares
        ssq(k,1) = (sum(resX(:).^2)')*100/ssqx;
        ssq(k,2) = (sum(resY(:).^2)')*100/ssqy;

end

if nargout > 6
    ssqdif            = zeros(ncomp,2);
    ssqdif(1,1)       = 100 - ssq(1,1);
    ssqdif(1,2)       = 100 - ssq(1,2);
    for k1=2:ncomp
        for k2=1:2        
            ssqdif(k1,k2)   = -ssq(k1,k2) + ssq(k1-1,k2);
        end
    end
    ssq  = [ssqdif(:,1) ssqdif(:,2)];
end

%matrix of regression coefficients
%---------------------------------
if nargout > 7
    
    reg=zeros(ncomp*py,m);
    C = zeros(py*ncomp,m);

    if py == 1
        C = (W*pinv(P'*W)*diag(B))';    
    else
        Cp = (W*pinv(P'*W)*diag(B))';    
        for k = 1:ncomp    
            Cpp = Cp(k,:);        
            C((k-1)*py+1:k*py,:) = diag(Q(:,k))*Cpp(ones(py,1),:);       
        end
    end

    reg(1:ncomp*py,:)=C;

    if py>1
        for k=2:ncomp    
            jj       = (k-1)*py+1;        
            i0       = jj-py;        
            reg(jj:k*py,:) = reg(jj:k*py,:) + reg(i0:(k-1)*py,:);        
        end
    else
        reg  = cumsum(reg,1);    
    end

    %Regression coefficients
        reg=reg((1:py)+(py*(ncomp-1)),:)';

end
