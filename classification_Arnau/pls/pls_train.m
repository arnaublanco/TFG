function model=pls_train(x,labels,options)
% 
% function model=pls_train(x,labels,options)
% 
% PLS: Partial Least Squares or Projection to Latent Structures
% The algorithm finds the Latent Variables in the X-block which maximize 
% the variance that explains the most of the Y-block. NIPALS algorithm [1] is
% used in the current implementation
% 
% INPUT ARGUMENTS
% ================
% % x --> X-block. Matrix of responses. Matrix of features for the independent variables (m). Size n x m.
% labels --> Y-block. Matrix of concentrations. Matrix of features for the dependent variables (p). Size n x p.
%           options.npls --> Number of Latent Variables (LV) to be found. (p by default)
%           options.prep --> Preprocessing options to apply to training data (X-Y):
%                   0: No preprocessing (DEFAULT)
%                   1: Mean substraction (columns)
%                   2: Autoscaling to unit variance (columns)
%          options.cl --> Confidence level for outlier detection. (95% by default)
%
% 
%  X=T*P'+E
%  Y=U*Q'+F
% 
% Note 1: The function "preprocDMA" is used to perform preprocessing steps to
% Xtrain and Ytrain. The same preprocessing steps must be applied to new
% data (Xtest) which will be used for prediction.
%
% OUTPUT ARGUMENTS
% ================
% model --> Output data structure containing the following fields:
%              .reg--> Matrix of regression coefficients. Size m x p. (Y=X*reg)
%              .T --> X-scores matrix. Size n x ncomp
%              .W --> X inner weights matrix. Size m x ncomp
%              .P --> X-loadings matrix. Size m x ncomp
%              .U --> Y-scores matriz. Size n x ncomp
%              .B --> Y inner regression coefficients. Size ncomp x 1
%              .Q --> Y-loadings matrix. Size p x ncomp
%              .ssq --> Residuals sum of squares (%). First column: X-block, 2nd column: Y-block
%                   Each row contains the percentage of residual for X-block or
%                   Y-block after extracting k LV.
%              .prep --> Preprocessing options applied to training data (X-Y):
%                   0: No preprocessing (DEFAULT)
%                   1: Mean substraction (columns)
%                   2: Scale to unit variance (columns)
%              .mux --> a row vector containing the mean value per each
%                     column of the training dataset "Xt"
%              .sx --> a row vector containing the standard deviation per each
%                    column of the training dataset "Xt"
%              .muy --> a row vector containing the mean value per each
%                     column of the training dataset "Yt"
%              .sy --> a row vector containing the standard deviation per each
%                    column of the training dataset "Yt"
%
% 
% Projection of new data:
% ------------------------
% Tnew=Xnew*W*pinv(P'*W);
% 
% References:
% -----------
% [1] Geladi, P. & Kowlaski B. Partial least square regression: A tutorial.
% Analytica Chemica Acta, 35, pp. 1-17. 1986. (NIPALS algorithm)
% [2] Wold, S, Sjöström, M., Eriksson, L. PLS-regression: a basic tool of chemometrics.
% Chemometrics and Intelligent Laboratory Systems, 58, pp. 109–130. 2001.
% [3] Sijmen de Jong, SIMPLS: An alternative approach to partial least squares regression.
% Chemometrics and Intelligent Laboratory Systems, Volume 18, Issue 3,
% March 1993, Pages 251-263. (SIMPLS algorithm)
% 
% 
import classification.pls.*;

if nargin<1
    error('Not enough input arguments');    
end

%Load data from input model
    X=x;
    Y=labels;
if ~isfield(options,'npls') || options.npls == -1        
    npls=size(Y,2);
else    
    npls=options.npls;
end
% if ~isfield(options,'prep')        
%     prep=0;
% else
%     prep=0;
% end
prep=0;
if ~isfield(options,'cl')
    cl=95;
else
    cl=options.cl;
end
    
% %Confidence level for outlier detection    
%     cl=95;
%Tolerance: convergence criterion
    tol=1e-10;
%Original data without preprocessing    
%    X0=x;    
%    Y0=labels;
%Dimensions of the loaded data
    [nx,m]=size(X);
    [ny,py]=size(Y);
%Maximum number of PLS components
    max_npls=min([nx,m]);

if nx~=ny
    error('The number of samples (rows) must be equal in X-block and Y-block');    
end


%=========================================================================
% CROSS VALIDATION
%=========================================================================

% thr=1; %Threshold to select the number of LV (in %)

% %Leave-one-out cross-validation
% %-------------------------------
% disp('Please wait. Performing CV leave-one out...');
%     RMSECV_loo=pls_cv(X0,Y0,1,prep);  
%     aux=mean(RMSECV_loo,1);
% %    npls_loo=find(abs(diff(aux))./aux(1:end-1)*100<thr,1)+1;
%     npls_loo=find(aux==min(aux),1);
%     if isempty(npls_loo)
%         npls_loo=0;
%     end

%Random selection
%-----------------
%thr=1;  %threshold in percentage
%disp('Please wait. Performing CV Random selection (70% train - 30% validation)...');
%    RMSECV_rs=pls_cv(X0,Y0,2,prep);  
%    aux=mean(RMSECV_rs,1);
%    aux2=abs(diff(aux));
%    npls_rs=find(aux2./sum(aux2)*100<thr,1);
%    npls_rs=find(aux==min(aux),1);
%    npls_rs=min([find(min(aux)==aux),npls_rs]);
%    if isempty(npls_rs)
%        npls_rs=0;
%    end


%K-fold cross-validation
%--------------------------
%thr=1;  %threshold in percentage
%disp('Please wait. Performing CV Kfold (K=10)...');
%    RMSECV_kfold=pls_cv(X0,Y0,3,prep,10);
%    aux=mean(RMSECV_kfold,1);
%    aux2=abs(diff(aux));
%    npls_kfold=find(aux2./sum(aux2)*100<thr,1);
%     npls_kfold=find(aux==min(aux),1);    
%    npls_kfold=min([find(min(aux)==aux),npls_kfold]);
%    if isempty(npls_kfold)
%        npls_kfold=0;
%    end

    
%=========================================================================
% CROSS VALIDATION
%=========================================================================

%Important question:
%Number of components selected by the user or suggested by the program?
%     if npls<py
%         npls=py;   %Minimum number of LV = Number of substances
%     end

mux=mean(X,1);
sx=std(X,0,1);
muy=mean(Y,1);
sy=std(Y,0,1);

%Preprocessing
%--------------
if prep>0    
    X=X-repmat(mux,nx,1);    %mean substraction    
    Y=Y-repmat(muy,ny,1);    %mean substraction    
    if prep>1    
        X=X./repmat(sx,nx,1); %scale to unit variance        
        Y=Y./repmat(sy,ny,1); %scale to unit variance        
    end
end

%Preprocessed original data
    X0p=X; Y0p=Y;

%Initialization
%  X=T*P'+E
%  Y=U*Q'+F

T=zeros(nx,max_npls);  %X scores    
W=zeros(m,max_npls);   %Inner X weights    
P=zeros(m,max_npls);   %X loadings    
U=zeros(ny,max_npls);  %Y scores    
B=zeros(max_npls,1);   %Inner Y weights    
Q=zeros(py,max_npls);  %Y loadings
%     reg=zeros(m,py);    %Regression coefficients    
ssq     = zeros(max_npls,2); %Sum of squares    
ssqx    = (sum(X(:).^2)'); %Sum of squares X-block    
ssqy    = (sum(Y(:).^2)'); %Sum of squares Y-block

%Basic NIPALS algorithm [1]
%--------------------------
%It is assumed that X and Y are mean-centered or previously preprocess
       
for k=1:max_npls
  
    %------------------------------------------------------------------------
    % NIPALS for 1 LV starts here
    %------------------------------------------------------------------------
    
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

ssqdif            = zeros(max_npls,2);
ssqdif(1,1)       = 100 - ssq(1,1);
ssqdif(1,2)       = 100 - ssq(1,2);

for k1=2:max_npls
    for k2=1:2        
        ssqdif(k1,k2) = -ssq(k1,k2) + ssq(k1-1,k2);
    end
end
ssq  = [ssqdif(:,1) ssqdif(:,2)];


%=========================================================================
% OSC
%=========================================================================
%Aqui hay que evaluar si es necesario utilizar OSC.
% Si la varianza recuperada por X es alta pero la varianza de Y es baja,
% significa que estamos introduciendo mucha variabilidad que no está
% relacionada con Y y que por lo tanto es ortogonal o próxima a ortogonal.
% OSC permitira eliminar esta variacion en X
% IMPORTANTE: Si OSC se aplica, hay que guardar el modelo a posteriori,
% para que los nuevos datos sean corregidos de la misma manera.

%Si aplicamos OSC , luego hay que volver a hacer el PLS y posteriormente
%seleccionar el numero de LV
%=========================================================================
% OSC
%=========================================================================


%=======================================================================
%Aqui tendriamos que estimar de manera inteligente donde estan los codos y
%cual es el nº optimo de variables latentes
%Deberia haber un campo de salida que fuera la sugerencia del programa
%=======================================================================
%     thr=1; %threshold to select NLV (%)
%     aux=abs(diff(ssq));
%     x_thr=find(aux(:,1)<thr,1);
%     y_thr=find(aux(:,2)<thr,1);
%     pmin_x=find(ssq(:,1)==min(ssq(:,1)));
%     pmin_y=find(ssq(:,2)==min(ssq(:,2)));
% 
%     x_thr=min([x_thr,pmin_x]);
%     y_thr=min([y_thr,pmin_y]);
%     %Number of latent variables suggested by the program
%         npls_sug=max([x_thr,y_thr]);


%=======================================================================
%Aqui hay que sobreescribir las matrices U, Q, P, T, W, B segun el nº de LV
%escogido
% npls a partir de aqui sera el nº de LV escogido o sugerido por el
% programa
%=======================================================================
%Nº variables sugerido por el programa o definido por el usuario?

%Output matrices
T=T(:,1:npls);
W=W(:,1:npls);
P=P(:,1:npls);
U=U(:,1:npls);
B=B(1:npls);
Q=Q(:,1:npls);

%matrix of regression coefficients
%---------------------------------
reg=zeros(npls*py,m);
C = zeros(py*npls,m);

if py == 1
    C = (W*pinv(P'*W)*diag(B))';    
else
    Cp = (W*pinv(P'*W)*diag(B))';    
    for k = 1:npls    
        Cpp = Cp(k,:);        
        C((k-1)*py+1:k*py,:) = diag(Q(:,k))*Cpp(ones(py,1),:);       
    end
end

reg(1:npls*py,:)=C;

if py>1
    for k=2:npls    
        jj       = (k-1)*py+1;        
        i0       = jj-py;        
        reg(jj:k*py,:) = reg(jj:k*py,:) + reg(i0:(k-1)*py,:);        
    end
else
    reg  = cumsum(reg,1);    
end

%Regression coefficients
reg=reg((1:py)+(py*(npls-1)),:)';


%Prediction of the training dataset
%------------------------------------
    Ypred=X0p*reg;

    F=Ypred-Y0p;
    RMSEC=sqrt(sum(F(:).^2)/ny/py);    
    Qy=sum(F.^2,2);
    resYnorm=(F-repmat(mean(F,1),size(F,1),1))./repmat(std(F,0,1),size(F,1),1);
    
%Undo preprocessing applied to Y-block
if prep==1
    Ypred=Ypred+repmat(muy,ny,1);
elseif prep==2
    Ypred=Ypred.*repmat(sy,ny,1)+repmat(muy,ny,1);
end


%===========================
% OUTLIER DETECTION
%===========================

%Residuals Q
%------------
E=X0p-T*P';
Qx=sum(E.^2,2);

sE=svd(E);
if nx>1
    sE=sE.^2/(nx-1);
else
    sE=sE.^2;
end

alfa = (2*cl-100)/100;
O1=sum(sE(:,1));
O2=sum(sE(:,1).^2);
O3=sum(sE(:,1).^3);

if O1 == 0
    Qalfa=NaN;    
else
    h0=1-2*O1*O3/(3*O2^2);    
    if h0<0.001    
        h0=0.001;        
    end
    Calfa = sqrt(2)*erfinv(alfa);    
    %Outlier detection   
    Qalfa=O1*(Calfa*sqrt(2*O2*h0^2)/O1 + 1 + O2*h0*(h0-1)/O1^2) ^(1/h0);    
end

% Q-statistics for a new preprocessed sample will be calculated as:
% Qval = xval*(eye(size(xval,2))-P*P')*xval'

%Hotelling's T^2 statistics
%---------------------------
if nx>1
    f=diag(sqrt(1./(diag(T'*T)/(nx-1))));
else
    f=diag(sqrt(1./(diag(T'*T))));
end

if npls>1
    tsq1=(sum((f*T').^2))';
%    tsq2=(sum((f*P').^2))';
else
    tsq1=((f*T').^2)';
%    tsq2=((f*P').^2)';
end

if nx>1
    leverage=tsq1/(nx-1);
else
    leverage=tsq1;
end

if npls >= nx
    Halfa = NaN;
else
    %Outlier detection
    alpha = (100-cl)/100;
    Halfa = npls*(nx-1)/(nx-npls)*ftest(alpha,npls,nx-npls);
end

if nx>1
    Halfa=Halfa/(nx-1);
end
    
% T^2 Hotellings-statistics for a new preprocessed sample will be calculated as:
% HT2val = xval*P*diag(1./v1(1:npc))*P'*xval'

%===========================
% END OF OUTLIER DETECTION
%===========================

%Output data structure
%======================
    model = struct();
    model.name = 'pls';
    model.reg=reg;
    model.T=T;
    model.W=W;
    model.P=P;
    model.U=U;
    model.B=B;
    model.Q=Q;
    model.ssq=ssq;
    model.prep=prep;
    model.mux=mux;
    model.sx=sx;
    model.muy=muy;
    model.sy=sy;
%     model.npls_sug=npls_sug;
%     model.npls_loo=npls_loo;
%    model.npls_rs=npls_rs;    
%    model.npls_kfold=npls_kfold;     
    model.Ypred=Ypred;
%RMSEC
    model.RMSEC=RMSEC;
%     model.RMSECV_loo=mean(RMSECV_loo,1);
%    model.RMSECV_rs=mean(RMSECV_rs,1);    
%    model.RMSECV_kfold=mean(RMSECV_kfold,1);      
%Confidence level in outlier identification    
    model.cl=cl;
%T^2 Hotelling statistics    
    model.leverage=leverage;
    model.Halfa=Halfa;
%Q residuals statistics    
    model.Qx=Qx;
    model.Qalfa=Qalfa;
    model.Qy=Qy;
    model.resYnorm=resYnorm;
%Outliers detection using Q residuals and T2 statistics
    model.Qout=double(Qx>Qalfa);
    model.levout=double(leverage>Halfa);
    model.outlier=double(Qx>Qalfa & leverage>Halfa);    
    model.resYout=double(resYnorm>3 | resYnorm<-3);

%=======================================    
%Representaciones graficas importantes:
%=======================================
%Ypred vs Y0
%Qx vs leverage con limites Qalfa, Halfa
%resYnorm vs leverage con limites en 3 y Halfa
%Scores con LV1 vs LV2, que se pueda escoger las componentes de cada eje
%Loadings LV1, LV2, etc...
%Varianza recuperada del bloque X y del Y
%Coeficientes de regresion


