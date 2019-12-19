function [Zpred] = pcog_eye_GLMAR(P,M,U)
% function [Zpred] = pcog_eye_GLMAR(P,M,U)
% Convolution model with AR

    
% P = M.P;
t = M.mspec.t;
nsess = M.mspec.ns; % Number of sessions
nr = length(M.mspec.vnames);

% Assign parameters
%--------------------------------------------------------------------------
Bsess = P.sess; % Session initial values
BAR = P.AR; % Autoregression coefficient
Breg = P.GLM; % Task factos

if strcmpi(M.mspec.cnv,'Dilation') % Single gamma convolution   
    [C] = pcog_eye_convgamma(U,[exp(P.g1(1:3))],t,'single'); % Create kernel for convolution 
elseif strcmpi(M.mspec.cnv,'null') % No convolution    
    [C] = pcog_eye_convgamma(U,[exp([0.5 0.5 0.5])],t,'single'); % Create dummy kernel for convenience
end

if size(U{1},2) % If some regressors to convolve
    Zpred = sum(C.*repmat(Breg,1,nsess),2);
else Zpred = zeros(size(C,1),1); 
end

% Add AR component
%--------------------------------------------------------------------------
ixc = 1; % Counter for data
for js = 1:nsess
    ixd = [ixc:ixc+size(U{js},1)-1]; % Index for data in session
    ixc = ixc+size(U{js},1); % For housekeeping
    Zpred(ixd(2:end)) = Zpred(ixd(2:end)) + BAR*( M.Z(ixd(1:end-1))-Zpred(ixd(1:end-1))); % AR component: PE
    Zpred(ixd(1)) = Zpred(ixd(1)) + Bsess(js); % fit initial value
end

function [C] = pcog_eye_convgamma(U,par,t,cnv)

 if strcmpi(cnv,'single')
     h = par(1);
     l = par(2);
     d = par(3);
     
     f = spm_Gpdf(t-d,h,l);
     
      for js=1:length(U)
            C{js} = [zeros(size(U{js}))]; % No constant term 
            for jc=1:size(C{js},2)
                tmp = conv(U{js}(:,jc),f);
                C{js}(:,jc)=tmp(1:size(C{js},1),1);
            end
        end
        C = blkdiag(C{:});
 end
