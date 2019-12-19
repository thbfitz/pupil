function [Y,Puv] = pcog_eye_GLMAR_cog_fstone(Puv,M,Y)
% [Y] = pcog_eye_convAR_cog_fstone(Puv,M,Y)
% Function to create design matrix Y for Francesco tone experiments


mode = 'newt'; %'newt'/'grad' 

%%  Transform and assign parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmpi(M.mspec.cogmod,'null') % No learning
    par = [];
elseif strcmpi(M.mspec.cogmod,'Q2') % Simple dual-update (mean and variance) - Gaussian-like
    par.a = sigmoid(Puv.cog(1)); % Learning rate (for precision)
    par.v0 = exp(Puv.cog(2)); % Initial value (for precision)
  %  par.a2 = sigmoid(Puv.cog(3)); % Learning rate (for mean)
    par.e = exp(Puv.cog(3)); % Learning rate (for mean)
    par.m0 = Puv.cog(4); % Initial value (for mean)
elseif strcmpi(M.mspec.cogmod,'Q') % Simple dual-update (mean and variance)
    par.a = sigmoid(Puv.cog(1)); % Learning rate (for precision)
    par.v0 = exp(Puv.cog(2)); % Initial value (for precision)
    par.e = sigmoid(Puv.cog(3)); % Learning rate (for mean)
    par.m0 = Puv.cog(4); % Initial value (for mean)
elseif strcmpi(M.mspec.cogmod,'PC')
    par.tk = exp(Puv.cog(1)); % Learning rate (for precision)
    par.t0 = exp(Puv.cog(2)); % Initial value (for precision)
    par.tm = exp(Puv.cog(3)); % Learning rate (for mean)
    par.m0 = Puv.cog(4); % Initial value (for mean)
end

%%  Calculate key quantities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for js=1:M.mspec.ns % Loop over sessions
    
    % Model-independent (stimulus bound)
    %----------------------------------------------------------------------
    Y{js}.u.ns = Y{js}.Y.ns;
    Y{js}.u.vnames = M.mspec.vnames;
    Y{js}.u.tx = Y{js}.Y.tx;
    tones = Y{js}.Y.y;
    Y{js}.u.targ = ~~Y{js}.Y.lat;
    try, Y{js}.u.ishigh = Y{js}.Y.ishigh; end
    Y{js}.u.prb = (Y{js}.Y.y==2000)+(Y{js}.Y.y==500);
    Y{js}.u.txp = Y{js}.u.tx(find(Y{js}.u.prb)); % probe tones
    Y{js}.u.txn=Y{js}.u.tx(find(~Y{js}.u.prb)); % Other tones
    
    % Learning
    %----------------------------------------------------------------------
    if strcmpi(M.mspec.cogmod,'Q2')||strcmpi(M.mspec.cogmod,'Q')
        % Dynamic belief updating
        m=zeros(1,length(tones)); m(1) = par.m0;
        v = zeros(1,length(tones)); v(1) = par.v0;
        for jt=2:length(tones)
            if strcmpi(M.mspec.cogmod,'Q2')
                m(jt) = ((log(tones(jt-1))/v(jt-1)) + (m(jt-1)/par.e))/(1/v(jt-1)+1/par.e); % Mean
            elseif strcmpi(M.mspec.cogmod,'Q')
                m(jt) = m(jt-1) + par.e*(log(tones(jt))-m(jt-1)); % Mean
            end
            sqe=(log(tones(jt))-m(jt)).^2;
            v(jt) = v(jt-1) + par.a*(sqe-v(jt-1)); % Variance

        end
    elseif strcmpi(M.mspec.cogmod,'PC')
        % Dynamic belief updating
        m=zeros(1,length(tones));
        k = zeros(1,length(tones)); 
        
        for jt=1:length(tones)
            % initialise
            if jt==1
                mpr(jt) = par.m0;
                kpr(jt) = log(par.t0);
                tm(jt) = par.tm;
                tk(jt) = par.tk;
            end
            m(jt) = mpr(jt); 
            k(jt) = kpr(jt); 
            for jit=1:16 %00
                
                % Mean terms
                %----------------------------------------------------------
                pe = log(tones(jt))-m(jt);
                dfm = exp(k(jt))*(pe) - (m(jt)-mpr(jt))*tm(jt);
                dfm2 = -exp(k(jt))-tm(jt);
                if strcmp(mode,'newt')
                    m(jt) = m(jt) - dfm/dfm2;
                elseif strcmp(mode,'grad')
                     m(jt) = m(jt) + dfm/16; % Rate is arbitrary (?)
                end
                
                % Variance/precision terms
                %----------------------------------------------------------
                pe = log(tones(jt))-m(jt); % recalculate with new mean
                dfk = 0.5 - 0.5*exp(k(jt))*pe^2 - (k(jt)-kpr(jt))*tk(jt);
                dfk2 = -exp(k(jt))*pe^2 - tk(jt);
                if strcmp(mode,'newt')
                     k(jt) = k(jt) - dfk/dfk2;
                elseif strcmp(mode,'grad')
                     k(jt) = k(jt) + dfk/16; %fppm; % Rate is arbitrary (?)
                end

                % Stopping criterion?
                %----------------------------------------------------------
%                 if jit>1
%                     if strcmp(mode,'newt')
%                       if ((abs(dfk/dfk2)<k(jt)/1e4))&&((abs(dfm/dfm2))<m(jt)/1e4) % If change by less than 0.1%
%                           break
%                       end
%                     end
%                 end
                % tmp(jit+1,:) = [m(jt) k(jt)];
                
                %??? Calculate variational free energy
         %       F(jit+1) = log(spm_Npdf(log(tones(jt)),m(jt),exp(kt(jt))));% - spm_kl_normal(m(jt),exp(kt(jt)),m(jt-1),exp(kt(jt-1)))
                % tmp4(jit,:) = [m(jt) m2(jt) k(jt) k2(jt)];
%                 if jit>1
%                     if max([abs(tmp4(jit,2)-tmp4(jit-1,2)) abs(tmp4(jit,4)-tmp4(jit-1,4))])<0.000001, 
%                         break, end
%                 end
               % F2(jit+1) = log(spm_Npdf(log(tones(jt)),m2(jt),exp(kt2(jt)))) - spm_kl_normal(m2(jt),exp(kt2(jt)),m(jt-1),par.vm);
                   
            end
            
            % Priors for next trial
            %--------------------------------------------------------------
            if jt<length(tones)
                mpr(jt+1) = m(jt);
                kpr(jt+1) = k(jt);
               if M.mspec.fix(1) % If estimates of variance are fixed
                    tk(jt+1) = par.tk; %
               else tk(jt+1) = 1/(1/par.tk - 1/dfk2); %
               end
               if M.mspec.fix(2) % If estimates of mean are fixed
                    tm(jt+1) = par.tm; %Check this
               else tm(jt+1) = 1/(1/par.tm - 1/dfm2); %Check this
               end
                % If want to used fixed undertainty about predictions
%                 vm(jt+1) = par.vm;
%                 vk(jt+1) = par.vv; 
            end  
        end
        tau=exp(k);
    elseif strcmpi(M.mspec.cogmod,'null') % Using true values
        v = 0.5 + Y{js}.u.ishigh;
        tau=1./v;
        m = log(500)*ones(size(Y{js}.u.ishigh));
        
    end
    
    try
        Y{js}.u.s=-log(spm_Npdf(log(tones),m,1./(tau+eps))); % Surprise
        % Store data
        Y{js}.u.m=m;
        Y{js}.u.tau=tau;
        Y{js}.u.k=k;
    end
    try
        Y{js}.u.mpr=mpr;
        Y{js}.u.kpr=kpr;
        Y{js}.u.tm=vm;
        Y{js}.u.tk=vk;
    end

    % Continuity across sessions
    %----------------------------------------------------------------------
    if M.mspec.mvcont&&(~strcmp(M.mspec.cogmod,'null'))
        par.m0 = m(end);
        par.t0 = tau(end);
    end

    Y{js}.u.tones=tones;
end
