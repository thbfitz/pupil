function [U,par] = pcog_eye_GLMAR_DM_fstone(Uin,mspec,sessn)
% [U,par] = pcog_eye_GLMAR_DM_fstone(Uin,mspec)
% Function to make design matrix for Francesco experiment 1 (Tones task)
% 
% Inputs
%--------------------------------------------------------------------------
% Uin - stimuli
% mspec - model details
% sessn - session number 
%
% Outputs
%--------------------------------------------------------------------------
% U - design matrix
% par - parameters (permits continuity across sessions)
%
% TF 09/19

%% Unpack inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fn = fieldnames(Uin.u);
for j=1:length(fn)
    eval([fn{j} ' = Uin.u.' fn{j} ';'])
end

%% Create design matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U = zeros(ns,length(vnames));

for jv=1:length(vnames)
    
    % Separate probe and non-probe tones (optionally)
    %----------------------------------------------------------------------
    usex = 1:length(tx); 
    if strcmp(vnames{jv}(end:end-1),'_p')
        usex = find(prb);
    elseif strcmp(vnames{jv}(end:end-2),'_np')
        usex = find(~prb);
    elseif strcmp(vnames{jv}(end:end-2),'_prb')
        usex = find(~prb);
    end
    uset = tx(usex);
    
    % Stimulus-bound regressors
    %----------------------------------------------------------------------
    if strcmp(vnames{jv},'Event')||strcmp(vnames{jv},'Event_p')||strcmp(vnames{jv},'Event_np')
        U(uset,jv)=1;
    elseif strcmp(vnames{jv},'Deviance')||strcmp(vnames{jv},'Deviance_p')||strcmp(vnames{jv},'Deviance_np')
        U(uset,jv)=abs(log(tones(usex))-log(500)); % 'Surprise'
    elseif strcmp(vnames{jv},'Pitch')||strcmp(vnames{jv},'Pitch_p')||strcmp(vnames{jv},'Pitch_np')
        U(uset,jv)=log(tones(usex))-log(500); % Pitch
    elseif strcmp(vnames{jv},'Target')
        U(uset,jv)=targ; % target
    elseif strcmp(vnames{jv},'Drift')||strcmp(vnames{jv},'Drift_p')||strcmp(vnames{jv},'Drift_np')
        U(uset,jv)=1/length(uset):1/length(uset):1; % target
    elseif strcmp(vnames{jv},'ExpDrift')||strcmp(vnames{jv},'ExpDrift_p')||strcmp(vnames{jv},'ExpDrift_np')
        if sessn==1, U(uset,jv)=exp(-[0:length(uset)-1]/20); end % 'novelty'
    elseif strcmp(vnames{jv},'Adapt')||strcmp(vnames{jv},'Adapt_p')||strcmp(vnames{jv},'Adapt_np')
        tmp = [0 abs(diff(log(tones)))];
        U(uset,jv)=tmp(usex); % target
        
%     % If magically know true task conditions
%     %------------------------------------------------------------------
%     elseif strcmp(vnames{jv}(1:min([4 length(vnames{jv})])),'TVar') 
%          if strcmp(vnames{jv},'TVar')||strcmp(vnames{jv},'TVar_p')||strcmp(vnames{jv},'TVar_np')
%             U(uset,jv)=0.5 + ishigh;
%         elseif strcmp(vnames{jv},'TVar*Deviance')||strcmp(vnames{jv},'TVar*Deviance_p')||strcmp(vnames{jv},'TVar*Deviance_np')
%              U(uset,jv)=abs(log(tones(usex))-log(500)).*(0.5 + ishigh);
%         elseif strcmp(vnames{jv},'TVarSurprise')||strcmp(vnames{jv},'TVarSurprise_p')||strcmp(vnames{jv},'TVarSurprise_np')  
%             U(uset,jv)=-log(spm_Npdf(log(tones(usex)),log(500),(0.5 + ishigh)+eps));
%         % Bimodal stuff
%         elseif strcmp(vnames{jv},'TVarU_Surprise')
%             U(uset,jv)=-log(spm_Npdf(log(tones(usex)),mean(log(tones(usex))),var(log(tones(usex)))));
%         elseif strcmp(vnames{jv},'TVarU*Deviance')
%              U(uset,jv)=abs(log(tones(usex))-mean(log(tones(usex))));
%         elseif strcmp(vnames{jv},'TVarB_Surprise')
%             U(uset,jv)=-log(  0.5*spm_Npdf(log(tones(usex)),5,0.5) + 0.5*spm_Npdf(log(tones(usex)),7.5,0.5)); %TESTING
%         elseif strcmp(vnames{jv},'TVarB*Deviance')
%          end
        
    % learn mean and precision
    %------------------------------------------------------------------
    elseif strcmp(vnames{jv}(1:min([3 length(vnames{jv})])),'Prc')%strcmp(vnames{jv}(1:min([3 length(vnames{jv})])),'Var')%||strcmp(vnames{jv}(1:min([4 length(vnames{jv})])),'Var')%||strcmp(vnames{jv},'QVar*Deviance')||strcmp(vnames{jv},'QVarSurprise') 


%         if strcmp(vnames{jv},'Var')||strcmp(vnames{jv},'Var_p')||strcmp(vnames{jv},'Var_np')
%             U(uset,jv)=v(usex);
%         elseif strcmp(vnames{jv},'Var*Deviance')||strcmp(vnames{jv},'Var*Deviance_p')||strcmp(vnames{jv},'Var*Deviance_np')
%              U(uset,jv)=abs(log(tones)-m(usex)).*v(usex);
%         elseif strcmp(vnames{jv},'VarSurprise')||strcmp(vnames{jv},'VarSurprise_p')||strcmp(vnames{jv},'VarSurprise_np')   
%             U(uset,jv)=s(usex); %-log(spm_Npdf(log(tones(usex)),m(usex),v(usex)+eps)+eps);
%         elseif strcmp(vnames{jv},'VarLong')||strcmp(vnames{jv},'VarLong_p')||strcmp(vnames{jv},'VarLong_np')
%             for jx=1:length(uset)-1
%                 U(uset(jx):uset(jx+1)-1,jv)=v(usex(jx));
%             end
%             U(tx(end):min([tx(end)+mspec.SR*1 size(U(:,jv),1)]),jv)=v(end); % One second response after last trial
        if strcmp(vnames{jv},'Prc')||strcmp(vnames{jv},'Prc_p')||strcmp(vnames{jv},'Prc_np')
             U(uset,jv)=tau(usex);
        elseif strcmp(vnames{jv},'Prc*Deviance')||strcmp(vnames{jv},'Prc*Deviance_p')||strcmp(vnames{jv},'Prc*Deviance_np')
             U(uset,jv)=abs(log(tones)-m(usex)).*tau(usex);
        elseif strcmp(vnames{jv},'PrcSurprise')||strcmp(vnames{jv},'PrcSurprise_p')||strcmp(vnames{jv},'PrcSurprise_np')   
            U(uset,jv)=s(usex); %-log(spm_Npdf(log(tones(usex)),m(usex),v(usex)+eps)+eps);
        elseif strcmp(vnames{jv},'PrcLong')||strcmp(vnames{jv},'PrcLong_p')||strcmp(vnames{jv},'PrcLong_np')
            for jx=1:length(uset)-1
                U(uset(jx):uset(jx+1)-1,jv)=tau(usex(jx));
            end
            U(tx(end):min([tx(end)+mspec.SR*1 size(U(:,jv),1)]),jv)=tau(end); % One second response after last trial
        elseif strcmp(vnames{jv},'PrcLog')||strcmp(vnames{jv},'PrcLog_p')||strcmp(vnames{jv},'PrcLog_np')
             U(uset,jv)=k(usex); %log(tau(usex)+eps);
        end  

%         % For continuity across sesions
%         %--------------------------------------------------------------
%         vend=tmpv(end);
%         par.v0 = tmpv(end); % Better but not used
%         par.m0 = mu(end);

    end
    
    % Do normalisation if required
    if mspec.unorm
        if length(unique(U(uset,jv)))>2 % If a binary regressor
            U(uset,jv) = (U(uset,jv)-mean(U(uset,jv)))./var(U(uset,jv));
        end
    end
    
end
% %     % Do normalisation if required
%     if mspec.unorm
%        U = U./var(U);
%     end

