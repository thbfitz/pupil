function [hac,M0] = pcog_eye_GLMAR_Estep(pr,prH,mspec,sY,sZ,Mpre)
% [hac,M0] = pcog_eye_GLMAR_Estep(pr,prH,mspec,sY,sZ)
% Function to batch estimate parameters for a set of subjects, given set of
% priors. (Can be combined with M step to give hierarchical modelling).
%
% Inputs
%--------------------------------------------------------------------------
% pr - priors on parameters
% hpr - priors on noise
% mspec - details of model
% sY - inputs {nsub,nsess}
% sZ - data {nsub}
% Mpre - previous model fit (optional - used for initialising parameters)
% 
% Outputs
%--------------------------------------------------------------------------
% hac - data structure for GLMAR modelling
% M0 - model structure containing priors and other relevant information
%
% TF 09/19

nsub = length(sZ);

for s=1:nsub
    
    %% Subject-specific stuff
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ns = size(sY(s,:),2); % Number of sessions
    Uin = sY(s,:); % Task stimuli
    Z = sZ{s}; % Observed data
    
    %% Set up model
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [M0{s}] = pcog_eye_GLMAR_Estep_setp(pr,prH,ns);  %% Specify priors - repeated to allow for different numbers of sessions
     M0{s}.IS = 'pcog_eye_GLMAR_learn'; % Function to use
     M0{s}.cogfx = mspec.cogfx; % Cog function to use
     M0{s}.ufx = mspec.ufx; % Function for making DM
     M0{s}.Z = Z;
     M0{s}.mspec = mspec;
     M0{s}.mspec.ns = ns; % Number of sessions
     
     %% Fit model and get results
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [M{s}] = spm_mci_vl(M0{s},Uin,Z); % Fit from prior means
    if ~isempty(Mpre) % Also try fitting from previous posteriors
        M20 = M0{s}; 
        M20.P = Mpre{s}.Ep; % Initial values from previous fitting
        [M2] = spm_mci_vl(M20,Uin,Z); % Fit from prior means
        if M2.F>M{s}.F, M{s} = M2; M0{s} = M20; end
    end
    [pve(s,:)] = pcog_eye_GLMAR_pve(M0{s},M{s}.Ep,Uin,Z,length(pr.E.g1)+length(pr.E.g2)+1); % Calculte pve
    [Zpred,U{s}] = feval(M0{s}.IS,M{s}.Ep,M0{s},Uin); % Predicted time series and 
    
    %% Store outputs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sEp(s,:) = M{s}.Ep;
    sCp(s,:,:) = full(M{s}.Cp);
    sEh(s) = M{s}.Eh;
    sCh(s) = M{s}.Ch;
    sF(s) = M{s}.F; % Subject-specific components of model evidence
    hac = struct('sEp',sEp,'sCp',sCp,'sEh',sEh,'sCh',sCh,'sF',sF,'pve',pve);
    hac.M = M;
    hac.U = U;
    hac.M0 = M0;
end

function [M] = pcog_eye_GLMAR_Estep_setp(pr,prH,ns)

%% Specify priors - repeated to allow for different numbers of sessions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pE0 = spm_vec(pr.E)'; % Mean
if ~isnumeric(pr.C) % Covariance
    pC0 = eye(length(pE0)).*spm_vec(pr.C); % Where only diagonal terms specified
else pC0 = pr.C; % Where full matrix specified
end
hE0 = prH.E; % Log noise precision
hC0 = prH.C; % Covariance of log noise precision 

% Add in session-specific priors
pE0 = [pE0 0*ones(1,ns)]; % Session means are not hierachical
pC0 = blkdiag(pC0,eye(ns).*(2*ones(1,ns))); % Covariance over parameter estimates

%% Set up model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M.pE = pE0;
M.pC = pC0;
M.hE = hE0;
M.hC = hC0;
M.P = M.pE; % Initial values
M.prH = prH;
M.pr = pr;
