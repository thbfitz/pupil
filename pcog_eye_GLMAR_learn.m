function [Zpred,U,UDM] = pcog_eye_GLMAR_learn(P,M,U)
% [Zpred,U,UDM] = pcog_eye_GLMAR_learn(P,M,U)
% Augment GLMAR model to allow for cognitive modelling
%
%  Inputs
%--------------------------------------------------------------------------
% P - parameters
% M - Model structure (see spm_mci_vl)
% U - Inputs
%
% Outputs
%--------------------------------------------------------------------------
% Zpred - predicted time series
% U - Inputs (augmented with information from cognitive model)
% UDM - Design matrix
%
% TF 09/19

% Label parameters
%--------------------------------------------------------------------------
Puv = spm_unvec(P(1:end-M.mspec.ns+1),M.pr.E); % Parameters for convolution 
Puv.sess=(P(end-M.mspec.ns+1:end)); % Session-specific parameters

% Calculate inputs using cog modelling function
%--------------------------------------------------------------------------
[U] = eval([M.cogfx '(Puv,M,U)']);

% Make design matrix
%--------------------------------------------------------------------------
for js=1:M.mspec.ns
    [UDM{js}] = eval([M.ufx '(U{js},M.mspec,js)']);
end

% Predict data
%--------------------------------------------------------------------------
[Zpred] = pcog_eye_GLMAR(Puv,M,UDM); 
