function [pve] = pcog_eye_GLMAR_pve(M0,Ep,y,z,arx)
% [pve] = pcog_eye_GLMAR_pve(M0,Ep,y,z,arx)
% Calculate percentage variance explained with and without AR component

[Zp,U] = feval(M0.IS,Ep,M0,y); % Full prediction

% Simple Percentage variance explained
%--------------------------------------------------------------------------
pve(1) = 1-var(z-Zp)/var(z); 

% Percentage variance explained excluding AR component
%--------------------------------------------------------------------------
ePx = Ep; ePx(arx) = 0; % Remove AR component (set to zero)
[Zpw] = feval(M0.IS,ePx,M0,y); % Prediction without AR
Zpa = Zp-Zpw; % Identify AR predictions
Znw = z-Zpa; % Regress AR out of data
pve(2) = (var(Znw)-var(Znw-Zpw))/var(Znw); % More meaningful pve