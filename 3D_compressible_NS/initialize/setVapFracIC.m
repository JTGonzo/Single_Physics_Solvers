%% initializeVar
% This function initializes the intial vapor fraction of the 
% mutli-phase fluid domain 
function [Sol] = setVapFracIC(Sol, fluid, crd, Flag3D)

% 2D bubble
	rBubble = 0.0004;
	xBubCent = 0.0;
	yBubCent = 0.0;

%% find the domain node indicies of the bubble boundary    
	if(~Flag3D)
		ind1 = find((sqrt( (crd(:,1)-xBubCent).^2 + (crd(:,2)-yBubCent).^2 ) - rBubble) <= 0.0 ) ; % edit - new
		ind2 = find((sqrt( (crd(:,1)-xBubCent).^2 + (crd(:,2)-yBubCent).^2 ) - rBubble) > 0.0 ) ; % edit - new	
	else
		zBubCent = 0.0;
		ind1 = find((sqrt( (crd(:,1)-xBubCent).^2 + (crd(:,2)-yBubCent).^2 + (crd(:,3)-zBubCent).^2 ) - rBubble) <= 0.0 ) ; % edit - new
		ind2 = find((sqrt( (crd(:,1)-xBubCent).^2 + (crd(:,2)-yBubCent).^2 + (crd(:,3)-zBubCent).^2  ) - rBubble) > 0.0 ) ; % edit - new
    end
	% ind1 = find((crd(:,2) <= 0.00011) & (crd(:,2) >= 0.00009) ) ; % edit - new (	% horizontal strip)
	% ind2 = find((crd(:,2) > 0.00011) | (crd(:,2) < 0.00009) ) ; % edit - new	

%% Set the boundary/initial values at the boundary nodes    
    Sol.vapFrac(ind1) = 0.01; % edit - new
	Sol.vapFrac(ind2) = 1.0;
	Sol.p(ind1) = fluid.pVap-50;
	Sol.p(ind2) = fluid.pInf;  % fluid.pres0;
    
    % Sol.p(ind2) = -(rBubble./(sqrt( (crd(ind2,1)-xBubCent).^2 + (crd(ind2,2)-yBubCent).^2 )))*...
					% (fluid.pInf - fluid.pVap) + fluid.pInf;
    
	clear ind1 ind2
end