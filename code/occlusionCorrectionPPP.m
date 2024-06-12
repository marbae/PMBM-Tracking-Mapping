function [Pd_PPP,corrFactor_PPP] = occlusionCorrectionPPP(angles,varAngles,radii,existProb,PPP,filter,model,egopos)

nPPP = length(PPP.w);
corrFactor_PPP = zeros(1,nPPP);
Pd_PPP = zeros(1,nPPP);
PPPangles = zeros(nPPP,2);
PPPvarangles = zeros(nPPP,2);
PPPradii = zeros(nPPP,1);
egoN = egopos.N;
egoE = egopos.E;
for iTar = 1:nPPP
   currComponent = PPP.GPState(iTar);
   state = currComponent.state;
   cov = currComponent.cov;
   state(1) = state(1)-egoN;
   state(2) = state(2)-egoE;
   [PPPangles(iTar,:),PPPvarangles(iTar,:)] = filter.getVisibilityAngle(state,cov);
    PPPradii(iTar,:) = sqrt(state(1).^2+state(2).^2);
end
nOccl = size(angles,1);
correctedAngles = PPPangles;
for i = 1:nPPP
   currComponent = PPP.GPState(iTar);
   nOcclusionProbability = 1;
   visibilityRatios = [];
   partialOcclusionProbabilities = [];
   for j = 1:nOccl
       distancecond = PPPradii(i)>radii(j)+min(currComponent.state(7:end));
       if(distancecond && existProb(j)>0 && angleDistance(angles(j,1),PPPangles(i,1))<pi/6)
           probFullOcl = mvncdf((PPPangles(i,1)-angles(j,1))/sqrt(varAngles(j,1)+PPPvarangles(i,1)))*...
               mvncdf((angles(j,2)-PPPangles(i,2))/sqrt(varAngles(j,2)+PPPvarangles(i,2)));

           probPartOclMin = mvncdf((PPPangles(i,1)-angles(j,1))/sqrt(varAngles(j,1)+PPPvarangles(i,1)))*...
               mvncdf((angles(j,2)-PPPangles(i,1))/sqrt(varAngles(j,2)+PPPvarangles(i,1)));

           probPartOclMax = mvncdf((PPPangles(i,2)-angles(j,1))/sqrt(varAngles(j,1)+PPPvarangles(i,2)))*...
               mvncdf((angles(j,2)-PPPangles(i,2))/sqrt(varAngles(j,2)+PPPvarangles(i,2)));

            condmin = (PPPangles(i,1) >= angles(j,1)) && (angles(j,2) >= PPPangles(i,1));
            condmax = (PPPangles(i,2) >= angles(j,1)) && (angles(j,2) >= PPPangles(i,2));
            nOcclusionProbability = (1-(existProb(j)*probFullOcl))*nOcclusionProbability;
            if(condmin && condmax)
                visibilityRatios = [visibilityRatios 0];
                partialOcclusionProbabilities = [partialOcclusionProbabilities existProb(j)*probFullOcl];
            elseif(condmin)
                correctedAngles(i,1)=max([angles(j,2) correctedAngles(i,1)]);
                currentVisibilityRatio = (correctedAngles(i,2)-correctedAngles(i,1))/(PPPangles(i,2)-PPPangles(i,1));
                visibilityRatios = [visibilityRatios currentVisibilityRatio];
                partialOcclusionProbabilities = [partialOcclusionProbabilities probPartOclMin*existProb(j)];
            elseif(condmax)
                correctedAngles(i,2)=min([angles(j,1) correctedAngles(i,2)]);
                currentVisibilityRatio = (correctedAngles(i,2)-correctedAngles(i,1))/(PPPangles(i,2)-PPPangles(i,1));
                visibilityRatios = [visibilityRatios currentVisibilityRatio];
                partialOcclusionProbabilities = [partialOcclusionProbabilities probPartOclMax*existProb(j)];
            end
       end
   end
  if(isempty(visibilityRatios) || min(visibilityRatios) == 0)
        corrFactor_PPP(i) =  1;
   else
        corrFactor_PPP(i) =  1-max(partialOcclusionProbabilities)*(1-min(visibilityRatios));
        if(corrFactor_PPP(i)<0)
            corrFactor_PPP(i) = 0;
        end
  end
   Pd_PPP(i) = PPP.Pd(i)*nOcclusionProbability;
   if(Pd_PPP(i) == 0)
       Pd_PPP(i) = 0.05;
   end
end
