function [PPP,MBM] = occlusionCorrection(MBM,PPP,time,filter,model,egopos)
    %Based on MBM, correct detection probability and gamma parameters based
    %on one target occluding another
    nHypotheses = length(MBM.w);
    egoN = egopos.N;
    egoE = egopos.E;
    for nh = 1:nHypotheses
       track_indices = find(MBM.table(nh,:)>0);
       nTar = length(track_indices);
       angles = zeros(nTar,2);
       varAngles = zeros(nTar,2);
       radii = zeros(nTar,1);
       existProb = zeros(nTar,1);
       for iTar = 1:nTar
           currTrack = MBM.track{track_indices(iTar)}(MBM.table(nh,track_indices(iTar)));
           if currTrack.Bern.t_death(end) == time && currTrack.Bern.w_death(end) > 0
               state = currTrack.Bern.GPState(end).state;
               cov = currTrack.Bern.GPState(end).cov;
               state(1) = state(1)-egoN;
               state(2) = state(2)-egoE;
               [angles(iTar,:),varAngles(iTar,:)] = filter.getVisibilityAngle(state,cov);
               radii(iTar,:) = sqrt(state(1).^2+state(2).^2);
               existProb(iTar) = currTrack.Bern.w_death(end);
           end
       end
       correctedAngles = angles;
       for i = 1:nTar
           currTrack = MBM.track{track_indices(i)}(MBM.table(nh,track_indices(i)));
           nOcclusionProbability = 1;
           visibilityRatios = [];
           partialOcclusionProbabilities = [];
           for j = 1:nTar
               distancecond = radii(i)>radii(j)+min(currTrack.Bern.GPState(end).state(7:end));
               if(distancecond && i~=j && existProb(j)>0)
                   probFullOcl = mvncdf((angles(i,1)-angles(j,1))/sqrt(varAngles(j,1)+varAngles(i,1)))*...
                       mvncdf((angles(j,2)-angles(i,2))/sqrt(varAngles(j,2)+varAngles(i,2)));

                   probPartOclMin = mvncdf((angles(i,1)-angles(j,1))/sqrt(varAngles(j,1)+varAngles(i,1)))*...
                       mvncdf((angles(j,2)-angles(i,1))/sqrt(varAngles(j,2)+varAngles(i,1)));

                   probPartOclMax = mvncdf((angles(i,2)-angles(j,1))/sqrt(varAngles(j,1)+varAngles(i,2)))*...
                       mvncdf((angles(j,2)-angles(i,2))/sqrt(varAngles(j,2)+varAngles(i,2)));

                    condmin = (correctedAngles(i,1) >= angles(j,1)) && (angles(j,2) >= correctedAngles(i,1));
                    condmax = (correctedAngles(i,2) >= angles(j,1)) && (angles(j,2) >= correctedAngles(i,2));
                    nOcclusionProbability = (1-(existProb(j)*probFullOcl))*nOcclusionProbability;
                    if(condmin && condmax)
                        visibilityRatios = [visibilityRatios 0];
                        partialOcclusionProbabilities = [partialOcclusionProbabilities existProb(j)*probFullOcl];
                    elseif(condmin)
                        correctedAngles(i,1)=max([angles(j,2) correctedAngles(i,1)]);
                        currentVisibilityRatio = (correctedAngles(i,2)-correctedAngles(i,1))/(angles(i,2)-angles(i,1));
                        visibilityRatios = [visibilityRatios currentVisibilityRatio];
                        partialOcclusionProbabilities = [partialOcclusionProbabilities probPartOclMin*existProb(j)];
                    elseif(condmax)
                        correctedAngles(i,2)=min([angles(j,1) correctedAngles(i,2)]);
                        currentVisibilityRatio = (correctedAngles(i,2)-correctedAngles(i,1))/(angles(i,2)-angles(i,1));
                        visibilityRatios = [visibilityRatios currentVisibilityRatio];
                        partialOcclusionProbabilities = [partialOcclusionProbabilities probPartOclMax*existProb(j)];
                    end
               end
%                 if(currentOcclusionProbability>occlusionProbability)
%                     occlusionProbability = currentOcclusionProbability;
%                     if(currentVisibilityRatio<visibilityRatio)
%                         visibilityRatio = currentVisibilityRatio;
%                     end
%                 end
           end
           if(isempty(visibilityRatios) || min(visibilityRatios) == 0)
                currTrack.Bern.gammaCorrectionFactor =  1;
           else
                corrFactor =  1-max(partialOcclusionProbabilities)*(1-min(visibilityRatios));
                if(corrFactor<0)
                    currTrack.Bern.gammaCorrectionFactor = 0;
                else
                    currTrack.Bern.gammaCorrectionFactor = corrFactor;
                end
           end
           currTrack.Bern.visibleAngles = correctedAngles(i,:);
           currTrack.Bern.Pd = model.Pd*nOcclusionProbability;
           if(currTrack.Bern.Pd == 0)
               currTrack.Bern.Pd = 0.05;
           end
           MBM.track{track_indices(i)}(MBM.table(nh,track_indices(i))) = currTrack;
       end
       [Pd_PPP(nh,:),corrFactor_PPP(nh,:)] = occlusionCorrectionPPP(angles,varAngles,radii,existProb,PPP,filter,model,egopos);
    end
    if(nHypotheses>0)
        PPP.Pd = (exp(MBM.w)'*Pd_PPP)';
        PPP.gammaCorrectionFactor = (exp(MBM.w)'*corrFactor_PPP)';
    else
        nComponents = length(PPP.w);
        %PPP.Pd = model.Pd*ones(nComponents,1);
        PPP.gammaCorrectionFactor = ones(nComponents,1);
    end
end


                    


