function updatePPP = misdetectionPPP(predictPPP,model)

if model.occlusion
    Pd = predictPPP.Pd;
else
    Pd = model.Pd;
end
Qd = 1-Pd;
w1 = predictPPP.w + log(Qd);
w2 = predictPPP.w + log(Pd) + arrayfun(@(x) x.alpha*log(x.beta/(x.beta+1)), predictPPP.GPState);
updatePPP.w = [w1;w2];

updatePPP.GPState = [predictPPP.GPState;arrayfun(@(x) bplus(x), predictPPP.GPState)];


    function GPState = bplus(GPState)
        GPState.beta = GPState.beta + 1;
    end

end

