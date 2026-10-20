cd /dssg/home/acct-bioxsyy/bioxsyy-zy/toolbox/cobratoolbox

initCobraToolbox
solverOK=changeCobraSolver('ibm_cplex','LP');
numWorkers = 64;
resPath = [pwd filesep 'Results'];
constrModPath = [resPath filesep 'Diet'];
metList = {'retinol'};
[minFluxes,maxFluxes,fluxSpans] = predictMicrobeContributions(constrModPath, 'metList', metList, 'numWorkers', numWorkers);
