cd /dssg/home/acct-bioxsyy/bioxsyy-zy/toolbox/cobratoolbox
%initialization
initCobraToolbox
solverOK=changeCobraSolver('ibm_cplex','LP');
modPath = [pwd filesep 'AGORA2'];
numWorkers = 32;
taxonLevel='strain';
dietFilePath='chinese_test_title.txt';
abunFilePath='normalizedCoverage 2.csv';
computeProfiles = true;
infoFilePath = '';
[init, netSecretionFluxes,netUptakeFluxes, Y, modelStats, summary] = initMgPipe(modPath, abunFilePath, computeProfiles,'dietFilePath', dietFilePath, 'infoFilePath', infoFilePath, 'numWorkers', numWorkers);
taxInfo = '/dssg/home/acct-bioxsyy/bioxsyy-zy/toolbox/cobratoolbox/AGORA2_infoFile.xlsx';
fluxPath = [pwd filesep 'Results' filesep 'inputDiet_net_secretion_fluxes.csv'];
corrMethod = 'Pearson';
[FluxCorrelations, PValues, TaxonomyInfo] = correlateFluxWithTaxonAbundance(abunFilePath, fluxPath, taxInfo, corrMethod);
save('my_workspace.mat')
