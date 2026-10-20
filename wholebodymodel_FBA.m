initCobraToolbox;
solverOK = changeCobraSolver('ibm_cplex', 'LP');

% location
inputFolderPath = '/dssg/home/acct-bioxsyy/bioxsyy-zy/toolbox/cobratoolbox/combine_WB_model_new';
outputFolderPath = '/dssg/home/acct-bioxsyy/bioxsyy-zy/toolbox/cobratoolbox/WB_fba';  
filePrefix = 'microbiota_model_diet_WB_';

% make sure exist
if ~exist(outputFolderPath, 'dir')
    mkdir(outputFolderPath);
end

% check all files
files = dir(fullfile(inputFolderPath, [filePrefix, '*.mat']));

% for files
for k = 1:length(files)
    fileName = files(k).name;
    sampleName = extractBetween(fileName, filePrefix, '.mat');

    % load models
    data = load(fullfile(inputFolderPath, fileName));
    model = data.combinedModel;

    % FBA 
    param = '';
    solution = optimizeWBModel(model, param);

    % extract flux and reaction
    fluxValues = solution.v;
    reactionNames = model.rxns;

    % check length
    if length(fluxValues) ~= length(reactionNames)
        error('reaction and flux number not same');
    end

    % create a table
    T = table(reactionNames, fluxValues, 'VariableNames', {'Reaction', sampleName{1}});

    % write a new table
    csvFileName = fullfile(outputFolderPath, sprintf('reaction_fluxes_%s.csv', sampleName{1}));
    writetable(T, csvFileName);

    % print finish
    fprintf('flux analysis done %s\n', csvFileName);
end