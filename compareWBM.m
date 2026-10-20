% Initialize COBRA Toolbox
initCobraToolbox;
solverOK = changeCobraSolver('ibm_cplex', 'LP');

% Define reactions for example
Rxn = {
    'Colon_RE1233C';
};

% Print reaction formulas (for debugging)
%printRxnFormula(modelHM, Rxn);

% Read sample names from text files
SCZ_samples = readSampleNames('male_SCZ.txt');
HC_samples = readSampleNames('male_HC.txt');
samples = [SCZ_samples; HC_samples];

% Initialize results container
results = zeros(length(Rxn), length(samples));

% Path to model files
model_path = '/dssg/home/acct-bioxsyy/bioxsyy-zy/toolbox/cobratoolbox/combine_WB_model_new';

% Process each sample
for i = 1:length(samples)
    sample = char(strtrim(samples{i}));
    model_file = fullfile(model_path, ['microbiota_model_diet_WB_', sample, '.mat']);
    if exist(char(model_file), 'file')
        load(model_file, 'combinedModel');
        if exist('combinedModel', 'var')
            model_S = combinedModel;
            sample_results = zeros(length(Rxn), 1);

            for k = 1:length(Rxn)
                model_S = changeObjective(model_S, Rxn{k});
                model_S.osense = -1; % Maximize
                tic;
                sol = optimizeCbModel(model_S);
                toc;

                if isempty(sol.x) || sol.stat == 0
                    sample_results(k) = NaN; % Mark as NaN if optimization failed
                    break;
                else
                    sample_results(k) = sol.x(find(ismember(model_S.rxns, Rxn{k})));
                end
            end
            results(:, i) = sample_results;
        else
            warning(['Variable "combinedModel" not found in file: ', model_file]);
            results(:, i) = NaN;
        end
    else
        warning(['Model file not found for sample: ', char(sample)]);
        results(:, i) = NaN;
    end
end

% Convert results to table and write to CSV
results_table = array2table(results, 'VariableNames', cellstr(samples));
writetable(results_table, '/dssg/home/acct-bioxsyy/bioxsyy-zy/toolbox/cobratoolbox/myWBM_germfree_compare/male_result_Colon_RE1233C.csv');

disp('Results have been written to male_results.csv');

% Helper function to read sample names
function samples = readSampleNames(filename)
    fid = fopen(filename, 'r');
    if fid == -1
        error(['Cannot open file: ', filename]);
    end
    samples = textscan(fid, '%s', 'Delimiter', '\n');
    samples = samples{1};
    fclose(fid);
end