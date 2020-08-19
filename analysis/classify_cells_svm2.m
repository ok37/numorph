function [ct,feature_table, ftable] = classify_cells_svm2(centroids, path_table, config)
%--------------------------------------------------------------------------
% Classify cell-types using a trained SVM classifier
%--------------------------------------------------------------------------

s = config.patch_size;
if length(s) < 2
    s = 6;
else
    s = s(2);
end

outlierTypes = 4:10;
centroids(:,1:3) = centroids(:,1:3)+1;
include_annotation = true;

% Read classifications + patch feautre info
itable = readtable(fullfile(config.output_directory,'classifier',...
    sprintf('%s_patch_info.csv',config.sample_name)));
ctable = readtable(fullfile(config.output_directory,'classifier',...
    sprintf('%s_classifications.csv',config.sample_name)));
ftable = readtable(fullfile(config.output_directory,'classifier',...
    sprintf('%s_patch_features.csv',config.sample_name)));


% Check table sizes
assert(size(ctable,2) == 1,"Classification table contains more than 1 column")
n_classified = height(ctable);
fprintf("Using %d cells for training SVM model \n",n_classified)
assert(height(ftable) == n_classified, "Number of rows in the feature table do not match")
assert(height(itable) == n_classified, "Number of rows in the info table do not match")

% Combine manual classifications with patch features
ctable.Properties.VariableNames = {'Type'};

if include_annotation
    catLayer = categorical("L"+string(ftable.Layer));
    orderedLayer = categories(catLayer);
    dummyLayer = dummyvar(catLayer);
    tblLayer = array2table(dummyLayer,'VariableNames',orderedLayer);
    ftable = [ftable tblLayer(:,1:end)];

    catStruct = categorical("S"+string(ftable.Structure));
    orderedStruct = categories(catStruct);
    dummyStruct = dummyvar(catStruct);
    tblStruct = array2table(dummyStruct,'VariableNames',orderedStruct);
    ftable = [ftable tblStruct(:,1:end)]; 
end

ftable.Layer = [];
ftable.Structure = [];
trainingData = horzcat(ctable,ftable);

% Train the model
trainedClassifier = trainClassifierSVM(trainingData, outlierTypes);


% Classify centroids in centroid list
ct = zeros(size(centroids,1),1);

% Remove low intensity cells
k_idx = zeros(size(centroids,1),1);
for i = 1:length(config.markers)-1
   idx = i + 5;
   thresh = prctile(centroids(:,idx),config.low_thresh*100);
   k_idx = k_idx | centroids(:,idx)>thresh;
end
ct(~k_idx) = 1;

% Get layer and structure info
structures = bin_annotation_structures(centroids(:,4),'cortex_large');
layers = bin_annotation_structures(centroids(:,4),'layers');

catLayer = categorical("L"+string(layers));
orderedLayer = categories(catLayer);
dummyLayer = dummyvar(catLayer);
tblLayer = array2table(dummyLayer,'VariableNames',orderedLayer);

catLayer = categorical("S"+string(structures));
orderedLayer = categories(catLayer);
dummyLayer = dummyvar(catLayer);
tblLayer2 = array2table(dummyLayer,'VariableNames',orderedLayer);
ls_table = [tblLayer(:,1:end) tblLayer2(:,1:end)];

%ls_table = table(structures,layers,'VariableNames',{'Layer','Structure'});

% Get patch features for remaining centroids
z_pos = unique(centroids(:,3));
pos(:,1) = centroids(:,1)-s;
pos(:,2) = centroids(:,1)+s;
pos(:,3) = centroids(:,2)-s;
pos(:,4) = centroids(:,2)+s;

% Generate vector of column names
feature_table = cell(1,length(config.markers));
feature_names = ["Max","Mean","Std","Middle","MoC","Solidity","FilledArea","CenDist"];
for i = 1:length(config.markers)
    col_names{i} = feature_names + sprintf('_%s',config.markers(i));
end
col_names = [col_names{:}];

z = centroids(:,3);
for i = 1:length(z_pos)
    fprintf("Working on z position %d \n",z_pos(i))
    idx = find(z == z_pos(i) & k_idx);
    if isempty(idx)
        continue
    end
    
    cen_sub = centroids(idx,:);
    p = pos(idx,:);
    
    n_patches = size(cen_sub,1);
    patches = cell(n_patches,length(config.markers));
    files = path_table(path_table.z == z_pos(i),:).file;
    
    % Load images and measure features
    parfor j = 1:length(config.markers)
        img = loadtiff(files{j});
        for k = 1:n_patches
            patch = img(p(k,1):p(k,2),p(k,3):p(k,4));
            patches{k,j} = measure_patch_features(patch, s, false);
        end
    end

    feature_table = array2table(cell2mat(patches),'VariableNames',col_names);
    feature_table = horzcat(ls_table(idx,:),feature_table);
    ct(idx) = trainedClassifier.predictFcn(feature_table);    
end

% Calculate sums
ct(ct>3) = 4;
max_ct = unique(ct(ct>0));
sums = zeros(1,length(max_ct));
for i = 1:length(max_ct)
    sums(i) = sum(ct==max_ct(i));
end
pct = sums/sum(sums(:));
disp(pct)

end


function [trainedClassifier, validationAccuracy] = trainClassifierSVM(trainingData, outlierTypes)

warning('off','stats:cvpartition:KFoldMissingGrp')

% Extract predictors and response
% This code processes the data into the right shape for training the
% model.
inputTable = trainingData;
predictorNames = inputTable.Properties.VariableNames(2:end);
predictors = inputTable(:, predictorNames);
response = inputTable.Type;

% Outliers
%response(ismember(response,outlierTypes)) = min(outlierTypes);
%classnames = unique(response);
classnames = min(response):max(response);
nClasses = length(classnames);

% Create cost matrix for outliers
costM = ones(nClasses,nClasses);
outlierTypes = outlierTypes(ismember(outlierTypes,classnames));
outIdx = nchoosek(outlierTypes,2);
costM(outIdx(:,1),outIdx(:,2)) = 1E-5;
costM(outIdx(:,2),outIdx(:,1)) = 1E-5;
costM = costM - diag(diag(costM));
if size(costM,1)>6
    costM(7,2) = 2;
end

% Train a classifier
% This code specifies all the classifier options and trains the classifier.
template = templateSVM(...
    'KernelFunction', 'linear', ...
    'PolynomialOrder', [], ...
    'KernelScale', 'auto', ...
    'BoxConstraint', 1, ...
    'Standardize', true);
classificationSVM = fitcecoc(...
    predictors, ...
    response, ...
    'Learners', template, ...
    'Coding', 'onevsone', ...
    'Cost', costM, ...
    'ClassNames', classnames);

% Create the result struct with predict function
predictorExtractionFcn = @(t) t(:, predictorNames);
svmPredictFcn = @(x) predict(classificationSVM, x);
trainedClassifier.predictFcn = @(x) svmPredictFcn(predictorExtractionFcn(x));

% Add additional fields to the result struct
trainedClassifier.RequiredVariables = predictorNames;
trainedClassifier.ClassificationSVM = classificationSVM;

% Perform cross-validation
partitionedModel = crossval(trainedClassifier.ClassificationSVM, 'KFold', 5);

% Compute validation predictions
[validationPredictions, ~] = kfoldPredict(partitionedModel);

% Display validation accuracy
validationPredictions(ismember(validationPredictions,outlierTypes)) = min(outlierTypes);
response(ismember(response,outlierTypes)) = min(outlierTypes);
fprintf("Correctly classified %.2f%% of cells after cross validation \n",...
    sum(validationPredictions==response)*100/length(response))

disp(sum(validationPredictions==1))
disp(sum(validationPredictions==2))
disp(sum(validationPredictions==3))
disp(sum(validationPredictions==4))

% Compute validation accuracy
validationAccuracy = 1 - kfoldLoss(partitionedModel, 'LossFun', 'ClassifError');

end

function [trainedClassifier, validationAccuracy] = trainClassifierTrees(trainingData, outlierTypes)

% Extract predictors and response
% This code processes the data into the right shape for training the
% model.
inputTable = trainingData;
predictorNames = {'Layer', 'Structure', 'Max_ToPro', 'Mean_ToPro', 'Std_ToPro', 'Middle_ToPro', 'MoC_ToPro', 'Solidity_ToPro', 'FilledArea_ToPro', 'CenDist_ToPro', 'Max_Ctip2', 'Mean_Ctip2', 'Std_Ctip2', 'Middle_Ctip2', 'MoC_Ctip2', 'Solidity_Ctip2', 'FilledArea_Ctip2', 'CenDist_Ctip2', 'Max_Cux1', 'Mean_Cux1', 'Std_Cux1', 'Middle_Cux1', 'MoC_Cux1', 'Solidity_Cux1', 'FilledArea_Cux1', 'CenDist_Cux1'};
predictors = inputTable(:, predictorNames);
response = inputTable.Type;
isCategoricalPredictor = [false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false];

% Train a classifier
% This code specifies all the classifier options and trains the classifier.
template = templateTree(...
    'MaxNumSplits', 20);
classificationEnsemble = fitcensemble(...
    predictors, ...
    response, ...
    'Method', 'AdaBoostM2', ...
    'NumLearningCycles', 60, ...
    'Learners', template, ...
    'LearnRate', 0.1, ...
    'ClassNames', [1; 2; 3; 4]);

% Create the result struct with predict function
predictorExtractionFcn = @(t) t(:, predictorNames);
ensemblePredictFcn = @(x) predict(classificationEnsemble, x);
trainedClassifier.predictFcn = @(x) ensemblePredictFcn(predictorExtractionFcn(x));

% Add additional fields to the result struct
trainedClassifier.RequiredVariables = {'CenDist_Ctip2', 'CenDist_Cux1', 'CenDist_ToPro', 'FilledArea_Ctip2', 'FilledArea_Cux1', 'FilledArea_ToPro', 'Layer', 'Max_Ctip2', 'Max_Cux1', 'Max_ToPro', 'Mean_Ctip2', 'Mean_Cux1', 'Mean_ToPro', 'Middle_Ctip2', 'Middle_Cux1', 'Middle_ToPro', 'MoC_Ctip2', 'MoC_Cux1', 'MoC_ToPro', 'Solidity_Ctip2', 'Solidity_Cux1', 'Solidity_ToPro', 'Std_Ctip2', 'Std_Cux1', 'Std_ToPro', 'Structure'};
trainedClassifier.ClassificationEnsemble = classificationEnsemble;
trainedClassifier.About = 'This struct is a trained model exported from Classification Learner R2020a.';
trainedClassifier.HowToPredict = sprintf('To make predictions on a new table, T, use: \n  yfit = c.predictFcn(T) \nreplacing ''c'' with the name of the variable that is this struct, e.g. ''trainedModel''. \n \nThe table, T, must contain the variables returned by: \n  c.RequiredVariables \nVariable formats (e.g. matrix/vector, datatype) must match the original training data. \nAdditional variables are ignored. \n \nFor more information, see <a href="matlab:helpview(fullfile(docroot, ''stats'', ''stats.map''), ''appclassification_exportmodeltoworkspace'')">How to predict using an exported model</a>.');

% Extract predictors and response
% This code processes the data into the right shape for training the
% model.
inputTable = trainingData;
predictorNames = {'Layer', 'Structure', 'Max_ToPro', 'Mean_ToPro', 'Std_ToPro', 'Middle_ToPro', 'MoC_ToPro', 'Solidity_ToPro', 'FilledArea_ToPro', 'CenDist_ToPro', 'Max_Ctip2', 'Mean_Ctip2', 'Std_Ctip2', 'Middle_Ctip2', 'MoC_Ctip2', 'Solidity_Ctip2', 'FilledArea_Ctip2', 'CenDist_Ctip2', 'Max_Cux1', 'Mean_Cux1', 'Std_Cux1', 'Middle_Cux1', 'MoC_Cux1', 'Solidity_Cux1', 'FilledArea_Cux1', 'CenDist_Cux1'};
predictors = inputTable(:, predictorNames);
response = inputTable.Type;
isCategoricalPredictor = [false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false];

% Perform cross-validation
partitionedModel = crossval(trainedClassifier.ClassificationEnsemble, 'KFold', 5);

% Compute validation predictions
[validationPredictions, validationScores] = kfoldPredict(partitionedModel);

% Compute validation accuracy
validationAccuracy = 1 - kfoldLoss(partitionedModel, 'LossFun', 'ClassifError');
end

