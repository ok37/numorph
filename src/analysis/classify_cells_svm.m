function ct = classify_cells_svm(centroids, path_table, config)
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

% Read classifications + patch feautre info
itable = readtable(fullfile(config.output_directory,'classifier',...
    sprintf('%s_patch_info.csv',config.sample_name)));
ftable = readtable(fullfile(config.output_directory,'classifier',...
    sprintf('%s_patch_features.csv',config.sample_name)));
ctable = readtable(fullfile(config.output_directory,'classifier',...
    sprintf('%s_classifications.csv',config.sample_name)));

% Check table sizes
assert(size(ctable,2) == 1,"Classification table contains more than 1 column")
n_classified = height(ctable);
fprintf("Using %d cells for training SVM model \n",n_classified)
assert(height(ftable) == n_classified, "Number of rows in the feature table do not match")
assert(height(itable) == n_classified, "Number of rows in the info table do not match")

% Combine manual classifications with patch features
ctable.Properties.VariableNames = {'Type'};
ftable.Layer = num2str(ftable.Layer);
ftable.Structure = num2str(ftable.Structure);
trainingData = horzcat(ctable,ftable);
[ncols, nrows] = size(ftable);

% Train the model
trainedClassifier = trainClassifier(trainingData, outlierTypes);

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
%structures = num2str(bin_annotation_structures(centroids(:,4),'cortex_large'));
%layers = num2str(bin_annotation_structures(centroids(:,4),'layers'));

% Get patch features for remaining centroids
z_pos = unique(centroids(:,3));
stable = {(zeros([sum(ct==0),8]))};
stable = repmat(stable,1,length(config.markers));
%stable = cell(1,length(config.markers));
a = 1;
for i = 1:length(z_pos)
    fprintf("Working on z position %d \n",z_pos(i))
    idx = centroids(:,3) == z_pos(i);
    if ~any(ct(idx))
        continue
    end
    cen_sub = centroids(idx,:);
    % Load images and get features
    for j = 1:length(config.markers)
        file = path_table(path_table.z == z_pos(i) &...
            path_table.markers == config.markers(j),:);
        file = file.file{1};
        f = zeros([size(cen_sub,1),8]);
        parfor k = 1:size(cen_sub,1)
            pos = cen_sub(k,[1,2]);
            ranges = {[pos(1)-s,pos(1)+s], [pos(2)-s,pos(2)+s]};
            img = imread(file,'PixelRegion',ranges);
            f(k,:) = measure_patch_features(img, s, false);
        end
        b = a+size(cen_sub,1)-1;
        stable{j}(a:b,:) = f;
        a = a+1;
        %if isempty(stable{j})
        %    stable{j} = f;
        %else
        %    stable{j} = vertcat(stable{j},f);
        %end
    end
end

% Run prediction
feature_table = array2table([layers,structures],'VariableNames',{'Layer', 'Structure'});
for i = 1:length(stable)
   feature_table = horzcat(feature_table,stable{i});
end
predictions = trainedClassifier.predictFcn(feature_table);
ct(~k_idx) = predictions;

end


function [trainedClassifier, validationAccuracy] = trainClassifier(trainingData, outlierTypes)

warning('off','stats:cvpartition:KFoldMissingGrp')

% Extract predictors and response
% This code processes the data into the right shape for training the
% model.
inputTable = trainingData;
predictorNames = inputTable.Properties.VariableNames(2:end);
predictors = inputTable(:, predictorNames);
response = inputTable.Type;
isCategoricalPredictor = logical(varfun(@ischar,inputTable,'output','uniform'));
classnames = unique(inputTable.Type);

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
    'Coding', 'onevsall', ...
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

% Compute validation accuracy
validationAccuracy = 1 - kfoldLoss(partitionedModel, 'LossFun', 'ClassifError');

end



