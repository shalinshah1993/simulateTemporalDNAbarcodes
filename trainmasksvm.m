%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Shalin Shah
% Affiliation: Dept. of Electrical & Computer Engineering, Duke University
% Email: shalin.shah@duke.edu
% Last modified: 01/08/2019
% Matlab version used: R2017a
%
% Description: This code takes data as input, divides in training/ test set
% and trains an SVM model. The output of test data set is analyzed using
% a confusion matrix plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;

nStaple = 5;
nSample = 600;
label = labeldata(nStaple, nSample);

% UN-COMMENT THIS IF A DEVICE TYPE SHOULD BE DISCARDED. FOR EXAMPLE, 7 nt.
% delete devices with combination of 7
% maskSeven = (any(ismember(label, '7'), 2));
% label(maskSeven, :) = [];

% load the csv dataclear;
addpath('libs/confusionmatStats');

% read input data from CSV file
data = csvread('traindata600.csv');
data = data(:, 1:nStaple);
data(maskSeven, :) = [];

% randomly shuffle input data rows
N = size(data, 1);
shuffleIdx = randperm(N);
data = data(shuffleIdx, :);
label = label(shuffleIdx, :);

% split input data into train/test set (70/30)
splitIdx = ceil(0.7 * N);

% split data into features (X), label (Y)
Xtrain = log(data(1:splitIdx, :));
Xtest = log(data(splitIdx+1:end, :));
Ytrain = label(1:splitIdx, :);
Ytest = label(splitIdx+1:end, :);

% set gaussian SVM kernal, enable parallel computing
options = statset('UseParallel',true);
svmTmp = templateSVM('KernelFunction','gaussian', 'KernelScale', 'auto' );

% fit model to ECOC multi-class SVM
svmMdl = fitcecoc(Xtrain, Ytrain, 'Coding', 'onevsall', 'Learners', svmTmp, ...
                                                            'Options', options);
CVMdl = crossval(svmMdl, 'Options', options);
loss = kfoldLoss(CVMdl);

% predict test data with trained model
[Ypred, score] = predict(svmMdl, Xtest);
stats = confusionmatStats(Ytest, Ypred);
imagesc(stats.confusionMat); colorbar;


function label = labeldata(maxStapleCount, nSample)
    % This method generates string labels for each device
    %
    % maxStapleCount is the number of dimensions for each data point
    %
    % nSample is the number of samples per device
    iDevice = 1;
    label = {};
    for iSeven = 0 : maxStapleCount
        for iEight = 0 : maxStapleCount
            for iNine = 0 : maxStapleCount
                for iTen = 0 : maxStapleCount
                    if iSeven + iEight + iNine + iTen ~= maxStapleCount
                        continue
                    end

                    label{iDevice, 1} = ([repmat(07, [1 iSeven]) ...
                            repmat(08, [1 iEight]) ...
                            repmat(09, [1 iNine]) ...
                            repmat(10, [1 iTen])]);

                    iDevice = iDevice + 1;
                end
            end
        end
    end
    label = repmat(label, nSample, 1);
    label = num2str(cell2mat(label));
end