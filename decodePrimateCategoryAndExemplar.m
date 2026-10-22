function primateInfo = decodePrimateCategoryAndExemplar( ...
    rasters,allimages,goodch,rollingwindowstart,binsz,psthstart,psthend,varargin)
% decodePrimateCategoryAndExemplar
%
% Trial-level analysis of:
%   1) 5-way primate category decoding
%   2) exemplar decoding across all retained primate images
%   3) within-category exemplar decoding
%
% Uses independent trials. By default, exemplars must have >=5 presentations.
% For each repeat, 4 trials/exemplar form the training template and 1 independent
% trial/exemplar is used for testing. The split is redrawn NRepeats times.
%
% allimages must correspond one-for-one to size(rasters,3).
%
% The rolling windows are EXACTLY:
%   rollingwindowstart(w) : rollingwindowstart(w)+binsz
%
% Example:
% primateInfo = decodePrimateCategoryAndExemplar( ...
%     rasters,allimages,goodch,rollingwindowstart,binsz,psthstart,psthend, ...
%     'PopulationName','red_220227');

p = inputParser;

root = '/n/data2/hms/neurobio/livingstone/marge/margemonkeys';

addParameter(p,'HFaceFile',fullfile(root,'hfaceims.mat'));
addParameter(p,'MFaceFile',fullfile(root,'mfaceims.mat'));
addParameter(p,'PFaceFile',fullfile(root,'pfaceims.mat'));
addParameter(p,'HumanFile',fullfile(root,'humanims.mat'));
addParameter(p,'MonkeyFile',fullfile(root,'monkeyims.mat'));

addParameter(p,'MinTrialsPerExemplar',5,@isscalar);
addParameter(p,'NTrainTrials',4,@isscalar);
addParameter(p,'NRepeats',200,@isscalar);
addParameter(p,'NMatchedSets',100,@isscalar);
addParameter(p,'NMatchedExemplars',5,@isscalar);
addParameter(p,'Seed',260825,@isscalar);
addParameter(p,'PopulationName','population',@(x)ischar(x)||isstring(x));
addParameter(p,'MakeFigure',true,@islogical);
addParameter(p,'FigureFolder', ...
    '/n/data2/hms/neurobio/livingstone/marge/figimages/temp', ...
    @(x)ischar(x)||isstring(x));

parse(p,varargin{:});
S = p.Results;

if size(rasters,3) ~= numel(allimages)
    error('size(rasters,3) must equal numel(allimages).');
end

%% Load curated primate sets
A = load(S.HFaceFile);
B = load(S.MFaceFile);
C = load(S.PFaceFile);
D = load(S.HumanFile);
E = load(S.MonkeyFile);

hface  = getNamedVariable(A,["hface","hfaceims"]);
mface  = getNamedVariable(B,["mface","mfaceims"]);
pface  = getNamedVariable(C,["Pface","pface","pfaceims"]);
human  = getNamedVariable(D,["human","humanims"]);
monkey = getNamedVariable(E,["monkey","monkeyims"]);

sets = {hface,mface,pface,human,monkey};
categoryNames = ["hface","mface","pface","human","monkey"];

trialNames = cleanNamesPrimate(allimages);

% Assign category to each presentation.
trialCategory = zeros(numel(trialNames),1,'uint8');

for cc = 1:5
    nm = cleanNamesPrimate(sets{cc});
    hit = ismember(trialNames,nm);

   overlap = hit & trialCategory>0;

% Mark ambiguous category membership with 255
trialCategory(overlap) = 255;

% Assign only truly unassigned images
newhit = hit & trialCategory==0;
trialCategory(newhit) = cc;
end

primateTrialIdx = find(trialCategory>=1 & trialCategory<=5);

fprintf('\n============================================================\n');
fprintf('PRIMATE CATEGORY + EXEMPLAR DECODING: %s\n',string(S.PopulationName));
fprintf('============================================================\n');
fprintf('All presentations       = %d\n',numel(allimages));
fprintf('Primate presentations   = %d\n',numel(primateTrialIdx));

%% Build exemplar table
pNames = trialNames(primateTrialIdx);
pCat   = double(trialCategory(primateTrialIdx));

[exemplarNames,~,exemplarLocalID] = unique(pNames,'stable');
nEx = numel(exemplarNames);

exemplarCounts = accumarray(exemplarLocalID,1);

exemplarCategory = zeros(nEx,1);
for ee = 1:nEx
    cats = unique(pCat(exemplarLocalID==ee));
    if numel(cats)~=1
        error('Exemplar %s maps to multiple categories.',exemplarNames(ee));
    end
    exemplarCategory(ee) = cats;
end

keepEx = exemplarCounts >= S.MinTrialsPerExemplar;
keptOldIDs = find(keepEx);

keptNames = exemplarNames(keepEx);
keptCounts = exemplarCounts(keepEx);
keptCategory = exemplarCategory(keepEx);
nKept = numel(keptNames);

fprintf('Unique primate exemplars = %d\n',nEx);
fprintf('Exemplars >=%d trials    = %d\n',S.MinTrialsPerExemplar,nKept);

for cc = 1:5
    fprintf('  %-7s retained = %d\n',categoryNames(cc),sum(keptCategory==cc));
end
myflag.montage=0;
if myflag.montage
%% ------------------------------------------------------------------------
% OPTIONAL CHECK: MONTAGES OF THE ACTUAL RETAINED EXEMPLARS
% -------------------------------------------------------------------------

imageRoot = ...
    '/n/data2/hms/neurobio/livingstone/Stimuli/manyOO';

imds = imageDatastore(imageRoot, ...
    'IncludeSubfolders',true, ...
    'FileExtensions',{'.jpg','.jpeg','.png','.bmp','.tif','.tiff'});

% Basenames of every image in manyOO
allFileNames = strings(numel(imds.Files),1);

for ii = 1:numel(imds.Files)
    [~,nm,ext] = fileparts(imds.Files{ii});
    allFileNames(ii) = string(nm) + string(ext);
end

% Make one or more montage pages for each of the five categories
nPerPage = 30;

for cc = 1:5

    theseNames = keptNames(keptCategory==cc);

    fprintf('\nMaking montage for %s: %d retained exemplars\n', ...
        categoryNames(cc),numel(theseNames));

    % Find each retained exemplar in manyOO
    theseFiles = strings(0,1);
    foundNames = strings(0,1);

    for ii = 1:numel(theseNames)

        hit = find(strcmp(allFileNames,theseNames(ii)));

        if isempty(hit)

            fprintf('NOT FOUND: %s\n',theseNames(ii));

        elseif numel(hit) > 1

            fprintf('WARNING: multiple files named %s; using first\n', ...
                theseNames(ii));

            theseFiles(end+1,1) = string(imds.Files{hit(1)});
            foundNames(end+1,1) = theseNames(ii);

        else

            theseFiles(end+1,1) = string(imds.Files{hit});
            foundNames(end+1,1) = theseNames(ii);

        end
    end

    nPages = ceil(numel(theseFiles)/nPerPage);

    for pg = 1:nPages

        first = (pg-1)*nPerPage + 1;
        last  = min(pg*nPerPage,numel(theseFiles));

        figure('Color','w');

        tiledlayout(5,6, ...
            'TileSpacing','compact', ...
            'Padding','compact');

        for jj = first:last

            nexttile

            im = imread(theseFiles(jj));
            imshow(im)

            title(foundNames(jj), ...
                'Interpreter','none', ...
                'FontSize',7);

        end

        sgtitle(sprintf('%s: %s retained exemplars, page %d/%d', ...
            string(S.PopulationName), ...
            categoryNames(cc), ...
            pg,nPages), ...
            'Interpreter','none');

    end
end
end


if nKept < 2
    error('Too few retained exemplars.');
end

if any(keptCounts < S.NTrainTrials+1)
    error('Some retained exemplars lack NTrainTrials + 1 test trial.');
end

oldToNew = zeros(nEx,1);
oldToNew(keptOldIDs) = 1:nKept;

retainLocalTrial = keepEx(exemplarLocalID);
retainedPrimateLocalIdx = find(retainLocalTrial);

retainedOriginalTrialIdx = primateTrialIdx(retainedPrimateLocalIdx);
retainedExemplarID = oldToNew(exemplarLocalID(retainedPrimateLocalIdx));

trialsByExemplar = cell(nKept,1);
for ee = 1:nKept
    trialsByExemplar{ee} = ...
        retainedOriginalTrialIdx(retainedExemplarID==ee);
end

%% Preselect independent train/test trials for every repeat
rng(S.Seed,'twister');

trainTrials = cell(S.NRepeats,nKept);
testTrial   = zeros(S.NRepeats,nKept);

for rr = 1:S.NRepeats
    for ee = 1:nKept
        ii = trialsByExemplar{ee};
        ord = randperm(numel(ii),S.NTrainTrials+1);
        trainTrials{rr,ee} = ii(ord(1:S.NTrainTrials));
        testTrial(rr,ee) = ii(ord(S.NTrainTrials+1));
    end
end

%% Storage
nWindows = numel(rollingwindowstart);

categoryAccuracyRepeat = nan(S.NRepeats,nWindows);
exemplarAccuracyRepeat = nan(S.NRepeats,nWindows);
withinCategoryAccuracyRepeat = nan(S.NRepeats,nWindows);
matched5AccuracyRepeat = nan(S.NRepeats,nWindows);
exemplarMRRRepeat = nan(S.NRepeats,nWindows);
withinCategoryMRRRepeat = nan(S.NRepeats,nWindows);

fprintf('Rolling windows          = %d\n',nWindows);
fprintf('Train/test per exemplar  = %d + 1 independent trial\n',S.NTrainTrials);
fprintf('Random splits            = %d\n\n',S.NRepeats);

%% Main loop
for ww = 1:nWindows

    window = rollingwindowstart(ww):rollingwindowstart(ww)+binsz;

    Xwin = squeeze(mean(rasters(goodch,window,:),2,'omitnan'));

    if isvector(Xwin)
        Xwin = reshape(Xwin,numel(goodch),[]);
    end

    for rr = 1:S.NRepeats

        templates = nan(nKept,numel(goodch));
        tests = nan(nKept,numel(goodch));

        for ee = 1:nKept
            templates(ee,:) = mean( ...
                Xwin(:,trainTrials{rr,ee}),2,'omitnan')';
            tests(ee,:) = Xwin(:,testTrial(rr,ee))';
        end

        % Standardize neurons from training templates only.
        mu = mean(templates,1,'omitnan');
        sd = std(templates,0,1,'omitnan');
        keepNeuron = isfinite(mu) & isfinite(sd) & sd>0;

        T = templates(:,keepNeuron);
        Q = tests(:,keepNeuron);

        mu = mu(keepNeuron);
        sd = sd(keepNeuron);

        for jj = 1:size(T,2)
            T(~isfinite(T(:,jj)),jj) = mu(jj);
            Q(~isfinite(Q(:,jj)),jj) = mu(jj);
        end

        T = (T-mu)./sd;
        Q = (Q-mu)./sd;

        % 1) 5-way category decoding
        catCentroid = nan(5,size(T,2));
        for cc = 1:5
            catCentroid(cc,:) = mean(T(keptCategory==cc,:),1);
        end

        dCat = squaredDistanceFast(Q,catCentroid);
        [~,predCat] = min(dCat,[],2);

        categoryAccuracyRepeat(rr,ww) = mean(predCat==keptCategory);

        % 2) All-primate exemplar identification
        dAll = squaredDistanceFast(Q,T);

        [~,predEx] = min(dAll,[],2);
        exemplarAccuracyRepeat(rr,ww) = ...
            mean(predEx(:)==(1:nKept)');

        exemplarMRRRepeat(rr,ww) = ...
            meanReciprocalRankFromDistance(dAll,(1:nKept)');

        % 3) Within-category exemplar identification
        nCorrect = 0;
        rrValues = nan(nKept,1);

        for ee = 1:nKept
            cc = keptCategory(ee);
            candidates = find(keptCategory==cc);

            [~,ord] = sort(dAll(ee,candidates),'ascend');
            rankedCandidates = candidates(ord);

            rankCorrect = find(rankedCandidates==ee,1);

            if rankCorrect==1
                nCorrect = nCorrect+1;
            end

            rrValues(ee) = 1/rankCorrect;
        end

        withinCategoryAccuracyRepeat(rr,ww) = nCorrect/nKept;
        withinCategoryMRRRepeat(rr,ww) = mean(rrValues,'omitnan');

        %% 4) MATCHED 5-WAY WITHIN-CATEGORY EXEMPLAR DECODING
% Exactly 5 alternatives, just like category decoding.
% Therefore chance = 1/5 = 0.20.
%
% For each category, randomly draw 5 exemplars NMatchedSets times.
% Classify each held-out exemplar trial against the 5 corresponding
% exemplar templates. Average across subsets and categories.

matchedAcc = nan(5,S.NMatchedSets);

for cc = 1:5

    candidates = find(keptCategory==cc);

    if numel(candidates) < S.NMatchedExemplars
        continue
    end

    for ss = 1:S.NMatchedSets

        pick = candidates(randperm(numel(candidates), ...
                                  S.NMatchedExemplars));

        % Distances from the 5 held-out trials to the
        % 5 corresponding exemplar templates
        d5 = squaredDistanceFast(Q(pick,:),T(pick,:));

        [~,pred5] = min(d5,[],2);

        matchedAcc(cc,ss) = ...
            mean(pred5(:)==(1:S.NMatchedExemplars)');

    end
end

matched5AccuracyRepeat(rr,ww) = ...
    mean(matchedAcc(:),'omitnan');
    end

    if ww==1 || mod(ww,10)==0 || ww==nWindows
        fprintf('  window %d/%d complete\n',ww,nWindows);
    end
end

%% Average repeats
categoryAccuracy = mean(categoryAccuracyRepeat,1,'omitnan')';
exemplarAccuracy = mean(exemplarAccuracyRepeat,1,'omitnan')';
withinCategoryAccuracy = mean(withinCategoryAccuracyRepeat,1,'omitnan')';

matched5Accuracy = ...
    mean(matched5AccuracyRepeat,1,'omitnan')';

matched5SEM = ...
    semAcrossRepeats(matched5AccuracyRepeat);

exemplarMRR = mean(exemplarMRRRepeat,1,'omitnan')';
withinCategoryMRR = mean(withinCategoryMRRRepeat,1,'omitnan')';

categorySEM = semAcrossRepeats(categoryAccuracyRepeat);
exemplarSEM = semAcrossRepeats(exemplarAccuracyRepeat);
withinCategorySEM = semAcrossRepeats(withinCategoryAccuracyRepeat);

exemplarMRRSEM = semAcrossRepeats(exemplarMRRRepeat);
withinCategoryMRRSEM = semAcrossRepeats(withinCategoryMRRRepeat);

%% Chance levels
categoryChance = 1/5;
exemplarChance = 1/nKept;
matched5Chance = 1/S.NMatchedExemplars;
matched5OverChance = matched5Accuracy/matched5Chance;
nByCat = accumarray(keptCategory,1,[5 1]);
withinCategoryChance = mean(1 ./ nByCat(keptCategory));

mrrChanceByCat = zeros(5,1);
for cc = 1:5
    N = nByCat(cc);
    mrrChanceByCat(cc) = sum(1./(1:N))/N;
end
withinCategoryMRRChance = mean(mrrChanceByCat(keptCategory));

N = nKept;
exemplarMRRChance = sum(1./(1:N))/N;

%% Exact PSTH crop
crop = psthstart:psthend;

categoryPerPSTH = categoryAccuracy(crop);
exemplarPerPSTH = exemplarAccuracy(crop);
withinCategoryPerPSTH = withinCategoryAccuracy(crop);
matched5PerPSTH = matched5Accuracy(crop);
exemplarMRRPerPSTH = exemplarMRR(crop);
withinCategoryMRRPerPSTH = withinCategoryMRR(crop);

%% Chance-normalized versions
categoryOverChance = categoryAccuracy/categoryChance;
exemplarOverChance = exemplarAccuracy/exemplarChance;
withinCategoryOverChance = withinCategoryAccuracy/withinCategoryChance;

exemplarMRROverChance = exemplarMRR/exemplarMRRChance;
withinCategoryMRROverChance = withinCategoryMRR/withinCategoryMRRChance;

%% Slopes over exact response interval
tNorm = linspace(0,1,numel(crop))';

[catSlope,catSlopeCI,catSlopeP] = slopeStats(tNorm,categoryPerPSTH);
[exSlope,exSlopeCI,exSlopeP] = slopeStats(tNorm,exemplarPerPSTH);
[withinSlope,withinSlopeCI,withinSlopeP] = ...
    slopeStats(tNorm,withinCategoryPerPSTH);

[exMRRSlope,exMRRSlopeCI,exMRRSlopeP] = ...
    slopeStats(tNorm,exemplarMRRPerPSTH);
[withinMRRSlope,withinMRRSlopeCI,withinMRRSlopeP] = ...
    slopeStats(tNorm,withinCategoryMRRPerPSTH);

%% Package output
primateInfo = struct;

primateInfo.categoryNames = categoryNames;
primateInfo.exemplarNames = keptNames;
primateInfo.exemplarCategory = keptCategory;
primateInfo.exemplarTrialCounts = keptCounts;
primateInfo.nExemplars = nKept;
primateInfo.nExemplarsByCategory = nByCat;

primateInfo.matched5Accuracy = matched5Accuracy;
primateInfo.matched5SEM = matched5SEM;
primateInfo.matched5AccuracyRepeat = matched5AccuracyRepeat;
primateInfo.matched5Chance = matched5Chance;
primateInfo.matched5OverChance = matched5OverChance;
primateInfo.nMatchedSets = S.NMatchedSets;
primateInfo.nMatchedExemplars = S.NMatchedExemplars;

primateInfo.minTrialsPerExemplar = S.MinTrialsPerExemplar;
primateInfo.nTrainTrials = S.NTrainTrials;
primateInfo.nRepeats = S.NRepeats;
primateInfo.seed = S.Seed;

primateInfo.rollingwindowstart = rollingwindowstart;
primateInfo.binsz = binsz;
primateInfo.psthstart = psthstart;
primateInfo.psthend = psthend;

primateInfo.categoryAccuracy = categoryAccuracy;
primateInfo.exemplarAccuracy = exemplarAccuracy;
primateInfo.withinCategoryAccuracy = withinCategoryAccuracy;

primateInfo.exemplarMRR = exemplarMRR;
primateInfo.withinCategoryMRR = withinCategoryMRR;

primateInfo.categorySEM = categorySEM;
primateInfo.exemplarSEM = exemplarSEM;
primateInfo.withinCategorySEM = withinCategorySEM;
primateInfo.exemplarMRRSEM = exemplarMRRSEM;
primateInfo.withinCategoryMRRSEM = withinCategoryMRRSEM;

primateInfo.categoryAccuracyRepeat = categoryAccuracyRepeat;
primateInfo.exemplarAccuracyRepeat = exemplarAccuracyRepeat;
primateInfo.withinCategoryAccuracyRepeat = withinCategoryAccuracyRepeat;
primateInfo.exemplarMRRRepeat = exemplarMRRRepeat;
primateInfo.withinCategoryMRRRepeat = withinCategoryMRRRepeat;

primateInfo.categoryChance = categoryChance;
primateInfo.exemplarChance = exemplarChance;
primateInfo.withinCategoryChance = withinCategoryChance;
primateInfo.exemplarMRRChance = exemplarMRRChance;
primateInfo.withinCategoryMRRChance = withinCategoryMRRChance;

primateInfo.categoryOverChance = categoryOverChance;
primateInfo.exemplarOverChance = exemplarOverChance;
primateInfo.withinCategoryOverChance = withinCategoryOverChance;
primateInfo.exemplarMRROverChance = exemplarMRROverChance;
primateInfo.withinCategoryMRROverChance = withinCategoryMRROverChance;

primateInfo.categoryPerPSTH = categoryPerPSTH;
primateInfo.exemplarPerPSTH = exemplarPerPSTH;
primateInfo.withinCategoryPerPSTH = withinCategoryPerPSTH;
primateInfo.exemplarMRRPerPSTH = exemplarMRRPerPSTH;
primateInfo.withinCategoryMRRPerPSTH = withinCategoryMRRPerPSTH;
primateInfo.matched5PerPSTH = matched5PerPSTH;
primateInfo.categorySlope = catSlope;
primateInfo.categorySlopeCI = catSlopeCI;
primateInfo.categorySlopeP = catSlopeP;

primateInfo.exemplarSlope = exSlope;
primateInfo.exemplarSlopeCI = exSlopeCI;
primateInfo.exemplarSlopeP = exSlopeP;

primateInfo.withinCategorySlope = withinSlope;
primateInfo.withinCategorySlopeCI = withinSlopeCI;
primateInfo.withinCategorySlopeP = withinSlopeP;

primateInfo.exemplarMRRSlope = exMRRSlope;
primateInfo.exemplarMRRSlopeCI = exMRRSlopeCI;
primateInfo.exemplarMRRSlopeP = exMRRSlopeP;

primateInfo.withinCategoryMRRSlope = withinMRRSlope;
primateInfo.withinCategoryMRRSlopeCI = withinMRRSlopeCI;
primateInfo.withinCategoryMRRSlopeP = withinMRRSlopeP;

%% Quick comparison figure
if S.MakeFigure
    if ~exist(S.FigureFolder,'dir')
        mkdir(S.FigureFolder);
    end

    figure; hold on

    plot(rollingwindowstart,categoryOverChance, ...
        'b-','LineWidth',2.5);

    plot(rollingwindowstart,withinCategoryMRROverChance, ...
        'r-','LineWidth',2.5);

    plot([rollingwindowstart(1) rollingwindowstart(end)], ...
        [1 1],'k--','LineWidth',1.5);

    xlabel('time after stimulus onset, ms');
    ylabel('information / chance');

    legend('5-way primate category', ...
           'within-category exemplar (MRR)', ...
           'chance','Location','best');

    title(sprintf('%s: category versus exemplar information', ...
        string(S.PopulationName)), ...
        'Interpreter','none');

    set(gca,'TickDir','out','LineWidth',2,'FontSize',13);
    box on

    filename = fullfile(S.FigureFolder, ...
        [char(S.PopulationName) '_primateCategory_vs_exemplar.png']);

    print(gcf,filename,'-dpng','-r300');
    fprintf('\nSaved comparison figure:\n%s\n',filename);
end

figure; hold on

plot(rollingwindowstart,categoryAccuracy, ...
    'b-','LineWidth',2.5);

plot(rollingwindowstart,matched5Accuracy, ...
    'r-','LineWidth',2.5);

yline(0.2,'k--','LineWidth',1.5);

xlabel('time after stimulus onset, ms');
ylabel('decoding accuracy');

legend('5-way category', ...
       '5-way within-category exemplar', ...
       'chance', ...
       'Location','best');

title(sprintf('%s: matched 5-way decoding', ...
    string(S.PopulationName)), ...
    'Interpreter','none');

set(gca,'TickDir','out','LineWidth',2,'FontSize',13);
box on


fprintf('\nChance levels:\n');
fprintf('  5-way category top-1           = %.5f\n',categoryChance);
fprintf('  all-exemplar top-1             = %.5f\n',exemplarChance);
fprintf('  within-category exemplar top-1 = %.5f\n',withinCategoryChance);
fprintf('  all-exemplar MRR               = %.5f\n',exemplarMRRChance);
fprintf('  within-category exemplar MRR   = %.5f\n',withinCategoryMRRChance);

fprintf('\nDONE\n');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOCAL FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function x = getNamedVariable(S,likelyNames)
f = string(fieldnames(S));
for k = 1:numel(likelyNames)
    hit = find(strcmp(f,likelyNames(k)),1);
    if ~isempty(hit)
        x = S.(char(f(hit)));
        return
    end
end
if numel(f)==1
    x = S.(char(f));
    return
end
error('Could not identify image-list variable in MAT file.');
end

function s = cleanNamesPrimate(x)
s = string(x(:));
for ii = 1:numel(s)
    thisName = s(ii);

    slash = strfind(thisName,'/');
    if ~isempty(slash)
        thisName = extractAfter(thisName,slash(end));
    end

    slash = strfind(thisName,'\');
    if ~isempty(slash)
        thisName = extractAfter(thisName,slash(end));
    end

    thisName = replace(thisName,'..','.');
    s(ii) = thisName;
end
end

function d = squaredDistanceFast(A,B)
a2 = sum(A.^2,2);
b2 = sum(B.^2,2)';
d = a2 + b2 - 2*(A*B');
d(d<0 & d>-1e-10) = 0;
end

function mrr = meanReciprocalRankFromDistance(d,correctIndex)
n = size(d,1);
rr = nan(n,1);
for ii = 1:n
    [~,ord] = sort(d(ii,:),'ascend');
    rankCorrect = find(ord==correctIndex(ii),1);
    rr(ii) = 1/rankCorrect;
end
mrr = mean(rr,'omitnan');
end

function sem = semAcrossRepeats(X)
sem = std(X,0,1,'omitnan')' ./ sqrt(sum(isfinite(X),1))';
end

function [slope,ciSlope,pSlope] = slopeStats(t,y)
good = isfinite(t) & isfinite(y);
if sum(good)<3
    slope = NaN;
    ciSlope = [NaN NaN];
    pSlope = NaN;
    return
end
mdl = fitlm(t(good),y(good));
slope = mdl.Coefficients.Estimate(2);
ci = coefCI(mdl,0.05);
ciSlope = ci(2,:);
pSlope = mdl.Coefficients.pValue(2);
end
