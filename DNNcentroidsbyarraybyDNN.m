%% ========================================================================
% CURRENT DNN TEMPORAL CENTROIDS FOR ALL 6 ARRAYS
%
% Recreates the old 6-panel centroid-vs-layer figure using the
% CURRENT centroidByNetwork values used in the final DNN/metric analysis.
%
% IMPORTANT:
% centroidByNetwork is used RAW here -- NO 0-to-1 normalization.
% Therefore y-axis is centroid time in ms.
% ========================================================================

close all
clearvars
clc


%% ------------------------------------------------------------------------
% Current bridge-result files
% -------------------------------------------------------------------------

base = ...
'/n/data2/hms/neurobio/livingstone/marge/margemonkeys/chatmetricsEdges/data';

files = { ...
    '/n/data2/hms/neurobio/livingstone/marge/margemonkeys/paul/data/paul230204_DNN_METRIC_NEURALTIME_bridge_results.mat', ...
    '/n/data2/hms/neurobio/livingstone/marge/margemonkeys/paul/data/paul220428_DNN_METRIC_NEURALTIME_bridge_results.mat', ...
    '/n/data2/hms/neurobio/livingstone/marge/margemonkeys/baby5/data/baby5221104_DNN_METRIC_NEURALTIME_bridge_results.mat', ...
    '/n/data2/hms/neurobio/livingstone/marge/margemonkeys/baby1/data/baby1230629_DNN_METRIC_NEURALTIME_bridge_results.mat', ...
    '/n/data2/hms/neurobio/livingstone/marge/margemonkeys/red/data/red240302_DNN_METRIC_NEURALTIME_bridge_results.mat', ...
    '/n/data2/hms/neurobio/livingstone/marge/margemonkeys/baby1/data/baby1220226_DNN_METRIC_NEURALTIME_bridge_results.mat' ...
    };

popNames = { ...
    'P-PIT', ...
    'P-CIT', ...
    'B5-CIT', ...
    'B1-AIT-NPx', ...
    'R-AIT', ...
    'B1-AIT-MWA'};


%% ------------------------------------------------------------------------
% DNN definitions
% -------------------------------------------------------------------------

netNames = {'AlexNet','ResNet-50','ViT-L'};

% Same colors as your old plot
netColors = [ ...
    0.25 0.30 0.90; ...   % AlexNet blue
    0.00 0.90 0.10; ...   % ResNet green
    0.95 0.05 0.05];      % ViT red


nPop = numel(files);
nNet = 3;


%% ------------------------------------------------------------------------
% Figure
% -------------------------------------------------------------------------

fig = figure( ...
    'Color','w', ...
    'Position',[100 80 1500 760]);

tl = tiledlayout(fig,2,3, ...
    'TileSpacing','compact', ...
    'Padding','compact');


%% ------------------------------------------------------------------------
% Load and plot CURRENT RAW temporal centroids
% -------------------------------------------------------------------------

for pop = 1:nPop

    S = load(files{pop},'centroidByNetwork');

    ax = nexttile(tl,pop);
    hold(ax,'on');


    for net = 1:nNet

        % RAW centroid values
        % NO min-max normalization here
        c = double(S.centroidByNetwork{net}(:));

        % Each DNN can have a different number of analyzed layers
        layer = 1:numel(c);

        plot(ax,layer,c, ...
            '-o', ...
            'Color',netColors(net,:), ...
            'MarkerFaceColor',netColors(net,:), ...
            'MarkerEdgeColor',netColors(net,:), ...
            'MarkerSize',8, ...
            'LineWidth',2);

    end


    %% Formatting

    % title(ax,popNames{pop}, ...
    %     'FontSize',18, ...
    %     'FontWeight','bold', ...
    %     'HorizontalAlignment','left');

    box(ax,'on');

    set(ax, ...
        'FontSize',15, ...
        'LineWidth',2, ...
        'TickDir','out');

    if  pop == 4
        ylabel(ax,'Neural-DNN RSA temporal centroid, ms');
    end

    if pop >= 4
        xlabel(ax,'DNN layer');
    end

end


%% ------------------------------------------------------------------------
% Legend
% -------------------------------------------------------------------------
% 
% lgd = legend( ...
%     netNames, ...
%     'Box','off', ...
%     'FontSize',13);
% 
% lgd.Layout.Tile = 1;


%% ------------------------------------------------------------------------
% Save
% -------------------------------------------------------------------------

saveFolder = ...
'/n/data2/hms/neurobio/livingstone/marge/figimages/temp';

exportgraphics( ...
    fig, ...
    fullfile(saveFolder,'CURRENT_DNN_centroids_all6arrays.png'), ...
    'Resolution',600);

exportgraphics( ...
    fig, ...
    fullfile(saveFolder,'CURRENT_DNN_centroids_all6arrays.pdf'), ...
    'ContentType','vector');

fprintf('\nSaved current DNN centroid figure.\n');


%% ============================================================
% STATISTICS: Does neural temporal centroid increase with DNN depth?
%
% For each population x DNN:
%   fit centroid(ms) = intercept + slope * normalized DNN depth
%
% Positive slope = deeper DNN layers correspond to later neural times.
% Depth is normalized 0-1 so slopes are comparable across networks.
% ============================================================

nPop = numel(files);
nNet = 3;

networkNames = {'AlexNet','ResNet-50','ViT-L'};

allSlopes = nan(nPop,nNet);

for pop = 1:nPop

    S = load(files{pop},'centroidByNetwork');

    for net = 1:nNet

        y = S.centroidByNetwork{net};
        y = y(:);

        % Remove NaNs
        good = isfinite(y);
        y = y(good);

        % Normalize DNN depth from shallow=0 to deep=1
        x = linspace(0,1,numel(y))';

        if numel(y) >= 3
            p = polyfit(x,y,1);
            allSlopes(pop,net) = p(1);
        end

    end
end


%% Display individual slopes

fprintf('\n============================================================\n');
fprintf('TEMPORAL CENTROID vs DNN DEPTH\n');
fprintf('Positive slope = deeper DNN layers correspond to later neural times\n');
fprintf('Slope units = ms across normalized DNN depth (0-1)\n');
fprintf('============================================================\n');

for net = 1:nNet

    fprintf('\n%s\n',networkNames{net});

    for pop = 1:nPop
        fprintf('  Population %d: %+7.3f ms\n', ...
            pop,allSlopes(pop,net));
    end

    s = allSlopes(:,net);
    s = s(isfinite(s));

    fprintf('  Mean +/- SEM: %.3f +/- %.3f ms\n', ...
        mean(s),std(s)/sqrt(numel(s)));

    % One-tailed t-test: H1 slope > 0
    [~,p,ci,stats] = ttest(s,0,'Tail','right');

    fprintf('  One-tailed t-test across populations:\n');
    fprintf('  t(%d) = %.3f, p = %.4g\n', ...
        stats.df,stats.tstat,p);

    % Also give nonparametric result
    pSR = signrank(s,0,'tail','right');

    fprintf('  Wilcoxon signed-rank: p = %.4g\n',pSR);

end


%% ============================================================
% MIXED-EFFECTS MODEL ACROSS ALL 18 SLOPES
% Population treated as random effect
% ============================================================

Slope = allSlopes(:);

Population = repmat((1:nPop)',nNet,1);

Network = [ ...
    repmat({'AlexNet'},nPop,1); ...
    repmat({'ResNet-50'},nPop,1); ...
    repmat({'ViT-L'},nPop,1)];

T = table( ...
    Slope, ...
    categorical(Population), ...
    categorical(Network), ...
    'VariableNames',{'Slope','Population','Network'});

T = rmmissing(T);

lme = fitlme(T, ...
    'Slope ~ 1 + Network + (1|Population)');

fprintf('\n============================================================\n');
fprintf('MIXED-EFFECTS MODEL ACROSS ALL POPULATIONS AND DNNs\n');
fprintf('Slope ~ Network + (1|Population)\n');
fprintf('============================================================\n');

disp(lme);

fprintf('\nANOVA:\n');
disp(anova(lme));


%% Overall mean slope across networks for each population
% This is a particularly clean population-level test.

populationMeanSlope = mean(allSlopes,2,'omitnan');

fprintf('\n============================================================\n');
fprintf('POPULATION-LEVEL MEAN SLOPE ACROSS THREE DNNs\n');
fprintf('============================================================\n');

for pop = 1:nPop
    fprintf('Population %d: %.3f ms\n', ...
        pop,populationMeanSlope(pop));
end

[~,pOverall,ciOverall,statsOverall] = ...
    ttest(populationMeanSlope,0,'Tail','right');

fprintf('\nMean +/- SEM = %.3f +/- %.3f ms\n', ...
    mean(populationMeanSlope), ...
    std(populationMeanSlope)/sqrt(nPop));

fprintf('One-tailed t-test:\n');
fprintf('t(%d) = %.3f, p = %.4g\n', ...
    statsOverall.df,statsOverall.tstat,pOverall);

fprintf('95%% CI = [%.3f %.3f]\n', ...
    ciOverall(1),ciOverall(2));

pOverallSR = signrank(populationMeanSlope,0,'tail','right');

fprintf('Wilcoxon signed-rank: p = %.4g\n',pOverallSR);