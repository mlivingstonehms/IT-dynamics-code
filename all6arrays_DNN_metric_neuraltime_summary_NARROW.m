%% all6arrays_DNN_metric_neuraltime_summary.m
% Six independent arrays, using Red 220227 only.
% For each population x DNN:
%   1) linearly normalize neural temporal centroids to 0..1
%   2) map layer depth to 0..1
%   3) interpolate both centroid and metric-RSA trajectories onto a common
%      depth grid for that DNN (AlexNet 22, ResNet 17, ViT-L 24)
%   4) average across populations
% Five panels show the mean trajectories for the five image metrics.
% Marker shape = DNN; marker color = normalized DNN depth.

close all
clear all

resultFolder = ...
    '/n/data2/hms/neurobio/livingstone/marge/figimages/temp';

baby1Folder = ...
    '/n/data2/hms/neurobio/livingstone/marge/margemonkeys/baby1/data';

files = { ...
    '/n/data2/hms/neurobio/livingstone/marge/margemonkeys/red/data/red240302_DNN_METRIC_NEURALTIME_bridge_results.mat', ...
    '/n/data2/hms/neurobio/livingstone/marge/margemonkeys/chatmetricsEdges/data/paul230204_DNN_METRIC_NEURALTIME_bridge_results.mat', ...
    '/n/data2/hms/neurobio/livingstone/marge/margemonkeys/chatmetricsEdges/data/paul220428_DNN_METRIC_NEURALTIME_bridge_results.mat', ...
    '/n/data2/hms/neurobio/livingstone/marge/margemonkeys/chatmetricsEdges/data/baby5221104_DNN_METRIC_NEURALTIME_bridge_results.mat', ...
    '/n/data2/hms/neurobio/livingstone/marge/margemonkeys/chatmetricsEdges/data/baby1230629_DNN_METRIC_NEURALTIME_bridge_results.mat', ...
    '/n/data2/hms/neurobio/livingstone/marge/margemonkeys/chatmetricsEdges/data/baby1220226_DNN_METRIC_NEURALTIME_bridge_results.mat'};

nPop = numel(files);
nNet = 3;
nMetric = 5;

netNames = {'AlexNet','ResNet-50','ViT-L'};
targetN = [22 17 24];
markerTypes = {'o','s','^'};

metricNames = { ...
     'Area fraction', ...
     'Colorfulness', ...
     'Spatial frequency', ...
    'Number of parts', ...
     'Size-normalized perimeter'};

   metricOrder = [4 2 5 1 3];
    

% Storage: one cell per DNN
centroidAll = cell(nNet,1);
metricAll   = cell(nNet,1);

for net = 1:nNet
    centroidAll{net} = nan(nPop,targetN(net));
    metricAll{net}   = nan(nPop,nMetric,targetN(net));
end

for pop = 1:nPop
    S = load(files{pop}, ...
        'centroidByNetwork','metricRSAByNetwork');

    for net = 1:nNet
        c = double(S.centroidByNetwork{net}(:));
       M = double(S.metricRSAByNetwork{net});

if size(M,1) ~= nMetric
    M = M';
end

% Reorder from original order:
% parts, color, perimeter, area, SF
% to:
% area, color, SF, parts, perimeter
M = M([4 2 5 1 3],:);

        % Relative neural time within this population x DNN
        cNorm = (c-min(c)) ./ (max(c)-min(c));

        % Normalize layer depth so arrays with different layer counts can
        % still be averaged (important for Red 220227 ResNet).
        srcDepth = linspace(0,1,numel(c));
        dstDepth = linspace(0,1,targetN(net));

        centroidAll{net}(pop,:) = ...
            interp1(srcDepth,cNorm,dstDepth,'linear');

        for m = 1:nMetric
            metricAll{net}(pop,m,:) = ...
                interp1(srcDepth,M(m,:),dstDepth,'linear');
        end
    end
end

meanCentroid = cell(nNet,1);
meanMetric   = cell(nNet,1);

for net = 1:nNet
    meanCentroid{net} = mean(centroidAll{net},1,'omitnan');
    meanMetric{net}   = squeeze(mean(metricAll{net},1,'omitnan'));
end

% Blue -> cyan -> gold -> orange depth colormap
depthMap = interp1( ...
    [0 .35 .70 1], ...
    [0.14 0.34 0.90; ...
     0.07 0.75 0.78; ...
     0.95 0.69 0.20; ...
     0.90 0.42 0.12], ...
    linspace(0,1,256));

fig = figure('Color','w','Position',[80 60 1650 900]);
tl = tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');

for m = 1:nMetric
    ax = nexttile(tl,m);
    hold(ax,'on');

    for net = 1:nNet
        d = linspace(0,1,targetN(net));
        x = meanMetric{net}(m,:);
        y = meanCentroid{net};

        % Connect successive layers with depth-colored line segments.
        % This is purely a display change; no smoothing is applied.
        for k = 1:targetN(net)-1
            ci = 1 + round(mean(d(k:k+1))*255);
            plot(ax,x(k:k+1),y(k:k+1),'-', ...
                'Color',depthMap(ci,:), ...
                'LineWidth',2, ...
                'HandleVisibility','off');
        end

        % Larger symbols at every actual averaged layer-depth point
        for k = 1:targetN(net)
            ci = 1 + round(d(k)*255);
            plot(ax,x(k),y(k),markerTypes{net}, ...
                'MarkerSize',9.5, ...
                'MarkerFaceColor',depthMap(ci,:), ...
                'MarkerEdgeColor',depthMap(ci,:), ...
                'LineWidth',1.2, ...
                'LineStyle','none');
        end
    end

    ylim(ax,[-.04 1.04])
    xlabel(ax,'Mean metric RSA across 6 arrays')
    ylabel(ax,'Relative neural time (0 = early, 1 = late)')
    title(ax,metricNames{m},'FontWeight','bold')
    box(ax,'on')
    set(ax,'LineWidth',2,'FontSize',12,'TickDir','out')
end

%% ================================================================
% LEGEND / DEPTH SCALE IN SIXTH TILE
% ================================================================

axLeg = nexttile(tl,6);
axis(axLeg,'off')
hold(axLeg,'on')

% DNN architecture symbols
h1 = plot(axLeg,nan,nan,'o', ...
    'MarkerSize',10, ...
    'MarkerFaceColor',[0.35 0.35 0.35], ...
    'MarkerEdgeColor',[0.35 0.35 0.35], ...
    'LineStyle','none');

h2 = plot(axLeg,nan,nan,'s', ...
    'MarkerSize',10, ...
    'MarkerFaceColor',[0.35 0.35 0.35], ...
    'MarkerEdgeColor',[0.35 0.35 0.35], ...
    'LineStyle','none');

h3 = plot(axLeg,nan,nan,'^', ...
    'MarkerSize',10, ...
    'MarkerFaceColor',[0.35 0.35 0.35], ...
    'MarkerEdgeColor',[0.35 0.35 0.35], ...
    'LineStyle','none');

lgd = legend(axLeg,[h1 h2 h3], ...
    {'AlexNet','ResNet-50','ViT-L'});

set(lgd, ...
    'Location','northwest', ...
    'Box','off', ...
    'FontSize',12);

text(axLeg,0.02,0.58, ...
    sprintf(['Each trajectory follows DNN depth from shallow to deep.\n\n' ...
    'Temporal centroids are normalized within each array x DNN\n' ...
    'before averaging across the 6 independent arrays.']), ...
    'Units','normalized', ...
    'FontSize',11, ...
    'VerticalAlignment','top');

% title(tl, ...
%     'Image-metric correspondence weakens as DNN correspondence shifts later in neural time', ...
%     'FontWeight','bold', ...
%     'FontSize',16);


%% ================================================================
% COMPACT HORIZONTAL DEPTH SCALE INSIDE SIXTH TILE
% Avoid colorbar(tiledlayout,...) for MATLAB-version compatibility.
% ================================================================

colormap(fig,depthMap)

% Make a tiny axes manually in the sixth-tile region.
cbAx = axes(fig, ...
    'Position',[0.72 0.105 0.18 0.018], ...
    'Visible','off');

% Draw the blue-to-orange gradient directly.
gradientImage = repmat(linspace(0,1,256),2,1);

imagesc(cbAx,[0 1],[0 1],gradientImage);
colormap(cbAx,depthMap);

set(cbAx, ...
    'Visible','on', ...
    'YTick',[], ...
    'XTick',[0 1], ...
    'XTickLabel',{'shallow','deep'}, ...
    'FontSize',11, ...
    'LineWidth',1.5, ...
    'TickDir','out', ...
    'Box','on');

xlabel(cbAx,'DNN layer depth', ...
    'FontSize',11, ...
    'FontWeight','bold');

%% ================================================================
% SAVE PUBLICATION FIGURE
% ================================================================

saveFolder = ...
    '/n/data2/hms/neurobio/livingstone/marge/figimages/temp';

if ~exist(saveFolder,'dir')
    mkdir(saveFolder);
end

% High-resolution PNG
print(fig, ...
    fullfile(saveFolder,'centroidsxtimexDNN_CURRENT_METRICS.png'), ...
    '-dpng','-r600');

% Vector PDF
set(fig,'PaperPositionMode','auto');

print(fig, ...
    fullfile(saveFolder,'centroidsxtimexDNN_CURRENT_METRICS.pdf'), ...
    '-dpdf','-painters','-bestfit');

fprintf('\nSaved PNG and PDF versions of figure.\n');


%% ================================================================
% SAVE FIGURE
% ================================================================

saveFolder = ...
    '/n/data2/hms/neurobio/livingstone/marge/figimages/temp';

if ~exist(saveFolder,'dir')
    mkdir(saveFolder);
end

print(fig, ...
    fullfile(saveFolder,'centroidsxtimexDNN_CURRENT_METRICS.png'), ...
    '-dpng','-r300');

fprintf('\nSaved figure:\n%s\n', ...
    fullfile(saveFolder,'centroidsxtimexDNN_CURRENT_METRICS.png'));



%% ========================================================================
% SECOND VERSION OF ORIGINAL FIGURE WITH COMMON X-AXIS
% ========================================================================

% The original figure is still open as fig.
% Make a complete copy so NOTHING about the original is changed.

figCommon = copyobj(fig,0);

set(figCommon, ...
    'Name','Common x-axis version', ...
    'NumberTitle','off');


%% Find all plotting axes in the copied figure

axCommon = findall(figCommon,'Type','axes');

% Remove any legend/colorbar axes if present
axCommon = axCommon(~strcmp(get(axCommon,'Tag'),'legend'));
axCommon = axCommon(~strcmp(get(axCommon,'Tag'),'Colorbar'));


%% Find global x range from the DATA ALREADY PLOTTED

globalXmax = -Inf;
globalXmin = Inf;

for a = 1:numel(axCommon)

    children = axCommon(a).Children;

    for c = 1:numel(children)

        if isprop(children(c),'XData')

            xd = children(c).XData;
            xd = xd(isfinite(xd));

            if ~isempty(xd)
                globalXmin = min(globalXmin,min(xd));
                globalXmax = max(globalXmax,max(xd));
            end
        end
    end
end


%% Use common scale

% Since RSA values are positive here, start at zero
commonXmin = 0;

% Round upper limit upward to nearest 0.1
commonXmax = 0.76;

fprintf('\nCommon x-axis limits = %.2f to %.2f\n', ...
    commonXmin,commonXmax);


%% Apply identical x limits to every panel

for a = 1:numel(axCommon)
    xlim(axCommon(a),[commonXmin commonXmax]);
end
for a = 1:numel(axCommon)
    xlim(axCommon(a),[0 0.7]);
    xticks(axCommon(a),0:0.1:0.7);
end

%% Save as DIFFERENT files

exportgraphics( ...
    figCommon, ...
    fullfile(saveFolder, ...
    'centroidsxtimexDNN_CURRENT_METRICS_COMMON_X.png'), ...
    'Resolution',300);

exportgraphics( ...
    figCommon, ...
    fullfile(saveFolder, ...
    'centroidsxtimexDNN_CURRENT_METRICS_COMMON_X.pdf'), ...
    'ContentType','vector');




%% ========================================================================
% QUANTIFY DNN-DEPTH DEPENDENCE OF EACH IMAGE METRIC
%
% x = normalized DNN layer depth (0 = shallow, 1 = deep)
% y = mean DNN layer-image metric RSA across 6 populations
%
% Negative slope = metric correspondence decreases with DNN depth
% ========================================================================

depthSlope = nan(nNet,nMetric);
depthIntercept = nan(nNet,nMetric);
depthRho = nan(nNet,nMetric);
depthP = nan(nNet,nMetric);

for net = 1:nNet

    depth = linspace(0,1,targetN(net))';

    for m = 1:nMetric

        rsa = meanMetric{net}(m,:)';

        good = isfinite(depth) & isfinite(rsa);

        % Linear slope: RSA change from shallow -> deep
        p = polyfit(depth(good),rsa(good),1);

        depthSlope(net,m) = p(1);
        depthIntercept(net,m) = p(2);

        % Rank-based monotonic relationship with depth
        [depthRho(net,m),depthP(net,m)] = ...
            corr(depth(good),rsa(good), ...
            'Type','Spearman', ...
            'Rows','complete');
    end
end


%% Display results

fprintf('\n\n============================================================\n');
fprintf('DNN DEPTH -> IMAGE-METRIC RSA\n');
fprintf('Negative slope = stronger correspondence in shallow layers\n');
fprintf('============================================================\n');

for m = 1:nMetric

    fprintf('\n%s\n',metricNames{m});

    for net = 1:nNet

        fprintf('  %-10s slope = % .4f   rho = % .3f   p = %.4g\n', ...
            netNames{net}, ...
            depthSlope(net,m), ...
            depthRho(net,m), ...
            depthP(net,m));
    end
end


%% ========================================================================
% SUMMARY FIGURE: ONE SLOPE PER DNN FOR EACH METRIC
% ========================================================================

figSlope = figure( ...
    'Color','w', ...
    'Position',[200 150 900 500]);

hold on

xpos = 1:nMetric;

offsets = [-0.18 0 0.18];

for net = 1:nNet

    plot( ...
        xpos + offsets(net), ...
        depthSlope(net,:), ...
        markerTypes{net}, ...
        'MarkerSize',10, ...
        'MarkerFaceColor',[0.4 0.4 0.4], ...
        'MarkerEdgeColor',[0.4 0.4 0.4], ...
        'LineStyle','none', ...
        'LineWidth',1.5);
end

yline(0,'k--','LineWidth',1.2);

xlim([0.5 nMetric+0.5]);

xticks(1:nMetric);
xticklabels(metricNames);
xtickangle(25);

ylabel('Slope of metric RSA vs DNN depth');

legend(netNames, ...
    'Location','best', ...
    'Box','off');

set(gca, ...
    'FontSize',12, ...
    'LineWidth',1.5, ...
    'TickDir','out');

box off

title('Change in image-metric correspondence across DNN depth');


%% Save

exportgraphics( ...
    figSlope, ...
    fullfile(saveFolder, ...
    'DNN_depth_metricRSA_slopes.png'), ...
    'Resolution',300);



%% ========================================================================
% POPULATION-LEVEL DNN DEPTH SLOPES + MIXED-EFFECTS STATISTICS
%
% One slope for every:
%       population x DNN x image metric
%
% slope < 0:
% metric correspondence is stronger in shallow DNN layers
%
% We then test whether depth slopes differ among image metrics while
% accounting for repeated observations from populations and DNNs.
% ========================================================================

nPop = 6;

%% ------------------------------------------------------------
% Calculate one depth slope per population x DNN x metric
% ------------------------------------------------------------

nPop = 6;

popDepthSlope = nan(nPop,nNet,nMetric);

for pop = 1:nPop

    for net = 1:nNet

        % normalized DNN depth: 0 = shallow, 1 = deep
        depth = linspace(0,1,targetN(net))';

        for m = 1:nMetric

            % metricAll is a CELL ARRAY:
            % metricAll{net} = population x metric x DNN depth
            rsa = squeeze(metricAll{net}(pop,m,:));
            rsa = rsa(:);

            good = isfinite(depth) & isfinite(rsa);

            if sum(good) >= 3
                p = polyfit(depth(good),rsa(good),1);
                popDepthSlope(pop,net,m) = p(1);
            end

        end
    end
end


%% ========================================================================
% Put everything into a table
% ========================================================================

Slope = [];
Metric = strings(0,1);
Network = strings(0,1);
Population = strings(0,1);

for pop = 1:nPop

    for net = 1:nNet

        for m = 1:nMetric

            thisSlope = popDepthSlope(pop,net,m);

            if isfinite(thisSlope)

                Slope(end+1,1) = thisSlope;

                Metric(end+1,1) = string(metricNames{m});
                Network(end+1,1) = string(netNames{net});
                Population(end+1,1) = "Pop" + pop;

            end
        end
    end
end

T = table( ...
    Slope, ...
    categorical(Metric), ...
    categorical(Network), ...
    categorical(Population), ...
    'VariableNames', ...
    {'Slope','Metric','Network','Population'});


%% ========================================================================
% MIXED-EFFECTS MODEL
%
% Main question:
% Do the depth slopes differ among image metrics?
%
% Population is the biological replicate.
% Network is included as a fixed effect because there are only 3 networks.
% Random intercept accounts for repeated measures within population.
% ========================================================================

lme = fitlme(T, ...
    'Slope ~ Metric + Network + (1|Population)');

disp(lme)


%% Overall test of fixed effects

anovaResults = anova(lme);

fprintf('\n============================================================\n');
fprintf('MIXED-EFFECTS MODEL\n');
fprintf('Slope ~ Metric + Network + (1|Population)\n');
fprintf('============================================================\n');

disp(anovaResults)


%% ========================================================================
% DESCRIPTIVE MEANS ACROSS POPULATIONS AND NETWORKS
% ========================================================================

fprintf('\n============================================================\n');
fprintf('DEPTH SLOPES BY METRIC\n');
fprintf('Mean +/- SEM across 6 populations x 3 networks\n');
fprintf('============================================================\n');

metricMean = nan(nMetric,1);
metricSEM  = nan(nMetric,1);

for m = 1:nMetric

    vals = popDepthSlope(:,:,m);
    vals = vals(:);
    vals = vals(isfinite(vals));

    metricMean(m) = mean(vals);
    metricSEM(m) = std(vals) / sqrt(numel(vals));

    fprintf('%-28s  %.4f +/- %.4f\n', ...
        metricNames{m}, ...
        metricMean(m), ...
        metricSEM(m));
end


%% ========================================================================
% PLANNED CONTRAST
%
% Hypothesis:
%
% "shallow-feature" metrics:
%       Colorfulness
%       Area fraction
%       Spatial frequency
%
% versus
%
% "structural" metrics:
%       Number of parts
%       Size-normalized perimeter
%
% More negative slope = stronger shallow-layer bias
% ========================================================================
% NEW order:
% 1 area
% 2 color
% 3 SF
% 4 parts
% 5 normalized perimeter

shallowMetrics = [1 2 3];   % area, color, SF
structMetrics  = [4 5];     % parts, normalized perimeter

shallowSlope = nan(nPop,nNet);
structSlope  = nan(nPop,nNet);

for pop = 1:nPop
    for net = 1:nNet

        shallowSlope(pop,net) = ...
            mean(popDepthSlope(pop,net,shallowMetrics),3,'omitnan');

        structSlope(pop,net) = ...
            mean(popDepthSlope(pop,net,structMetrics),3,'omitnan');

    end
end


%% Difference within each population x network

differenceSlope = shallowSlope - structSlope;

fprintf('\n============================================================\n');
fprintf('PLANNED CONTRAST\n');
fprintf('Color + Area + SF versus Parts + Normalized Perimeter\n');
fprintf('Negative difference = stronger shallow bias for Color/Area/SF\n');
fprintf('============================================================\n');

fprintf('Mean shallow-feature slope = %.4f\n', ...
    mean(shallowSlope(:),'omitnan'));

fprintf('Mean structural slope      = %.4f\n', ...
    mean(structSlope(:),'omitnan'));

fprintf('Mean difference            = %.4f\n', ...
    mean(differenceSlope(:),'omitnan'));


%% ========================================================================
% Mixed model specifically for planned contrast
% ========================================================================

ContrastSlope = [];
Group = strings(0,1);
Network2 = strings(0,1);
Population2 = strings(0,1);

for pop = 1:nPop

    for net = 1:nNet

        ContrastSlope(end+1,1) = shallowSlope(pop,net);
        Group(end+1,1) = "Color_Area_SF";
        Network2(end+1,1) = string(netNames{net});
        Population2(end+1,1) = "Pop" + pop;

        ContrastSlope(end+1,1) = structSlope(pop,net);
        Group(end+1,1) = "Parts_Perimeter";
        Network2(end+1,1) = string(netNames{net});
        Population2(end+1,1) = "Pop" + pop;

    end
end

TC = table( ...
    ContrastSlope, ...
    categorical(Group), ...
    categorical(Network2), ...
    categorical(Population2), ...
    'VariableNames', ...
    {'Slope','Group','Network','Population'});

lmeContrast = fitlme(TC, ...
    'Slope ~ Group + Network + (1|Population)');

fprintf('\nPLANNED-CONTRAST MIXED MODEL\n');
disp(lmeContrast);

fprintf('\nANOVA:\n');
disp(anova(lmeContrast));


%% ========================================================================
% FIGURE: population-level slopes
% ========================================================================

figure('Color','w','Position',[200 150 950 550]);
hold on

% individual population x network observations
for m = 1:nMetric

    vals = popDepthSlope(:,:,m);
    vals = vals(:);
    vals = vals(isfinite(vals));

    % jitter slightly just so individual observations can be seen
    xj = m + linspace(-0.12,0.12,numel(vals))';

    plot(xj,vals,'o', ...
        'MarkerSize',4, ...
        'MarkerEdgeColor',[0.65 0.65 0.65], ...
        'MarkerFaceColor',[0.65 0.65 0.65], ...
        'LineStyle','none');

end

% means +/- SEM
errorbar(1:nMetric, ...
    metricMean,metricSEM, ...
    'ko', ...
    'MarkerFaceColor','k', ...
    'MarkerSize',8, ...
    'LineWidth',2, ...
    'CapSize',12);

yline(0,'k--');

xlim([0.5 nMetric+0.5]);

xticks(1:nMetric);
xticklabels(metricNames);
xtickangle(25);

ylabel('Slope of metric RSA vs DNN depth');

title('Image-metric correspondence across DNN depth');

set(gca, ...
    'FontSize',12, ...
    'LineWidth',1.5, ...
    'TickDir','out');

box off
%% ========================================================================
% CONSERVATIVE POPULATION-LEVEL TEST
%
% One independent contrast value per neural population.
% First average across DNN architectures within each population.
%
% Negative = Area/Color/SF have stronger shallow-layer bias
%            than Parts/Perimeter
% ========================================================================

popContrast = nan(nPop,1);

for pop = 1:nPop

    % Mean across the 3 DNN architectures
    meanAcrossDNN = squeeze( ...
        mean(popDepthSlope(pop,:,:),2,'omitnan') );
    
    % meanAcrossDNN is now metric x 1

    shallowMean = mean( ...
        meanAcrossDNN(shallowMetrics), ...
        'omitnan');

    structuralMean = mean( ...
        meanAcrossDNN(structMetrics), ...
        'omitnan');

    popContrast(pop) = shallowMean - structuralMean;
end


fprintf('\n============================================================\n');
fprintf('POPULATION-LEVEL CONTRAST\n');
fprintf('Area/Color/SF minus Parts/Perimeter\n');
fprintf('Negative = stronger shallow bias for Area/Color/SF\n');
fprintf('============================================================\n');

for pop = 1:nPop
    fprintf('Population %d: % .4f\n',pop,popContrast(pop));
end

fprintf('\nMean difference = %.4f\n',mean(popContrast));
fprintf('SEM             = %.4f\n', ...
    std(popContrast)/sqrt(nPop));


%% Parametric test across the 6 independent populations

[~,p_t,ci_t,stats_t] = ttest(popContrast,0);

fprintf('\nOne-sample t-test across populations:\n');
fprintf('t(%d) = %.3f, p = %.6g\n', ...
    stats_t.df,stats_t.tstat,p_t);

fprintf('95%% CI = [%.4f %.4f]\n', ...
    ci_t(1),ci_t(2));


%% Nonparametric signed-rank test as robustness check

p_signrank = signrank(popContrast,0);

fprintf('\nWilcoxon signed-rank test:\n');
fprintf('p = %.6g\n',p_signrank);


%% Show the six independent observations

figure('Color','w','Position',[300 200 450 500]);

plot(ones(nPop,1),popContrast,'ko', ...
    'MarkerFaceColor','k', ...
    'MarkerSize',8);

hold on

yline(0,'k--');

errorbar(1, ...
    mean(popContrast), ...
    std(popContrast)/sqrt(nPop), ...
    'ro', ...
    'MarkerFaceColor','r', ...
    'MarkerSize',10, ...
    'LineWidth',2, ...
    'CapSize',15);

xlim([0.5 1.5]);
xticks(1);
xticklabels({'Area/Color/SF - Parts/Perimeter'});

ylabel('Difference in DNN-depth slope');

set(gca,'FontSize',12,'LineWidth',1.5,'TickDir','out');
box off

%% ========================================================================
% NARROW PUBLICATION FIGURE: DNN-DEPTH SLOPES WITH SEM
%
% One mean slope per DNN x metric, with SEM across the 6 populations.
% AlexNet = blue circles
% ResNet-50 = green squares
% ViT-L = red triangles
%
% This is a NEW figure and does not replace the earlier slope figure.
% ========================================================================

% Mean and SEM across the 6 independent populations, separately
% for each DNN and each metric.
slopeMeanByNet = nan(nNet,nMetric);
slopeSEMByNet  = nan(nNet,nMetric);

for net = 1:nNet
    for m = 1:nMetric

        vals = popDepthSlope(:,net,m);
        vals = vals(isfinite(vals));

        slopeMeanByNet(net,m) = mean(vals,'omitnan');
        slopeSEMByNet(net,m)  = std(vals,'omitnan') / sqrt(numel(vals));

    end
end

% Requested colors:
% AlexNet blue, ResNet-50 green, ViT-L red
dnnColors = [ ...
    0.00 0.35 1.00; ...   % AlexNet
    0.00 0.55 0.00; ...   % ResNet-50
    1.00 0.00 0.00];      % ViT-L

% Tall/narrow figure: width is about half the height
figSlopeNarrow = figure( ...
    'Color','w', ...
    'Position',[250 80 460 880]);

ax = axes(figSlopeNarrow);
hold(ax,'on');

xpos = 1:nMetric;

% Plot each DNN as one connected colored trajectory with SEM error bars
h = gobjects(nNet,1);

for net = 1:nNet

    % connecting line
    plot(ax, ...
        xpos, ...
        slopeMeanByNet(net,:), ...
        '-', ...
        'Color',dnnColors(net,:), ...
        'LineWidth',4, ...
        'HandleVisibility','off');

    % points + SEM
    h(net) = errorbar(ax, ...
        xpos, ...
        slopeMeanByNet(net,:), ...
        slopeSEMByNet(net,:), ...
        markerTypes{net}, ...
        'Color',dnnColors(net,:), ...
        'MarkerFaceColor',dnnColors(net,:), ...
        'MarkerEdgeColor',dnnColors(net,:), ...
        'MarkerSize',24, ...
        'LineWidth',4, ...
        'CapSize',8, ...
        'LineStyle','none');

end

yline(ax,0,'k--','LineWidth',4,'HandleVisibility','off');

xlim(ax,[0.6 nMetric+0.4]);
ylim(ax,[-0.70 0.02]);

xticks(ax,1:nMetric);
xticklabels(ax,metricNames);
xtickangle(ax,48);

ylabel(ax,'Slope of metric RSA vs DNN depth');

% title(ax, ...
%     {'Change in image-metric correspondence', ...
%      'across DNN depth'}, ...
%     'FontWeight','bold');

legend(ax,h,netNames, ...
    'Location','southeast', ...
    'Box','on');

set(ax, ...
    'FontSize',11, ...
    'LineWidth',4, ...
    'TickDir','out');

box(ax,'on');


%% Save narrow figure separately

exportgraphics( ...
    figSlopeNarrow, ...
    fullfile(saveFolder, ...
    'DNN_depth_metricRSA_slopes_NARROW_SEM.png'), ...
    'Resolution',600);

exportgraphics( ...
    figSlopeNarrow, ...
    fullfile(saveFolder, ...
    'DNN_depth_metricRSA_slopes_NARROW_SEM.pdf'), ...
    'ContentType','vector');

fprintf('\nSaved narrow DNN-depth slope figure with SEM:\n');
fprintf('%s\n', ...
    fullfile(saveFolder,'DNN_depth_metricRSA_slopes_NARROW_SEM.png'));

