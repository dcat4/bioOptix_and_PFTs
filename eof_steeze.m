% this is a function for performing more convenient EOF analysis on
% phytoplankton pigment data. it wraps matlabs built-in pca function to do
% the EOF analysis, automates the manipulation of the data prior to
% calc'ing EOFs, and spits out a nice plot (maybe more in the future) for 
% analyzing results.

% function written by Dylan Catlett
% if used please cite:
% Catlett and Siegel, 2018, JGR-Oceans, 
% https://doi.org/10.1002/2017JC013195

% it's operational, eventually add correlation coefficients to the plots 
% you output, and possibly additional plots (ordination-style, plots of the
% amplitudes, etc...).

%% inputs:

% (1) [pigs_in] a numeric matrix of phytoplankton pigment concentrations, 
% where columns are variables (pigments) and rows are observations 
% (samples)

% (2) [varz] a cell array of variable names corresponding to each column of 
% pigs_in

% (3) [put_me_in] a cell array of strings indicating the mames of pigments
% (in [varz]) you'd like to consider in your EOFs. 

% (4) [stdize] a string indicating how you want your data manipulated prior
% to computing EOFs. 
% Options:
% ** 'none' == EOFs are computed from raw data. you can supply an empty
% array for the same effect
% ** 'Z' == data is Z-scored (each variable is scaled to 0 mean and unit
% variance; aka each variable has it's mean subtracted and is normalized to
% it's standard deviation) prior to computing EOF's
% ** 'center' == each variable has it's mean removed prior to computing
% EOF's

% (5) [n_pcs] a number indicating the number of modes you want returned

% (6) [modes2plot] a vector of the mode numbers you'd like plotted in
% [load_plot]

%% outputs:

% (1) [loadz] matrix of EOF loadings. each COLUMN??? is a mode, and col
% indices correspond to mode numbers. each ROW??? is a pigment and is in the
% order pigments provided in [varz]

% (2) [afz] a matrix of the EOF amplitudes (a.k.a. principal component scores).
% each column is a mode, and col indices correspond to mode numbers. each
% row is an observation and are in the same order as those provided in 
% [pigs_in]

% (3) [varexp] a vector containing the % variance explained by each
% more. indices again correspond to mode numbers (eg, var_by_mode(1) is the 
% variance explained by mode 1).

% (4) [load_plot] a figure showing the loadings for up to 6 modes specified 
% [modes2plot]

function [loadz,afz,varexp,load_plot] = eof_steeze(pigs_in,varz,put_me_in,stdize,n_pcs,modes2plot)

if isempty(put_me_in) == 1
    error('I cant do EOFs on no pigment data, dawg! Modify put_me_in please');
else
end

%% find the pigments in the supplied data matrix:
for i = 1:length(put_me_in)
    idx = find(ismember(varz,put_me_in{i}) == 1);
    if idx == 0
        error(['We cant find ', put_me_in{i}, 'in the variable names you supplied. Please check the format of your pigment data and variable names and try again'])
    elseif idx > 0
        pigcols(i) = idx;
    end
    
end

pigs4eof = pigs_in(:,pigcols);
varz4eof = varz(pigcols);

%% preprocess your data for EOF calculation:

pigs4eof(any(isnan(pigs4eof),2),:) = []; % remove NaN's

if strcmp(stdize, 'none') == 1 || isempty(stdize) == 1
    % do nothing
elseif strcmp(stdize, 'Z') == 1
    % z score:
    pigs4eof = zscore(pigs4eof, 0, 1);
elseif strcmp(stdize, 'center') == 1
    % subtract column mean:
    pigs4eof = pigs4eof - mean(pigs4eof,1);
else
    error('Please supply [stdize] when calling this function');
end

%% calc EOFs:

[loadz,afz,~,~,varexp] = pca(pigs4eof,'Centered',false,'NumComponents',n_pcs,'Rows','complete'); 

%% create your plot:

npanel = length(modes2plot);
nrowpanel = ceil(npanel/2);
ncolpanel = ceil(npanel/nrowpanel);

counter = 1;
load_plot = figure(); hold on; 
for i = 1:npanel
    subplot(nrowpanel,ncolpanel,counter); hold on; box on;
    title(['Mode ',num2str(modes2plot(i)),': ', num2str(varexp(modes2plot(i))),'%']);
    bar(1:size(loadz,1),loadz(:,modes2plot(i)),'k');
    set(gca,'YLim',[-1 1],'XTick',1:size(loadz,1),'XTickLabel',varz4eof,'XTickLabelRotation',45,...
        'FontSize',14);
    if rem(counter,2) ~= 0
        ylabel('EOF Loadings');
    else
        % do nothing
    end
    counter = counter+1;
end

