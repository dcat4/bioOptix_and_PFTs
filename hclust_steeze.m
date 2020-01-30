
% this is a convenient function for doing heirarchical cluster analysis
% using a series of matlabs built in functions. makes a pretty plot of the
% results and can also return indices to work up the cluster results

% operational for getting the cluster tree plots; working on fleshing out other outputs


%% inputs:

% (1) [pigs_in] a numeric matrix of phytoplankton pigment concentrations, 
% where columns are variables (pigments) and rows are observations 
% (samples)

% (2) [varz] a cell array of variable names corresponding to each column of 
% pigs_in

% (3) [put_me_in] a cell array of strings indicating the mames of pigments
% (in [varz]) you'd like to consider in your cluster analysis. 

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

% (5) [clust_params] a cell array of strings indicating parameters
% (distance & linkage metrics) to use in cluster analysis
% ** one entry should include 'dist=X', where X is the word specifying
% your distance and should match one of those supported by matlab's pdist
% function or the Fathom toolbox's f_dis function
% ** one entry should include 'link=Y', where Y is the word specifying your
% linkage method and match one of those supported by matlab's linkage
% function.

%% outputs:

% (1) [tree_plot] a figure showing the results of the cluster analysis.

% (2) [clustidx] NOT IN THERE YET - eventually will want to return an array
% containing the clusters each variable belongs to.. with options to color
% the cluster tree accordingly 

function [tree_plot] = hclust_steeze(pigs_in,varz,put_me_in,stdize,clust_params)

%% find the pigments in the supplied data matrix:
for i = 1:length(put_me_in)
    idx = find(ismember(varz,put_me_in{i}) == 1);
    if idx == 0
        error(['We cant find ', put_me_in{i}, 'in the variable names you supplied. Please check the format of your pigment data and variable names and try again'])
    elseif idx > 0
        pigcols(i) = idx;
    end
    
end

pigs4clust = pigs_in(:,pigcols);
varz4clust = varz(pigcols);

%% preprocess your data for cluster calculation:

pigs4clust(any(isnan(pigs4clust),2),:) = []; % remove NaN's

if strcmp(stdize, 'none') == 1 || isempty(stdize) == 1
    % do nothing
elseif strcmp(stdize, 'Z') == 1
    % z score:
    pigs4clust = zscore(pigs4clust, 0, 1);
elseif strcmp(stdize, 'center') == 1
    % subtract column mean:
    pigs4clust = pigs4clust - mean(pigs4clust,1);
else
    error('Please supply [stdize] when calling this function');
end

%% below is copy/pasted from clustering scripts but needs editing:

% extract distance and linkage method from clust_params
for i = 1:length(clust_params)
    pp = clust_params{i};
    if contains(pp,'dist') == 1
        idx = strfind(pp,'=');
        pigdistmetric = pp(idx+1:end);
    elseif contains(pp,'link') == 1
        idx = strfind(pp,'=');
        piglinkmetric = pp(idx+1:end);
    end
end

if strcmp(pigdistmetric,'bc') == 1
    dd = f_dis(pigs4clust',pigdistmetric,0,1,1);    
else
    dd = pdist(pigs4clust',pigdistmetric);
    dd = squareform(dd);
end

ll = linkage(dd, piglinkmetric);

[tree_plot] = dendrogram(ll,'Labels',varz4clust); hold on; box on; 
set(gca,'XTickLabelRotation',45);
set(tree_plot,'linewidth',2, 'Color','k');
ylabel('Linkage Distance'); 
set(gca, 'FontSize', 14); 



