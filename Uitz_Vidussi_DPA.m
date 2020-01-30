% performs Diagnostic Pigment Analysis according to Uitz et al., 2006,
% JGR-Oceans: https://doi.org/10.1029/2005JC003207
% OR
% Vidussi et al, 2001, JGR-Oceans: https://doi.org/10.1029/1999JC000308
% please cite the above paper(s) if you use the corresponding DPA
% approaches below

% function written by Dylan Catlett

%% inputs:

% (1) [pigs_in] a numeric matrix of phytoplankton pigment concentrations, 
% where columns are variables (pigments) and rows are observations 
% (samples)

% (2) [varz] a cell array of variable names corresponding to each column of 
% pigs_in

% (3) [weight_scheme] a string indicating whether the weighting scheme
% you'd like to use in your DPA. Options are: 'Uitz' or 'local'; 
% ** use 'Uitz' for the multiple-linear-regression-derived coefficients 
% provided by Uitz 
% ** use 'Vidussi' to follow the method of Vidussi et al., 2001 (see above
% citation). Note that if the standard DPA pigment suite is altered by
% specifying pigments in leave_me_out or put_me_in, these pigments are
% included in the calculation of biomass proportions and tchla associated
% with each PFT without warning. 
% ** use 'local' to run an MLR on your data set and use the derived
% coefficients as weights for your DPA

% (4) [leave_me_out] a cell array of strings indicating the name or 
% abbreviation of pigments you'd like to remove from the standard DPA 
% analysis. provide an empty cell array is you do not wish to leave any of
% the standard diagnostic pigments out.

% (5) [put_me_in] a cell array of strings indicating the mames of any other
% pigments you'd like to add to the DPA. provide an empty cell array is you
% do not wish to add any pigments to the analysis.

% (6) [what_tot] a string indicating what indicator of biomass you want to
% use. in other words, the response variable in the multiple linear 
% regression employed by the Uitz DPA approach. the standard and default 
% is to use total chlorophyll a. observations of this variable should be 
% included in [pigs_in] and the name of the string you supply should be
% found in [varz]

%% outputs:

% (1) [varz_out] a cell array of strings indicating the order of the output 
% variables in pft_chl and pft_frac

% (2) [pft_data] a matrix of the TChla (same units as input pigment
% concentrations) and the TChla fraction (%) associated with each PFT; 
% rows are in the same order as those provided, columns are in the order of
% the [varz_out] cell array

% (3) [weight_out] a vector of the coefficients used to derive
% contributions of each PFT to TChla. This is a vector of 1's if 'Vidussi'
% is specified as the [weight_scheme]

% (4) [linmdl] 
% ** if weight_scheme == 'local', the struct w/ linear model statistics 
% returned by matlab's fitlm function
% ** if weight_scheme == 'Uitz', the struct w/ linear model statistics 
% returned by matlab's fitlm function using the sum of the considered 
% diagnostic pigments as the predictor and [what_tot] as the response
% variable
% ** if weight_scheme == 'Vidussi', this is 


function [varz_out,pft_data,weight_out,linmdl] = Uitz_Vidussi_DPA(pigs_in,varz,weight_scheme,leave_me_out,put_me_in,what_tot)


%% if weight_scheme == Uitz & put_me_in or leave_me_out exist, break...
if isequal(weight_scheme, 'Uitz') == 1 && (isempty(leave_me_out) == 0 || isempty(put_me_in) == 0)
    error(['You cannot add or remove pigments in DPA when using coefficients from Uitz et al., 2006. ' ...
        'Please set weight_scheme to ''local'' or remove ''leave_me_out'' and/or ''put_me_in''']);
else
    % do nothing
end

%% designate "diagnostic pigments" based on Uitz and user-specified variables:

% default pigments, taking into account variations:
pigs2get = {'Fuco', 'Perid', 'Hex', 'But', 'Allo', 'TChlb', 'Zea'};
pigs2get2 = {'Fucoxanthin', 'Peridinin', '19''-hexanoyloxyfucoxanthin',...
    '19''-butanoyloxyfucoxanthin', 'Alloxanthin', 'Total_chlorophyll_b', 'Zeaxanthin'};
pigs2get3 = {'Fuco', 'Perid', 'Hexfuco', 'Butfuco', 'Allo', 'TChlb', 'Zea'};
pigs2get4 = {'Fuco', 'Perid', 'HexFuco', 'ButFuco', 'Allo', 'TChlb', 'Zea'};

% pigments to remove specified by user
if isempty(leave_me_out) == 1
    % do nothing...
elseif isempty(leave_me_out) == 0
    % remove names in leave_me_out from pigs2get
    for i = 1:length(leave_me_out)
        ii = contains(pigs2get, leave_me_out);
        if sum(ii) == 0
            ii = contains(pigs2get2, leave_me_out);
            if sum(ii) == 0
                ii = contains(pigs2get3, leave_me_out);
                if sum(ii) == 0
                    ii = contains(pigs2get4, leave_me_out);
                end
            end
        end
        pigs2get(ii) = [];
        pigs2get2(ii) = [];
        pigs2get3(ii) = [];
        pigs2get4(ii) = [];
    end
end
% addtional DP's to consider from user:
if isempty(put_me_in) == 1
    % do nothing...
elseif isempty(put_me_in) == 0
    % add names in put_me_in to pigs2get
    if size(put_me_in,1) > size(put_me_in,2)
        % it's a row vector, so transpose first:
        pigs2get = cat(2, pigs2get, put_me_in');
        pigs2get2 = cat(2, pigs2get2, put_me_in');
        pigs2get3 = cat(2, pigs2get3, put_me_in');
        pigs2get4 = cat(2, pigs2get4, put_me_in');
    elseif size(put_me_in,1) <= size(put_me_in,2)
        pigs2get = cat(2, pigs2get, put_me_in);
        pigs2get2 = cat(2, pigs2get2, put_me_in);
        pigs2get3 = cat(2, pigs2get3, put_me_in);
        pigs2get4 = cat(2, pigs2get4, put_me_in);
    end
end

%% find the DP's in the supplied data matrix:
for i = 1:length(pigs2get)
    idx = find(ismember(varz,pigs2get{i}) == 1 | ismember(varz,pigs2get2{i}) == 1 ...
        | ismember(varz,pigs2get3{i}) == 1 | ismember(varz,pigs2get4{i}) == 1);
    if idx == 0
        error(['We cant find ', pigs2get2{i}, 'in the variable names you supplied. Please check the format of your pigment data and variable names and try again'])
    elseif idx > 0
        pigcols(i) = idx;
    end
    
    % also find total chlorophyll a column
    idx = find(ismember(varz,what_tot) == 1);
    totzcol = idx;
    
end

%% do the linear modeling (or application of Uitz coefficients and record stats and data:
if isequal(weight_scheme, 'Uitz') == 1
    weight_out = [1.41 1.41 1.27 0.35 0.6 1.01 0.86];
    % ^order = Fuco, perid, hex, but, allo, tchlb, zea
    DPsum = sum(pigs_in(:,pigcols) .* weight_out,2);
    pft_frac = (weight_out .* pigs_in(:,pigcols)) ./ DPsum;
    pft_chl = pft_frac .* pigs_in(:,totzcol);
    % calculate the linear model for DP sum and total chlorophyll:
    linmdl = fitlm(DPsum,pigs_in(:,totzcol),'linear','Intercept',false);
elseif isequal(weight_scheme, 'local') == 1
    linmdl = fitlm(pigs_in(:,pigcols),pigs_in(:,totzcol),'linear','Intercept',false);
    weight_out = linmdl.Coefficients.Estimate';
    DPsum = linmdl.Fitted;
    % get pft fractions and total chlorophyll:
    pft_frac = (weight_out .* pigs_in(:,pigcols)) ./ DPsum;
    pft_chl = pft_frac .* pigs_in(:,totzcol);
elseif isequal(weight_scheme, 'Vidussi') == 1
    linmdl = NaN;
    weight_out = ones(1, length(pigcols));
    DPsum = sum(pigs_in(:,pigcols),2);
    % get pft fractions and total chlorophyll:
    pft_frac = pigs_in(:,pigcols) ./ DPsum;
    pft_chl = pft_frac .* pigs_in(:,totzcol);
end

%% combine output data and make output variable names accordingly
pft_data = [pft_frac,pft_chl];
% set your variable names:
varz_out = cat(2,pigs2get,pigs2get);
for i = 1:length(varz_out)
    if i <= size(pft_frac,2)
        varz_out{i} = [varz_out{i},'_frac'];
    elseif i > size(pft_frac,2)
        varz_out{i} = [varz_out{i},'_chl'];
    end
end

