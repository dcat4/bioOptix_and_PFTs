% this is the aph derivative model cross-validation routine I developed.
% details in Catlett and Siegel, 2018, JGR-Oceans.

% you give it a set of aph derivative spectra (or any spectra),
% corresponding observations of a pft index that you want to model (HPLC
% pigment concentration or derived parameters), and specify some modeling
% parameters (see below), and it outputs a set of spectral coefficients 
% from the cross-validation exercise, intercepts (arbitrary, should be
% close to zero...), and goodness of fit statistics, and saves these to a
% .mat file with the name specified in the fcn call

%% input arguments:
% 1. daph: an array with all the spectra you want to use to train the model
% each spectrum/observation should be a row, and each daph(lambda) value a column
% (note that the model-training slices this data set so that each of the
% n_permutations blindly validates a model against an unknown validation
% data set)

% 2. pft: the corresponding values of the pft index you want to model with
% your input daph spectra. in other words, the order of observations/rows
% in pft should align with that in daph.

% 3. pft_index: specifies any constraints to apply to model outputs. 
% options: 
% (A)'pigment' --> model outputs are constrained to be >= 0 at
% each iteration of the model development. 
% (B) 'EOFs' --> model outputs are not constrained
% (C) 'compositions' --> model outputs are constrained to be >= 0 and <= 1

% 4. n_permutations: the number of times you want to do the entire
% cross-validation exercise. My go-to is 500

% 5. max_pcs: the maximum number of spectral principal components that can
% be used in model training. the model will test/use all principal components up
% to and including this number
% NOTE: because the model optimization chops up your data set, max_pcs needs to 
% be less than or equal to: 0.75 * (1-1/k) * X; where k is input # 6, X is the number
% of observations in your data set

% 6. k: the number of actual model trainings to do for each of the
% n_permutations cross-validation iterations. usually 5.

% 7. mdl_pick_metric: the gof statistic you want to optimize the model off
% of. Options:
% (A) 'R2' --> picks mdls with maximum R2
% (B) 'RMSE' --> picks mdls with minimum RMSE
% (C) 'avg' --> minimum mean % error
% (D) 'med' --> minimum median % error

% 8. output_file_name: what you want the output saved as (only available as
% .mat's right now).

%% output arguments:
% 1. coefficients: an array of model coefficients that is n_permutations x
% size(daph,2) [which should be the number of wavelengths you've have
% derivative values for]

% 2. intercepts: a vector of model intercepts (should be reasonably close
% to 0) that is 1xn_permutations and contains the correspong intercept for
% each set of coefficients in coefficients array

% 3. summary_gofs: a summary table of goodness of fit statistics across all
% n_permutations of model cross-validations

% 4. all_gofs: a struct with each field a g.o.f statistic and each entry an
% array detailing all statistics from each of n_permutation
% cross-validations

function [coefficients, intercepts, summary_gofs, all_gofs] = aphModelTrainCV(daph, pft, pft_index, n_permutations, max_pcs, k, mdl_pick_metric, output_file_name)

    %% basic QC...
    % check for NaN's in data and break if they are there...
    if sum(any(isnan(pft))) > 0
        disp('pft data has NaNs. remove them and try again plz');
        return
    elseif sum(any(isnan(daph))) > 0
        disp('spectral data has NaNs. remove them and try again plz');
        return
    end
    
    % check that pigment and aph data are the same number of rows
    if size(pft,1) ~= size(daph,1)
        disp('You dont have corresponding pft and spectral observations. fix that and try again');
        return
    end
    
    % check that the file output is a .mat
    if contains(output_file_name, '.mat') == 0
        disp('Im a robot who can only save .mat files. change your output file name please');
        return
    end
    
    % check that pft_index and mdl_pick_metric are specified properly
    if strcmp(pft_index, 'pigment') == 0 && strcmp(pft_index, 'EOFs') == 0 && strcmp(pft_index, 'compositions') == 0
        disp('I dont understand what youre trying to tell me with your pft_index');
        return
    elseif strcmp(mdl_pick_metric,'R2') == 0 && strcmp(mdl_pick_metric,'RMSE') == 0 && strcmp(mdl_pick_metric,'avg') == 0 && strcmp(mdl_pick_metric,'med') == 0  
        disp('I dont understand what youre trying to tell me with your mdl_pick_metric')
        return
    end

    max_components = max_pcs; % max derivative spectra components that can be used to formulate the model (standard is 100).
    spectra_4_mdl = daph;
                                                          
    % set the random number generator so you're random cross-validating is
    % reproducible...
    rng(1); 
    
    %% the model training loop:
    for i = 1:n_permutations
        % Create broad training data (75%) and validation data (25%)
        training_indices = randperm(length(pft),floor(length(pft) * 0.75))'; 
        pigs_training = pft(training_indices); % training data
        spectra_4_mdl_training = spectra_4_mdl(training_indices,:); % training spectra
        % validation data:
        pigs_validate = pft; 
        pigs_validate(training_indices) = []; % validation data
        spectra_4_mdl_validate = spectra_4_mdl; 
        spectra_4_mdl_validate(training_indices,:) = []; % validation spectra

        % Get set up for a k-fold cross validation:
        rand_ns = randperm(length(pigs_training));
        CV_indices = zeros(k,ceil(length(pigs_training)/k));
        n_leftovers = (length(pigs_training)/k - (floor(length(pigs_training)/k))) * k;
        counter_start = n_leftovers + 1;
        counter_end = n_leftovers + floor(length(pigs_training)/k);
        for j = 1:k
            CV_indices(j,1:floor(length(pigs_training)/k)) = rand_ns(counter_start:counter_end);
            counter_start = counter_start + floor(length(pigs_training)/k);
            counter_end = counter_end + floor(length(pigs_training)/k);
        end
        % add the leftovers to the CV_indices array, and put NaN's for the sets
        % where there's not enough leftovers to go around
        leftovers = rand_ns(1:n_leftovers);
        na_array = ones(1,k - length(leftovers)).*NaN;
        leftovers = [leftovers,na_array];
        CV_indices(:,length(CV_indices)) = leftovers;
        % so now CV indices contains k sets of random indices for k-fold CV

        % develop/optimize the model with k-fold cross-validation:
        for j = 1:k
            % split up CV data sets:
            these_CV_indices = CV_indices(j,:);
            these_CV_indices(isnan(these_CV_indices)) = [];
            CV_valid_pigs = pigs_training(these_CV_indices); % CV validation pigments
            CV_train_pigs = pigs_training;
            CV_train_pigs(these_CV_indices) = []; % CV training pigments
            CV_valid_spec = spectra_4_mdl_training(these_CV_indices,:); % CV validation spectra
            CV_train_spec = spectra_4_mdl_training;
            CV_train_spec(these_CV_indices,:) = []; % CV training spectra

            % Standardize spectra for pc's:
            CV_train_spec = (CV_train_spec - repmat(mean(CV_train_spec),size(CV_train_spec,1),1))./repmat(std(CV_train_spec),size(CV_train_spec,1),1);
            CV_valid_spec = (CV_valid_spec - repmat(mean(CV_valid_spec),size(CV_valid_spec,1),1))./repmat(std(CV_valid_spec),size(CV_valid_spec,1),1);

            % take principal components of this CV training set:
            [CV_EOFs_train,CV_AFs_train,~,~,~] = pca(CV_train_spec,'Centered',false,'NumComponents',max_components,'Rows','complete');
            % NOTE: "Centered" is set to false here because I've manually
            % centered and standardized the spectra just above

            % find the best possible model in terms of # components to use for
            % this CV:
            for l = 1:size(CV_AFs_train,2)
                the_lin_mdl = fitlm(CV_AFs_train(:,1:l), CV_train_pigs); % for an MLR model w/ derivative components

                these_betas = the_lin_mdl.Coefficients; % pull out linear coefficients of AF's
                these_betas = table2array(these_betas(:,1)); % convert data table to a matlab array
                this_alpha = these_betas(1); % pull out intercept
                these_betas(1) = []; % remove intercept so you only have multipliers in these_betas

                % now turn model coefficients for AF's into coefficients for the combined
                % derivative spectra:
                these_betas = CV_EOFs_train(:,1:l) * these_betas;

                % Model Validation:
                CV_modeled_pigs = (CV_valid_spec * these_betas) + this_alpha; 

                % constrain modeled pigments based on the pft_index you're
                % modeling
                if strcmp(pft_index, 'pigment') == 1
                    CV_modeled_pigs(CV_modeled_pigs < 0) = 0;
                elseif strcmp(pft_index, 'compositions') == 1
                    CV_modeled_pigs(CV_modeled_pigs < 0) = 0;
                    CV_modeled_pigs(CV_modeled_pigs > 1) = 1;
                elseif strcmp(pft_index, 'EOFs') == 1
                end

                percent_errors(1:length(CV_valid_pigs),l) = ((CV_valid_pigs - CV_modeled_pigs)./CV_valid_pigs).*100;
                mean_percent_error(l) = mean(abs(percent_errors(:,l)));
                median_percent_error(l) = median(abs(percent_errors(:,l)));

                % fit linear model to look at modeled vs. observed:
                lin_mdl_4_validation = fitlm(CV_valid_pigs,CV_modeled_pigs);
                R2s(l) = lin_mdl_4_validation.Rsquared.Ordinary;
                RMSEs(l) = lin_mdl_4_validation.RMSE;
            end

            if strcmp(mdl_pick_metric, 'R2') == 1
                ii = find(R2s == max(R2s));
                if length(ii) > 1
                    ii = min(ii);
                end
                n_modes_to_use(j) = ii; % for selecting based on R^2 of predicted vs. observed
            elseif strcmp(mdl_pick_metric, 'RMSE') == 1
                ii = find(RMSEs == min(RMSEs));
                if length(ii) > 1
                    ii = min(ii);
                end
                n_modes_to_use(j) = ii; % for selecting based on minimizing RMSE
            elseif strcmp(mdl_pick_metric, 'avg') == 1
                ii = find(mean_percent_error == min(mean_percent_error));
                if length(ii) > 1
                    ii = min(ii);
                end
                n_modes_to_use(j) = ii; % for selecting based on minimizing avg percent error
            elseif strcmp(mdl_pick_metric, 'med') == 1
                ii = find(median_percent_error == min(median_percent_error));
                if length(ii) > 1
                    ii = min(ii);
                end
                n_modes_to_use(j) = ii; % for selecting based on minimizing median percent error
            end

            % apply your optimized model and record the g.o.f. statistics for this k-th CV:
            the_lin_mdl = fitlm(CV_AFs_train(:,1:n_modes_to_use(j)), CV_train_pigs); % for an MLR model w/ derivative components
            these_betas = the_lin_mdl.Coefficients;
            these_betas = table2array(these_betas(:,1));
            alpha(j) = these_betas(1); %intercept
            these_betas(1) = []; % model coefficients for each spectral AF

            % now turn model coefficients for AF's into coefficients for the combined
            % derivative spectra:
            betas(:,j) = CV_EOFs_train(:,1:n_modes_to_use(j)) * these_betas;

            % 5-CV Model Validation:
            CV_modeled_pigs = (CV_valid_spec * betas(:,j)) + alpha(j);
            
            % another if statement for constraining mdl outputs by pft_index
            if strcmp(pft_index, 'pigment') == 1
                CV_modeled_pigs(CV_modeled_pigs < 0) = 0;
            elseif strcmp(pft_index, 'compositions') == 1
                CV_modeled_pigs(CV_modeled_pigs < 0) = 0;
                CV_modeled_pigs(CV_modeled_pigs > 1) = 1;
            elseif strcmp(pft_index, 'EOFs') == 1
            end
            
            % fit linear model to look at modeled vs. observed:
            lin_mdl_4_validation = fitlm(CV_modeled_pigs,CV_valid_pigs);
            CV_R2s(j) = lin_mdl_4_validation.Rsquared.Ordinary;
            CV_RMSEs(j) = lin_mdl_4_validation.RMSE;
        end

        % so now you have k sets of optimized coefficients. grab the
        % average of them and validate against your original 25% validation
        % data set. 

        % Store mean/std of each set of k-fold CV betas (the model coefficients for the ith run of the n_permutations):
        mean_betas(:,i) = mean(betas,2);
        mean_alphas(i) = mean(alpha);
        std_betas(:,i) = std(betas,0,2);
        std_alphas(i) = std(alpha);

        % Put beta's and alpha into non-standardized space so you can apply
        % them straight to the validation spectra (alternatively, you could put
        % the validation spectra in standardized space and apply it that way).
        % They work the same (I've tested) and this way's easier:
        mean_betas_nonstd(:,i) = mean_betas(:,i) ./ std(spectra_4_mdl_training)';    
        mean_alphas_nonstd(i) = mean_alphas(i) - sum(mean_betas(:,i) .* (mean(spectra_4_mdl_training)' ./ std(spectra_4_mdl_training)'));

        % Validate on the data you set aside previously for this ith run of the n_permutations using mean betas of this
        % permutation from cross-validation, and store g.o.f stats across permutations:
        modeled_pigs = (spectra_4_mdl_validate * mean_betas_nonstd(:,i)) + mean_alphas_nonstd(i);

        % constrain modeled pigments:
        if strcmp(pft_index, 'pigment') == 1
            modeled_pigs(modeled_pigs < 0) = 0;
        elseif strcmp(pft_index, 'compositions') == 1
            modeled_pigs(modeled_pigs < 0) = 0;
            modeled_pigs(modeled_pigs > 1) = 1;
        elseif strcmp(pft_index, 'EOFs') == 1
        end

        lin_mdl_4_validation = fitlm(modeled_pigs,pigs_validate);
        R2s_final(i) = lin_mdl_4_validation.Rsquared.Ordinary;
        RMSEs_final(i) = lin_mdl_4_validation.RMSE;

        pigs_validate(pigs_validate == 0) = 1e-4; % add a really small # to zeros to make % error calculations reasonable (pft_index-independent)
        % calc all the validation statistics
        pct_bias(i) = mean(((modeled_pigs - pigs_validate)./pigs_validate)*100);
        pct_errors(i,:) = abs(((modeled_pigs - pigs_validate)./pigs_validate)*100);
        med_pct_error(i) = median(pct_errors(i,:));
        avg_pct_error(i) = mean(pct_errors(i,:));
        sort_pct_errors = sort(pct_errors(i,:));
        CI_pct_error(i) = sort_pct_errors(ceil((0.95 * size(pct_errors(i,:),2))));
        std_pct_error(i) = std(pct_errors(i,:));
        
        % display permutation number here...
        disp(['hey dude, im doing good. still gonna be a while, im only on permutation number ',num2str(i)]);

    end
    
    % compile outputs and save them appropriately.
    coefficients = mean_betas_nonstd;
    intercepts = mean_alphas_nonstd;
    
    % get summary gof stats and pop them in a table:
    summary_gofs = [mean(R2s_final),std(R2s_final),mean(RMSEs_final),std(RMSEs_final),mean(avg_pct_error),std(avg_pct_error),...
        mean(med_pct_error),std(med_pct_error),mean(pct_bias),std(pct_bias)];
    summary_gofs = array2table(summary_gofs);
    summary_gofs.Properties.VariableNames = {'Mean_R2','SD_R2','Mean_RMSE','SD_RMSE',...
        'Mean_mean_pct_error','SD_mean_pct_error','Mean_median_pct_error','SD_median_pct_error',...
        'Mean_pct_bias','SD_pct_bias'};
    
    % get all gof stats from n_permutations CV's and pop them into a struct
    all_gofs = struct('R2s',R2s_final,'RMSEs',RMSEs_final,'mean_pct_error',avg_pct_error,...
        'median_pct_error',med_pct_error,'pct_bias',pct_bias, 'all_pct_errors', pct_errors);
    
    % save everything to a .mat...
    save(output_file_name,'coefficients','intercepts','summary_gofs','all_gofs');
end

