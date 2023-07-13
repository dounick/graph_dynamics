function sampling = create_Sampling_GraphCompletion(n,p,mode,...
    Hermitian_flag,parameters_sampling)
%create_Sampling_GraphCompletion This function creates the set of indices
%Omega containing the indices of [T,T^2,T^3,...,T^p] that are sampled.
%    Parameters
%    ----------
%       mode: sampling pattern to create for the given problem.
%            == 'row_steps_column_random': Entrywise sampling of T.
%               Fixed rows given by 'step_rows', nr. of 
%               'percentage_columns' random columns
%               to be sampled (both fixed for all powers of T)
%            == 'random_fixed_sensors': Entrywise sampling of T.
%               Sample uniform at random 'm_single' entries (or diagonal entries/
%               entry pairs if Hermitian_flag == 1), fixed for all powers
%               of T.
%            == 'random_indep_sensors': Entrywise sampling of T.
%               Sample uniform at random m_single entries (or diagonal entries/
%               entry pairs if Hermitian_flag == 1), independently at random
%               for all powers of T.
%            == 'random_adaptive_probs': Entrywise sampling of T.
%               Sample an expected total number of 'm_total_expected' entries
%               for different powers of T, distributed
%               across different powers of T based on a density prescribed
%               by 'coherences_dist'.
%            == 'random_fixed_sensors_DMD': Entrywise sampling of T^k*x_first
%               for k = 1,...,p (dynamic mode decomposition, DMD, setting).
%               Sample uniform at random 'm_single' entries of  
%               T^k*x_first, fixed sampling set for all  k = 1,...,p. 
%            == 'random_indep_sensors_DMD': Entrywise sampling of T^k*x_first
%               for k = 1,...,p (dynamic mode decomposition, DMD, setting).
%               Sample uniform at random 'm_single' entries of  
%               T^k*x_first, random independent sampling sets for all  k = 1,...,p.               
if strcmp(mode,'row_steps_column_random') || strcmp(mode,'random_fixed_sensors') ...
        || strcmp(mode,'random_indep_sensors') || strcmp(mode,'random_adaptive_probs')
   DMD_flag = 0;
elseif strcmp(mode,'random_fixed_sensors_DMD') || strcmp(mode,'random_indep_sensors_DMD')
   DMD_flag = 1;
else
    error('Mode is not well-specified.')
end

sampling = struct;
if size(p) == 1
    sampling.powers = [1:p];
else
    sampling.powers = p;
end
sampling.Omega_cell = cell(1,length(sampling.powers));

if strcmp(mode,'random_fixed_sensors') || strcmp(mode,'full_DMD') || strcmp(mode,'random_fixed_sensors_DMD') ||...
    strcmp(mode,'random_indep_sensors')
    if isfield(parameters_sampling,'m_single')
        m_single = parameters_sampling.m_single;
    else
        error("A sampling scheme was chosen that requires the parameter 'm_single' for the number of sensors per time step / matrix power.");
    end
end
if strcmp(mode,'random_adaptive_probs')
    if isfield(parameters_sampling,'m_total_expected')
        m_total_expected = parameters_sampling.m_total_expected;
    else
        error("Since the sampling scheme is 'random_adaptive_probs', please provide the total number of provided samples 'm_total_expected'.")
    end
    if isfield(parameters_sampling,'coherences_dist')
        coherences_dist = parameters_sampling.coherences_dist;
    else
        error("Since the sampling scheme is 'random_adaptive_probs', please provide a set of parameters 'coherences_dist' corresponding to the desired density of the sampling probability.")
    end
else
    
end

if strcmp(mode,'row_steps_column_random')
    if isfield(parameters_sampling,'step_rows')
        step_rows = parameters_sampling.step_rows;
    else
        error("Since the sampling scheme is 'row_steps_column_random', please set the parameter 'step_rows' in 'parameters_sampling' as an integer between 1 and n.");
    end
    if isfield(parameters_sampling,'percentage_columns')
        percentage_columns = parameters_sampling.percentage_columns;
    else
        error("Since the sampling scheme is 'row_steps_column_random', please set the parameter 'percentage_columns' in 'parameters_sampling' as a number between 0 and 1.")
    end
    col_inds = randperm(n,floor(percentage_columns*n));
    row_inds = 1:step_rows:n;
    row_col_kron = kron(ones(length(row_inds),1),(col_inds-1)*(n)+1)+kron(row_inds'-1,ones(1,length(col_inds))); % Mimic indexing: T(J,I);
    if Hermitian_flag
        Omega_single = get_symmetric_indices(sort(row_col_kron(:)),n);
        m_single = length(Omega_single);
    else
        col_row_kron = kron(ones(length(col_inds),1),(row_inds-1)*(n)+1)+kron(col_inds'-1,ones(1,length(row_inds)));
        Om1 = sort(row_col_kron(:));
        Om2 = sort(col_row_kron(:));
        tmp_Om = unique(reshape([Om1,Om2],1,[]),'stable');
        Omega_single = sort(tmp_Om(:));
        m_single = length(Omega_single);
    end
    for k=1:length(sampling.powers)
        sampling.Omega_cell{k} = Omega_single;
    end
elseif strcmp(mode,'random_fixed_sensors')
    if Hermitian_flag
        Omega_single = sort(randperm(n*(n+1)/2,m_single)');%get_symmetric_indices(sort(randperm(n*(n+1)/2,m_single)),n);  
    else
        Om1 = sort(randperm(n^2,m_single));
        [row_ind,col_ind] = ind2sub([n,n],Om1);
        Om2 = sort(sub2ind([n,n],col_ind,row_ind));
        tmp_Om = unique(reshape([Om1,Om2],1,[]),'stable');
        Omega_single = sort(tmp_Om(:));
    end
    for k=1:length(sampling.powers)
        sampling.Omega_cell{k} = Omega_single;
    end
    m_single = length(Omega_single);
elseif strcmp(mode,'full_DMD')
    if not(m_single == n)
       error('Full DMD sampling requires a number of n samples.') 
    end
    
elseif strcmp(mode,'random_fixed_sensors_DMD')
    if Hermitian_flag
        error('DMD only works right now for non-Hermitian matrix space.')
    else
        Omega_single = sort(randperm(n,m_single));
        for k=1:length(sampling.powers)
            sampling.Omega_cell{k} = Omega_single;
        end
    end
elseif strcmp(mode,'random_indep_sensors')
    sampling.Omega_cell = cell(1,length(sampling.powers));
    for k=1:length(sampling.powers)
        if Hermitian_flag
            sampling.Omega_cell{k} = sort(randperm(n*(n+1)/2,m_single)');%get_symmetric_indices(sort(randperm(n^2,m_single)),n);  
        else
            Om1 = sort(randperm(n^2,m_single));
            [row_ind,col_ind] = ind2sub([n,n],Om1);
            Om2 = sort(sub2ind([n,n],col_ind,row_ind));
            tmp_Om = unique(reshape([Om1,Om2],1,[]),'stable');
            sampling.Omega_cell{k} = sort(tmp_Om(:));
        end
    end
elseif strcmp(mode,'random_indep_sensors_DMD')
    if Hermitian_flag
        error('DMD only works right now for non-Hermitian matrix space.')
    else
        for k=1:length(sampling.powers)
            tmp_Om = sort(randperm(n,m_single));
            sampling.Omega_cell{k} = tmp_Om;
        end
    end
elseif strcmp(mode,'random_adaptive_probs')
    sampling.Omega_cell = cell(1,length(sampling.powers));
    coherences_values=coherences_dist(sampling.powers);
    sum_coh=sum(coherences_values);
    if Hermitian_flag
        normalization_fac = 2*m_total_expected/(n*(n+1)*sum_coh);
    else
        normalization_fac = m_total_expected/(n^2*sum_coh);
    end
    sampling.normalization_fac = normalization_fac;
    m_total = 0;
    for k=1:length(sampling.powers)
        if Hermitian_flag
            sampling_c = binornd(1,normalization_fac*coherences_dist(sampling.powers(k)),1,n*(n+1)/2);
            sampling.Omega_cell{k} = find(sampling_c);
        else
            sampling_c = binornd(1,normalization_fac*coherences_dist(sampling.powers(k)),1,n^2);
            Om1 = find(sampling_c);
            [row_ind,col_ind] = ind2sub([n,n],Om1);
            Om2 = sort(sub2ind([n,n],col_ind,row_ind));
            tmp_Om = unique(reshape([Om1,Om2],1,[]),'stable');
            sampling.Omega_cell{k} = sort(tmp_Om(:));
            normalization_fac = m_total_expected/(n^2*sum_coh);
        end
        m_total = m_total + length(sampling.Omega_cell{k});
    end
end

if strcmp(mode,'random_adaptive_probs')
    Omega = zeros(m_total,1);
else
    Omega = zeros(m_single*length(sampling.powers),1);
end

m_cumulative = 0;
if DMD_flag
    if Hermitian_flag
        
    else
        for k=1:length(sampling.powers)
            Omega(m_cumulative+1:m_cumulative+length(sampling.Omega_cell{k})) = ...
                sampling.Omega_cell{k}+(sampling.powers(k)-1)*n;
            m_cumulative = m_cumulative+length(sampling.Omega_cell{k});
        end
        Omega_complement = setdiff([1:n*sampling.powers(end)]',Omega);
    end
else
    if Hermitian_flag
        for k=1:length(sampling.powers)
            Omega(m_cumulative+1:m_cumulative+length(sampling.Omega_cell{k})) = ...
                sampling.Omega_cell{k}+(sampling.powers(k)-1)*((n+1)*n/2);
            m_cumulative = m_cumulative+length(sampling.Omega_cell{k});
            %Omega((k-1)*m_single+1:k*m_single) = sampling.Omega_cell{k}+(sampling.powers(k)-1)*((n+1)*n/2);
            %sampling.Omega_cell{k} = Omega_single;
        end
        Omega_complement = setdiff([1:(n*(n+1)/2)*sampling.powers(end)]',Omega);
    else
        for k=1:length(sampling.powers)
            Omega(m_cumulative+1:m_cumulative+length(sampling.Omega_cell{k})) = ...
                sampling.Omega_cell{k}+(sampling.powers(k)-1)*n^2;
            m_cumulative = m_cumulative+length(sampling.Omega_cell{k});
            %Omega((k-1)*m_single+1:k*m_single)=sampling.Omega_cell{k}+(sampling.powers(k)-1)*n^2;
            %sampling.Omega_cell{k} = Omega_single;
        end
        Omega_complement = setdiff([1:n^2*sampling.powers(end)]',Omega);
    end
end
sampling.Omega = Omega;
sampling.Omega_complement = Omega_complement;
sampling.Hermitian_flag = Hermitian_flag;
sampling.n = n;
end

