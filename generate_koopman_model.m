% This provides interface to open-source code for Koopman data-driven 
% model-based control. It can be used to build a lifted linear model

% Setup: Add this file to the base folder of following cloned github repo:
% https://github.com/ramvasudevan/soft-robot-koopman
% Save training and validation data in baselines/[robot_name]_koopman.py


% Goes through full process of getting a linear model from data
% 1. Gather training and validation data and save in req format
% 2. Load saved diamond_koopman and train this data
% 3. Validate data
% 4. Save data to as .mat file with functionality for interface with python

%% parameters to consider for tuning
Ts = 0.1; % sampling time to consider
obs_degree = 2;  % Degree of monioials to consider
delay = 1;  % Nbr of delays to consider in observables
lasso = [10]; % Lasso parameter value
includeConst = true;
truncate_model = false; fractionModes = 1.0;
hardware = true;

%% gather training data (need to prepare data file before running this)

% load in data file(s)
[ datafile_name , datafile_path ] = uigetfile( 'datafiles/koopman_paper_import/*.mat' ,...
    'Choose training data file...' );

[ valfile_name, valfile_path ] = uigetfile('datafiles/koopman_paper_import/*.mat' ,...
    'Choose validation data file...' );

training_data = load([datafile_path, datafile_name]);
% training_data.u = training_data.u(1:end-1,:);
if hardware
    training_data.t = training_data.t';
    training_data.u = training_data.u';
    training_data.y = training_data.y';
end
training_data = data.resample(training_data, Ts);

val_data = load([valfile_path, valfile_name]);
if hardware
    val_data.t = val_data.t';
    val_data.u = val_data.u';
    val_data.y = val_data.y';
end
val_data = data.resample(val_data, Ts);

data_inst = data();

data_inst.get_data4sysid(training_data, val_data,...
    'True', 'ee_pos_20Hz');


%%
[ datafile_name , datafile_path ] = uigetfile( 'datafiles/*.mat' , 'Choose data file for sysid...' );
data4sysid = load( [datafile_path , datafile_name] );

% data4sysid = load( [datafile_path , datafile_name] );


%% construct sysid class
ksysid_inst = ksysid( data4sysid, ...
                'model_type' , 'linear' ,...    % model type (linear or nonlinear)
                'obs_type' , { 'poly' } ,...    % type of basis functions
                'obs_degree' , [ obs_degree ] ,...       % "degree" of basis functions
                'snapshots' , Inf ,...          % Number of snapshot pairs
                'lasso' , lasso ,...           % L1 regularization term
                'delays' , delay, ...          % Numer of state/input delays
                'includeConst', includeConst);         % Include constant term (modify this for EDMD)


%% train model(s)
models = ksysid_inst.train_models;
 
%% Truncate Model
A = models.model.A;
A_H = ctranspose(A);

% Obtain the eigenvectors
[v, d] = eig(A); [w, ~] = eig(A_H);

% Scale eigenvectors w_i such that (v_i, w_j) = \delta_ij
w_norm = zeros(size(w));
for i=1:size(A)
    for j=1:size(A_H)
        dot_prod = dot(v(:,i), w(:,j));
        if norm(dot_prod) >= 1e-10
            w_norm(:, i) = w(:,j)/dot_prod;
        end
    end
end

% Organize training data

[x,~,u] = deal( models.snapshotPairs.alpha , models.snapshotPairs.beta , models.snapshotPairs.u );

% Build matrices
[~, m] = deal( models.params.n , models.params.m );
Nm = models.params.N + m;   % dimension of z plus input
N = models.params.N;       % dimension of z (i.e. lifted state)

Px = zeros(length(x), models.params.N);

for i = 1:length(x)
    psix = models.lift.full( [ x(i,:) ]' )';
    Px(i,:) = psix;
end

% Calculate the power of each mode i.e. |\psi(X)| and then sort in
% decreasing order
modePower = sum(abs(Px * w_norm),1);
[modePower_sorted, idx_sorted] = sort(modePower, 'descend');

% Maintain only a fraction of modes
% Operation here assumes that the eigenfunction mode power for conjugate
% transpose is equal to the original imaginary eigenvector and hence next
% to it in the sorted matrices
n = round(fractionModes * N);
kept_idxs = idx_sorted(1:n);
v = v(:, idx_sorted);

% Transform v into jordan form basis
% TODO: There's a bug here somewhere - T is being truncated somehow
T = zeros(size(v)); % Pre-allocate change of basis
ii = 1; jj = 1;
while ii <= length(v)
    if ~isreal(v(:, ii))
        T(:, jj) = real(v(:, ii));
        T(:, jj + 1) = imag(v(:, ii));
        ii = ii + 1;
        jj = jj + 1;
    else
        T(:, jj) = v(:, ii);
    end
    ii = ii + 1;
    jj = jj + 1;
end

% Check that we are including the complex conjugate
if isequal(v(:, n), conj(v(:, n-1)))
    V = T(:, 1:n);
else
    if n < length(T)
        n = n + 1;
    end
    V = T(:, 1:n);
end


T_inv = inv(T);
V_inv = T_inv(1:n, :);

% Store projection operators for use later
if truncate_model
    models.basis.V = V;
    models.basis.W = V_inv;
else    
    models.basis.V = eye(N);
    models.basis.W = eye(N);
end

A_r = models.basis.W * A * models.basis.V;
B_r = models.basis.W * models.model.B;
C_r = models.model.C * models.basis.V;

%% validate model(s)
% could also manually do this for one model at a time
results = cell( size(models.candidates) );    % store results in a cell array
err = cell( size(models.candidates) );    % store error in a cell array

if iscell(models.candidates)
    for i = 1 : length(models.candidates)
        [ results{i} , err{i} ] = models.valNplot_model( i );
    end
else
    [ results{1} , err{1} ] = models.valNplot_model;
end

%% If not saved whilst validating, you can opt to save it separately here
models.save_class( )
    
%% Save the aforementioned model (or select one if none is defined)
[ model_name , model_path] = uigetfile( 'systems/fromData/*.mat' ,...
    'Choose model to save to python...' );
model = load([model_path, model_name]);


% Store changed matrices and parameters
model.sysid_class.model.A = A_r + 10^-6 * eye(length(A_r));
model.sysid_class.model.B = B_r;
model.sysid_class.model.C = C_r;
model.sysid_class.params.N = length(A_r);

export_to_python(model);

%% save 
% You do this based on the validation results.
% Call this function:



function export_to_python(model)
    % model should be of ksysid type
    py_data = struct();
    lin_model = struct();
    lin_model.A = model.sysid_class.model.A;
    lin_model.B = model.sysid_class.model.B;
    lin_model.C = model.sysid_class.model.C;
    lin_model.M = model.sysid_class.model.M;
    lin_model.K = model.sysid_class.model.K;
    lin_model.V = model.sysid_class.basis.V;
    lin_model.W = model.sysid_class.basis.W;
    py_data.model = lin_model;
    
    params = struct();
    params.n = model.sysid_class.params.n;
    params.m = model.sysid_class.params.m;
    params.N = model.sysid_class.params.N;
    params.nzeta = model.sysid_class.params.nzeta;
    params.Ts = model.sysid_class.params.Ts;
    params.scale = model.sysid_class.params.scale;
    params.delays = model.sysid_class.delays;
    params.obs_type = model.sysid_class.obs_type;
    params.obs_degree = model.sysid_class.obs_degree;
    
    py_data.params = params;
    save('py_model.mat', 'py_data', '-v7');

end