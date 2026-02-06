function [GAINS, MODES, FREQS, EXPANSION] = bSPOD(data, Nf, varargin)
%   Band-ensemble SPOD (bSPOD) with frequency attribution
%
%   [GAINS, MODES, FREQS] = bSPOD(data, Nf)
%   [GAINS, MODES, FREQS, EXPANSION] = bSPOD(data, Nf, ...)
%
%   INPUTS
%   data : (Nt x Ndof) time-series data (TIME in rows!, DOF in columns!)
%   Nf   : integer, filter width / number of consecutive Fourier bins per band / number of POD modes (Nmodes = Nf)
%
%   Name-value optional parameters (with defaults):
%   'dt'      : time step (default 1)
%   'ell'     : vector of exponents for frequency attribution (default [inf 1 0])
%               Frequencies are returned for each ell 
%               ell = inf - frequency of most contributing Fourier mode
%               ell = 1   - frequencies weighted by respectective Fourier mode contribution
%               ell = 0   - band center frequency
%   'weights' : spatial weights vector (Ndof x 1) / (1 x Ndof), applied elementwise
%              (default: ones(Ndof,1))
%   'ovlp'    : integer, number of overlapping Fourier bins between bands (default 0)
%   'window'  : 0/1, if 1 apply Hann taper to time series before FFT (default 0)
%
%             -> the number of frequnecy bands is 
%             Nwin = floor((Nt − Nf)/HOP) + 1 with the band hop size HOP = Nf − ovlp      
%
%   OUTPUTS
%   GAINS : (Nwin x Nmodes) bSPOD gains per band 
%   MODES : (Nwin x Nmodes x Ndof) bSPOD mode shapes
%   FREQS : (Nwin x Nmodes x (numel(ell))) frequencies for each exponent in ell:
%           FREQS(:,:,k) corresponds to ell(k) exponent
%   EXPANSION (optional) : (Nwin x Nmodes x Nmodes) eigenvectors V per band
%
%   For feedback and questions please contact:
%   j.vonsaldern@tu-berlin.de
%   (c) Jakob G.R. von Saldern, February 5th 2026
%
%   This file is licensed under the GNU General Public License v3.0.

% -------------------------- input parsing ---------------------------------
ip = inputParser;
ip.FunctionName = mfilename;

addRequired(ip, 'data', @(x) isnumeric(x) && ismatrix(x) && ~isempty(x));
addRequired(ip, 'Nf',   @(x) isnumeric(x) && isscalar(x) && x == floor(x) && x > 0);

addParameter(ip, 'dt',      1,           @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(ip, 'ell',     [inf 1 0],   @(x) isnumeric(x) && isvector(x) && ~isempty(x));
addParameter(ip, 'weights', [],          @(x) isnumeric(x) && ~isempty(x));
addParameter(ip, 'ovlp',    0,           @(x) isnumeric(x) && isscalar(x) && x == floor(x) && x >= 0);
addParameter(ip, 'window',  0,           @(x) isnumeric(x) && isscalar(x) && any(x == [0 1]));

parse(ip, data, Nf, varargin{:});

dt      = ip.Results.dt;
ell     = ip.Results.ell(:).';      % row vector
ovlp    = ip.Results.ovlp;
useWin  = ip.Results.window;

% -------------------------- dimensions & checks ----------------------------
[Nt, Ndof] = size(data);
Nmodes = Nf;

assert(ovlp >= 0 && ovlp < Nf, 'bSPOD:Overlap', 'ovlp must satisfy 0 <= ovlp < Nf.');
assert(Nt >= Nf, 'bSPOD:DataLength', 'Nt must be >= Nf.');

% band-ensemble shift
HOP = Nf - ovlp; 

W = ip.Results.weights;
if isempty(W)
    W = ones(Ndof,1);
end
assert(isvector(W), 'bSPOD:WeightsType', 'weights must be a vector of length Ndof.');
W = W(:);
assert(numel(W) == Ndof, 'bSPOD:WeightsSize', 'weights must have length Ndof.');

% -------------------------- FFT (time -> frequency) ------------------------
if useWin == 1
    win = hann(Nt);                                  % modify for other taper function
    Uw  = sum(win.^2) / Nt;
    Q_hat = fft(bsxfun(@times, data, win), [], 1);   % (Nt x Ndof)
else
    Uw  = 1;
    Q_hat = fft(data, [], 1);                        % (Nt x Ndof)
end
Q_hat = (Q_hat.' / Nt);                              % (Ndof x Nt) forward normalization

% -------------------------- frequency vector -------------------------------
Fs = 1/dt;                                        % sampling freq 
df = Fs/Nt;                                       % freq. resolution dft
frequ = (0:Nt-1) * df;                            % dft frequencies
frequ(frequ > Fs/2) = frequ(frequ > Fs/2) - Fs;   % shift to negative frequencies

% -------------------------- banding in frequency ---------------------------
Nwin = floor((Nt - Nf)/HOP) + 1;

% -------------------------- preallocation ----------------------------------
GAINS = zeros(Nwin, Nmodes);
MODES = zeros(Nwin, Nmodes, Ndof);
FREQS = zeros(Nwin, Nmodes, numel(ell) + 1);

if nargout >= 4
    EXPANSION = zeros(Nwin, Nmodes, Nmodes);
else
    EXPANSION = [];
end

% constant factor (matches original scaling intent)
alpha = (Nt*dt) / (Nf*Uw);

% -------------------------- main loop --------------------------------------
for j = 1:Nwin
    idx0 = (j-1)*HOP + 1;
    idx1 = idx0 + Nf - 1;

    f_window = frequ(idx0:idx1);              % (1 x Nf) or (Nf x 1)
    Q_block  = Q_hat(:, idx0:idx1);           % (Ndof x Nf)

    % cross-spectral density in the band (frequency-focused)
    % W is (Ndof x 1)
    WQ = bsxfun(@times, Q_block, W);   % (Ndof x Nf)
    C  = alpha * (Q_block' * WQ);      % (Nf x Nf)

    % eigendecomposition
    [V, D] = eig(C, 'vector');
    [lambda, p] = sort(D, 'descend');   % gains are real for Hermitian C
    V = V(:, p);
    lambda = lambda(:);               

    % bSPOD modes
    denom = sqrt(lambda);
    U = sqrt(alpha) * (Q_block * V) * diag(1 ./ denom);  % (Ndof x Nf)
    U = U.';                                             % (Nf x Ndof)

    MODES(j,:,:) = U;
    GAINS(j,:)   = lambda(:);

    % keep original "real-signal one-sided doubling" behavior:
    if isreal(data)
        if j ~= 1 && j ~= Nwin
            GAINS(j,:) = 2 * GAINS(j,:);
        end
    end

    % Optional output: expansion eigenvectors
    if nargout >= 4
        EXPANSION(j,:,:) = V; % raw, unnormalized expansion coeffs.
    end

    % ---------- frequency attribution using |V|.^ell ----------
    Wabs = abs(V);                 
    
    % finite ell
    for k = 1:numel(ell)
        if ell(k)==inf
            % ell = inf: max of |V| 
            % (frequency of most contributing Fourier mode)
            [~, imax] = max(Wabs, [], 1);
            FREQS(j,:,k) = f_window(imax);
        else
            wk = Wabs.^ell(k);
            sk = sum(wk, 1);
            wk = bsxfun(@rdivide, wk, sk);
            FREQS(j,:,k) = (f_window(:).' * wk).';
        end
    end

    disp(sprintf('Frequency band %d out of %d', j, Nwin));
end

end
