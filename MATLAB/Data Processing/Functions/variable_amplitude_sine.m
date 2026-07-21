function [y, A, y_components] = variable_amplitude_sine(t, f, A_func, phase)
%VARIABLE_AMPLITUDE_SINE Generate multiple sinusoids with varying amplitudes.
%
% Inputs:
%   t      - Time vector, 1 x Nt or Nt x 1
%   f      - Frequency vector [f1, f2, ..., fN] in Hz
%
%   A_func - Amplitude specification:
%              1) Scalar:
%                   Same constant amplitude for all frequencies
%
%              2) Nf-element vector:
%                   Constant amplitude for each frequency
%
%              3) Nf-by-Nt matrix:
%                   Time-varying amplitude for each frequency
%
%              4) Function handle:
%                   Must return an Nf-by-Nt amplitude matrix
%
%              5) Cell array of function handles:
%                   One amplitude function per frequency
%
%   phase  - Phase offset in radians:
%              Scalar or Nf-element vector
%              Optional; default is zero
%
% Outputs:
%   y            - Sum of all sinusoidal components, same shape as t
%   A            - Amplitude envelope, Nf x Nt
%   y_components - Individual sinusoidal components, Nf x Nt

    % Save original orientation
    output_is_column = iscolumn(t);

    % Use row vectors internally
    t = t(:).';
    f = f(:);

    Nf = length(f);
    Nt = length(t);

    if nargin < 4 || isempty(phase)
        phase = zeros(Nf,1);
    elseif isscalar(phase)
        phase = phase * ones(Nf,1);
    else
        phase = phase(:);
    end

    if length(phase) ~= Nf
        error('phase must be scalar or have one element per frequency.');
    end

    % Construct amplitude matrix A: Nf x Nt
    if isa(A_func, 'function_handle')

        A = A_func(t);

    elseif iscell(A_func)

        if length(A_func) ~= Nf
            error('A_func must contain one function handle per frequency.');
        end

        A = zeros(Nf,Nt);

        for k = 1:Nf
            A(k,:) = A_func{k}(t);
        end

    elseif isscalar(A_func)

        A = A_func * ones(Nf,Nt);

    elseif isvector(A_func) && length(A_func) == Nf

        A = A_func(:) .* ones(1,Nt);

    elseif isequal(size(A_func), [Nf,Nt])

        A = A_func;

    else
        error(['A_func must be a scalar, an Nf-element vector, ', ...
               'an Nf-by-Nt matrix, a function handle, or a cell array.']);
    end

    % Allow a single returned envelope to be applied to every frequency
    if isvector(A) && length(A) == Nt
        A = repmat(A(:).',Nf,1);
    end

    if ~isequal(size(A), [Nf,Nt])
        error('The amplitude specification must produce an Nf-by-Nt matrix.');
    end

    % Generate each sinusoidal component
    angles = 2*pi*f*t + phase;
    y_components = A .* sin(angles);

    % Sum all frequencies
    y = sum(y_components,1);

    % Match the original orientation of t
    if output_is_column
        y = y.';
        A = A.';
        y_components = y_components.';
    end
end