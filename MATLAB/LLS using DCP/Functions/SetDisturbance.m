function [SensorNoise, ForceDist] = SetDisturbance(t_simu,Ts,Tu)
% Disturbance and noise profile for spiral writing
%% disturbance parameters
% ============================================================
%   Copyright (c) 2008-, Xu Chen, xuchen@cal.berkeley.edu
%   Author(s): Xu Chen 
% ============================================================



% Ts = PlantData.Ts;
% Tu = PlantData.Tu;

% Data from plant model
Fns = 1/Ts/2; % Nyquist freq. of PES sampling
Fnu = 1/Tu/2; % Nyquist freq. of plant input

% dFs = 6; % Frequency resolution of PES
% dFu = 6; % Frequency resolution of plant input
%% modified for any L
%% ======================== modified =====================================
div_s = divisors(Fns);
div_u = divisors(Fnu);
n_s = length(div_s);
n_u = length(div_u);
i_s = 1;
i_u = 1;
while i_s < n_s
    if div_s(i_s) <= round(Fns/100)
    i_s = i_s+1;
    else
        break
    end
end
dFs = div_s(i_s);

while i_u < n_u
    if div_s(i_u) <= round(Fnu/(Ts/Tu*100))
    i_u = i_u+1;
    else
        break
    end
end
dFu = div_u(i_u);
%% ======================== modified =====================================

% Check frequency resolution
if ( (Fns/dFs) ~= fix(Fns/dFs) )
    error('Frequency resolution dF0s is not appropriate')
end
if ( (Fnu/dFu) ~= fix(Fnu/dFu) )
    error('Frequency resolution dF0u is not appropriate')
end

Nums  = Fns/dFs;
Numu  = Fnu/dFu;
Freqs = dFs*(1:Nums)';
Frequ = dFu*(1:Numu)';

% Scaling parameter
scl = sqrt(2*Ts);

% Sensor noise
AmpSensorNoise = 1.5e-2*scl; 
%       unit:[Track], 1 sigma of Sensor Noise

% Force disturbance
AmpForceDist = 1.0e-4*scl; 
%       unit:[A] at control input, 1 sigma of Torque Noise with sampling Ts

% Seed for random signal, the same seed results in the same output from the
% rand function. if we don't set the seed, all dist will be the same.
Seed_ForceDist    = 2; 
Seed_SensorNoise  = 3;

% Frequency vector
DistParam.Nums              = Nums;
DistParam.Numu              = Numu;
DistParam.dFs               = dFs;
DistParam.dFu               = dFu;
DistParam.Freqs             = Freqs;
DistParam.Frequ             = Frequ;
% Sensor noise
DistParam.AmpSensorNoise    = AmpSensorNoise;
% Force disturbance
DistParam.AmpForceDist      = AmpForceDist;
% Seeds for random signal
DistParam.Seed_ForceDist    = Seed_ForceDist;
DistParam.Seed_SensorNoise  = Seed_SensorNoise;
DistParam.Tu                = Tu;
DistParam.Ts                = Ts;
DistParam.t_simu            = t_simu;

%% Data output
ForceDist = SetForceDist(DistParam);
SensorNoise = SetSensorNoise(DistParam);

%% set sensor noise
    function [SensorNoise] = SetSensorNoise(DistParam);
        %SetSensorNoise  Set sensor noise in time and frequency domain

        % simulation parameters
        Ts                 = DistParam.Ts; % PES Sampling
        % seeds for random signal
        Seed_SensorNoise  = DistParam.Seed_SensorNoise;

        % Frequency vector is defined by PES sampling frequency
        Freq = DistParam.Freqs;
        NUM  = DistParam.Nums;

        % Sensor noise in time domain
        NUM_SensorNoise = DistParam.t_simu/Ts+1;
        randn('state',Seed_SensorNoise);
        Data = whitenoise(NUM_SensorNoise,Ts)*DistParam.AmpSensorNoise;
        Time = (0:NUM_SensorNoise-1)'*Ts;

        % Sensor noise in frequency domain
        Spec = DistParam.AmpSensorNoise*ones(NUM,1);

        % Output parameters in time domain
        SensorNoise.Data = Data;
        SensorNoise.Time = Time;
        SensorNoise.Ts   = Ts;
        % Output parameters in frequency domain
        SensorNoise.Freq = Freq;
        SensorNoise.Spec = Spec;
    end
%% EOF of SetSensorNoise.m

%% set force dist
    function [ForceDist] = SetForceDist(DistParam);
        % Read parameters
        AmpForceDist       = DistParam.AmpForceDist;
        Ts                 = DistParam.Ts; % PES Sampling
        Tu                 = DistParam.Tu; % Plant input
        % Seeds for random signal
        Seed_ForceDist  = DistParam.Seed_ForceDist;
        % Frequency vector is defined by PES sampling frequency
        Freq = DistParam.Freqs;
        NUM  = DistParam.Nums;

        % Force dist. in time domain
        if fix(Ts/Tu) ~= (Ts/Tu)
            error('Ts/Tu must be integer');
        end
        NUM_ForceDist = DistParam.t_simu/Tu+1;
        randn('state',Seed_ForceDist);
        Data       = whitenoise(NUM_ForceDist,Tu)*DistParam.AmpForceDist;
        Time       = (0:NUM_ForceDist-1)'*Tu;
        % dist at PES
%         DistOut    = lsim(PlantData.Pn,Data,Time);
%         dcidx      = 1:(Ts/Tu):NUM_ForceDist;
%         DataAtPes  = DistOut(dcidx); % down sampling by (Ts/Tu)
%         TimeAtPes  = (0:length(DataAtPes)-1)'*Ts;

        % Force dist. in frequency domain
        Spec      = DistParam.AmpForceDist*ones(NUM,1);
%         [mag,phs] = bode( PlantData.Pn*DistParam.AmpForceDist, 2*pi*Freq );
%         SpecAtPes = mag(:);

        % Output parameters in time domain
        ForceDist.Data      = Data;
        ForceDist.Time      = Time;
        ForceDist.Ts        = Tu;
%         ForceDist.DataAtPes = DataAtPes;
%         ForceDist.TimeAtPes = TimeAtPes;
%         ForceDist.TsAtPes   = Ts;
        % Output parameters in frequency domain
        ForceDist.Freq      = Freq;
        ForceDist.Spec      = Spec;
%         ForceDist.SpecAtPes = SpecAtPes;
    end
%% EOF SetForceDist.m

end
%% EOF SetDisturbance.m
