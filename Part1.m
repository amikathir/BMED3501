function mtf_ramlak_shepp_hann
% MTF_RAMLAK_SHEPP_HANN
% Creates an MTF-like curve comparing three CT reconstruction filters:
%   1) Ram-Lak
%   2) Shepp-Logan
%   3) Hann
%
% Frequency axis is in LP/mm from 0 to 1. The y-axis is plotted as a percent.

    clear; clc; close all;

    % 1) Frequency axis: 0 to 1 LP/mm
    fmax = 1.0;           % max spatial frequency (LP/mm)
    df   = 0.01;          % step size
    f    = 0:df:fmax;     % frequency array

    % 2) Define ramp function (normalized from 0 to 1)
    rampNorm = f / fmax;  

    % 3) Compute filter responses
    ramLak      = rampNorm;  
    sheppLogan  = rampNorm .* sheppLoganWindow(f, fmax);
    hannFlt     = rampNorm .* hannWindow(f, fmax);

    % 4) Convert to percentage modulation (0-100%)
    ramLak_pct      = 100 * ramLak;
    sheppLogan_pct  = 100 * sheppLogan;
    hannFlt_pct     = 100 * hannFlt;

    % 5) Plot the MTF curves
    figure('Name','Filter MTF Comparison');
    hold on; grid on; box on;

    plot(f, ramLak_pct,      '-k','LineWidth',2, 'DisplayName','Ram-Lak filter');
    plot(f, sheppLogan_pct,  '-r','LineWidth',2, 'DisplayName','Shepp-Logan filter');
    plot(f, hannFlt_pct,     '-b','LineWidth',2, 'DisplayName','Hann filter');

    xlabel('Spatial Frequency (LP/mm)');
    ylabel('Modulation (%)');
    title('Comparison of Ram-Lak, Shepp-Logan, and Hann MTFs');
    legend('Location','northeast');
    ylim([-5, 105]);  % Add a small margin above and below
    hold off;
end

%% ---------- HELPER WINDOW FUNCTIONS ---------- %%
function w = sheppLoganWindow(f, fmax)
% SHEPP-LOGAN = sinc(pi*f/fmax). Zero out if f>fmax.
    w = zeros(size(f));
    mask = (f <= fmax);
    x = f(mask) / fmax; 
    w(mask) = sin(pi*x) ./ (pi*x);
    w(abs(x) < 1e-12) = 1;  % Handle x = 0
end

function w = hannWindow(f, fmax)
% HANN = 0.5(1 + cos(pi*f/fmax)) for |f|<fmax, else 0.
    w = zeros(size(f));
    mask = (f <= fmax);
    x = f(mask) / fmax;
    w(mask) = 0.5 * (1 + cos(pi*x));
end
