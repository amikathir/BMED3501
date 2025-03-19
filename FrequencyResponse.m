clear; clc; close all;

f = linspace(-0.5, 0.5, 1000);  % 1000 points for a smooth curve
ramp = abs(f);                 % Basic ramp (Ram-Lak)

fcut = 0.5;  % band limit
sheppLogan = ramp .* sheppLoganWindow(f, fcut);

hannFilter = ramp .* hannWindow(f, fcut);

figure('Name','Filter Frequency Responses');
plot(f, ramp,         'k','LineWidth',2, 'DisplayName','Ram-Lak');
hold on; grid on; box on;
plot(f, sheppLogan,   'r','LineWidth',2, 'DisplayName','Shepp-Logan');
plot(f, hannFilter,   'b','LineWidth',2, 'DisplayName','Hann');

xlabel('Spatial Frequency');
ylabel('Filter Amplitude');
title('Comparison of Filter Responses');
legend('Location','best');
xlim([-0.5, 0.5]);
ylim([0, 1.1*max(ramp(:))]);  

function w = sheppLoganWindow(f, fcut)
    w = zeros(size(f));
    mask = (abs(f) <= fcut);
    x = f(mask) / fcut;          
    w(mask) = sin(pi*x)./(pi*x);
    w(abs(x)<1e-10) = 1;         
end

function w = hannWindow(f, fcut)
    w = zeros(size(f));
    mask = (abs(f) <= fcut);
    x = f(mask) / fcut;
    w(mask) = 0.5*(1 + cos(pi*x));
end
