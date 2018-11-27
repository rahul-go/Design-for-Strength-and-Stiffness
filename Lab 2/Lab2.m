%% Reset
clear all;
close all;
clc;



%% Solid Element
stress = [7.889, 8.583, 8.645, 8.664];
deflection = [3.735, 3.924, 3.940, 3.948] * 10^-1;
global_size = [8, 3, 2, 1];

figure;
plot(global_size, stress, 'LineWidth', 2);
title('Stress vs. Global Size (Solid Element)');
xlabel('Global Size (in)');
ylabel('Stress (kpsi)');
figure;
plot(global_size, deflection, 'LineWidth', 2);
title('Deflection vs. Global Size (Solid Element)');
xlabel('Global Size (in)');
ylabel('Deflection (in)');



%% Shell Element
stress = [8.542, 8.647, 8.667, 8.676];
deflection = [3.839, 3.927, 3.937, 3.939] * 10^-1;
global_size = [8, 4, 2, 1];

figure;
plot(global_size, stress, 'LineWidth', 2);
title('Stress vs. Global Size (Shell Element)');
xlabel('Global Size (in)');
ylabel('Stress (kpsi)');
figure;
plot(global_size, deflection, 'LineWidth', 2);
title('Deflection vs. Global Size (Shell Element)');
xlabel('Global Size (in)');
ylabel('Deflection (in)');



%% Beam Element
stress = [8.640, 8.640, 8.640, 8.640];
deflection = [3.921, 3.921, 3.921, 3.921] * 10^-1;
global_size = [8, 4, 2, 1];

figure;
plot(global_size, stress, 'LineWidth', 2);
title('Stress vs. Global Size (Beam Element)');
xlabel('Global Size (in)');
ylabel('Stress (kpsi)');
figure;
plot(global_size, deflection, 'LineWidth', 2);
title('Deflection vs. Global Size (Beam Element)');
xlabel('Global Size (in)');
ylabel('Deflection (in)');