Settings_and_Defaults_v4;

ENA_gains = [15 18 20 23];
gains_1Q_m6 = ENA_gains - [0.17 0.65 1.4 2];

fig = figure; 
hold on; 

plot(ENA_gains, gains_1Q_m6, '.', 'MarkerSize', 30);
plot([0, 40], [0, 40], '--', 'Color', SD.myred)

xlim([0, 25]);
ylim([0, 25]);
xlabel('ENA gain [dB]');
ylabel('1Q gain - 6 [dB]');

hold off; 