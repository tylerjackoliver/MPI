% Generate results from time_per_iteration and plot the stronk scaling

% Plot runtimes198

cores = [1 2 4 8 16 48];
cores = cores';
[io, time, per, gather] = time_per_iteration(runtimes198);

plot(cores, io, '-kx');
xlabel('Number of cores []', 'FontSize', 14);
ylabel('I/O time [s]', 'FontSize', 14);

saveas(gca, 'cores_io_198.eps', 'eps');

plot(cores, time(1)./time, '-bx'); 
hold on;
plot(cores, speedup, '-bx');
hold off;
xlabel('Number of cores []', 'FontSize', 14);
ylabel('Speedup []');
legend('Empirical scaling', 'Ideal scaling');
saveas(gca, 'speedup_198.eps', 'eps');

plot(cores, time, '-kx');
xlabel('Number of cores []', 'FontSize', 14);
ylabel("Iteration time [s]", 'FontSize', 14);
saveas(gca, 'cores_itertime_198.eps', 'eps');

plot(cores, per, '-kx');
xlabel('Number of cores []', 'FontSize', 14);
ylabel("Time per iteration [s]", 'FontSize', 14);
saveas(gca, 'cores_timeper_198.eps', 'eps');

plot(cores, gather, '-kx');
xlabel('Number of cores []', 'FontSize', 14);
ylabel('Gathering time [s]', 'FontSize', 14);
saveas(gca, 'cores_gather_198.eps', 'eps');

% Plot runtimes256

[io, time, per, gather] = time_per_iteration(runtimes256);

plot(cores, io, '-kx');
xlabel('Number of cores []', 'FontSize', 14);
ylabel('I/O time [s]', 'FontSize', 14);
saveas(gca, 'cores_io_256.eps', 'eps');

plot(cores, time(1)./time, '-kx', cores, speedup, '-rx');
xlabel('Number of cores []', 'FontSize', 14);
ylabel('Speedup []');
legend('Empirical scaling', 'Ideal scaling');
saveas(gca, 'speedup_256.png', 'png');

plot(cores, time, '-kx');
xlabel('Number of cores []', 'FontSize', 14);
ylabel("Iteration time [s]", 'FontSize', 14);
saveas(gca, 'cores_itertime_256.eps', 'eps');

plot(cores, per, '-kx');
xlabel('Number of cores []', 'FontSize', 14);
ylabel("Time per iteration [s]", 'FontSize', 14);

saveas(gca, 'cores_timeper_256.eps', 'eps');

plot(cores, gather, '-kx');
xlabel('Number of cores []', 'FontSize', 14);
ylabel('Gathering time [s]', 'FontSize', 14);
saveas(gca, 'cores_gather_256.eps', 'eps');

% Plot runtimes384

[io, time, per, gather] = time_per_iteration(runtimes384);

plot(cores, io, '-kx');
xlabel('Number of cores []', 'FontSize', 14);
ylabel('I/O time [s]', 'FontSize', 14);
saveas(gca, 'cores_io_384.eps', 'eps');

plot(cores, time, '-kx');
xlabel('Number of cores []', 'FontSize', 14);
ylabel("Iteration time [s]", 'FontSize', 14);
saveas(gca, 'cores_itertime_384.eps', 'eps');

plot(cores, time(1)./time, '-kx', cores, speedup, '-rx');
xlabel('Number of cores []', 'FontSize', 14);
ylabel('Speedup []');
legend('Empirical scaling', 'Ideal scaling');
saveas(gca, 'speedup_384.png', 'png');

plot(cores, per, '-kx');
xlabel('Number of cores []', 'FontSize', 14);
ylabel("Time per iteration [s]", 'FontSize', 14);
saveas(gca, 'cores_timeper_384.eps', 'eps');

plot(cores, gather, '-kx');
xlabel('Number of cores []', 'FontSize', 14);
ylabel('Gathering time [s]', 'FontSize', 14);
saveas(gca, 'cores_gather_384.eps', 'eps');

% Plot runtimes768

[io, time, per, gather] = time_per_iteration(runtimes768);

plot(cores, io, '-kx');
xlabel('Number of cores []', 'FontSize', 14);
ylabel('I/O time [s]', 'FontSize', 14);
saveas(gca, 'cores_io_768.eps', 'eps');

plot(cores, time, '-kx');
xlabel('Number of cores []', 'FontSize', 14);
ylabel("Iteration time [s]", 'FontSize', 14);
saveas(gca, 'cores_itertime_768.eps', 'eps');

plot(cores, time(1)./time, '-kx', cores, speedup, '-rx');
xlabel('Number of cores []', 'FontSize', 14);
ylabel('Speedup []');
legend('Empirical scaling', 'Ideal scaling');
saveas(gca, 'speedup_768.png', 'png');

plot(cores, per, '-kx');
xlabel('Number of cores []', 'FontSize', 14);
ylabel("Time per iteration [s]", 'FontSize', 14);
saveas(gca, 'cores_timeper_768.eps', 'eps');

plot(cores, gather, '-kx');
xlabel('Number of cores []', 'FontSize', 14);
ylabel('Gathering time [s]', 'FontSize', 14);
saveas(gca, 'cores_gather_768.eps', 'eps');