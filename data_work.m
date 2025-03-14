
% Read the CSV file (replace 'yourfile.csv' with your actual file name)
data = readmatrix('iteration_382_pop_1300.csv');
% Logistic PDF
logisticPDF = @(x, mu, s) exp(-(x - mu)/s) ./ (s * (1 + exp(-(x - mu)/s)).^2);

% Logistic CDF
logisticCDF = @(x, mu, s) 1 ./ (1 + exp(-(x - mu)/s));


column2 = data(2:end, 2);
[histValues, edges] = histcounts(column2, 'Normalization', 'pdf');
%disp(histValues);
disp(S);
S = sum(histValues);
%mu = mean(column2);
mu = mean(column2);
disp(mu);
%mu = mean(histValues);
% Estimate the scale (s)
s = std(column2) * sqrt(3) / pi;
%s = std(histValues) * sqrt(3) / pi;
x = linspace(0, max(column2), 1000);
%x = linspace(min(histValues), max(histValues), 1000);
% Calculate the PDF
pdf_values = logisticPDF(x, mu, s);

% Calculate the CDF
cdf_values = logisticCDF(x, mu, s);

binWidth = edges(2) - edges(1); % Adjust the bin width if needed
%histogramArea = binWidth * sum(pdf_values);
%pdf_values = pdf_values;

%[histValues, edges] = histcounts(column2, 'Normalization', 'pdf');
binCenters = (edges(1:end-1) + edges(2:end))/2;
maxHistValue = max(histValues);
%maxHistValue = max(histValues);
pdf_values = (pdf_values/max(pdf_values))*maxHistValue;


figure;
% Create histogram plot for column 2
histogram(column2, 'Normalization', 'pdf', 'BinWidth', binWidth);
hold on;

% Plot the PDF
plot(x, pdf_values, 'r', 'LineWidth', 2);
xline(0.5, 'b--', {'Threshold = 0.5'}, 'LineWidth', 1.5);
xline(mu, 'g--', {'mean'}, 'LineWidth', 1.5);
% Labels and title
xlabel('QMP level');
ylabel('Node count x100');
title('Pop: 1300, Q: 1000, iter: 382');
%legend('Histogram of qmp values');
legend('Histogram of qmp values', 'Logistic pdf');
hold off;
grid on;
% %%
% % Create cumulative plot for column 2
% figure;
% [counts2, edges2] = histcounts(column2, 'Normalization', 'cdf');
% cdf_centers = (edges2(1:end-1) + edges2(2:end)) / 2;
% cdf2 = cumsum(counts2);
% %disp(cdf_centers);
% %disp(counts2);
% subplot(2,1,1);
% plot(cdf_centers, counts2, '-r');
% xlabel('QMP Values');
% ylabel('Cumulative Distribution');
% title('trace of the cumulative histogram');
% 
% subplot(2, 1, 2);
% histogram(column2, 'Normalization', 'cdf');
% hold on;
% % Plot the PDF
% plot(x, cdf_values, 'r', 'LineWidth', 2);
% xlabel('QMP Values');
% ylabel('Cumulative Distribution');
% title('Cumulative Histogram of values with logistic cdf');
% grid on;
% %%/
% Compute the cumulative distribution (CDF) for column 2
% Compute 1 - CDF values
one_minus_cdf = 1 - counts2;

% Plot 1 - CDF for column 2
%subplot(3, 1, 3);
%plot(cdf_centers, one_minus_cdf, '-r');
%xlabel('QMP Values');
%ylabel('1 - Cumulative Distribution');
%title('1 - CDF of values');