% Read the input file
fileID = fopen('output_1300.txt', 'r');
text = fscanf(fileID, '%c');
fclose(fileID);

% Regular expressions to find iteration numbers and quorum sizes
iterationPattern = 'Iteration:\s*(\d+)';
quorumPattern = 'quorum\s*size:\s*(\d+)';

% Extract iteration numbers
iterationMatches = regexp(text, iterationPattern, 'tokens');
iterationNumbers = cellfun(@str2double, [iterationMatches{:}]);
iterationNumbers = iterationNumbers+1;
% Extract quorum sizes
quorumMatches = regexp(text, quorumPattern, 'tokens');
quorumSizes = cellfun(@str2double, [quorumMatches{:}]);
quorumSizes = quorumSizes-1;

% Display the extracted iteration numbers and quorum sizes
figure;
plot(iterationNumbers,quorumSizes);
yline(130, 'b--', {'required = 130'}, 'LineWidth', 1.5);
% Labels and title
xlabel('Iteration');
ylabel('Quorum Size');
title('Pop: 1300, Q: 1000');
%legend('Histogram of qmp values', 'Logistic pdf');
grid on;
%disp('Iteration Numbers:');
%disp(iterationNumbers);
%disp('Quorum Sizes:');
%disp(quorumSizes);
