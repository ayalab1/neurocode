function a_prime = regress_out(a, b)
% regress_out Regress variable b from variable a while preserving the original mean of a.
%
%   adapted from the seaborn function of the same name
%       https://github.com/mwaskom/seaborn/blob/824c102525e6a29cde9bca1ce0096d50588fda6b/seaborn/regression.py#L337
%
%   a_prime = regress_out(a, b) performs regression of variable b from variable a
%   while preserving the original mean of a. The function calculates the residual
%   component of a that remains after removing the effect of b. The regression is
%   performed using the ordinary least squares method.
%
%   Inputs:
%   - a: Input array representing the main variable of interest.
%   - b: Input array representing the variable to be regressed out from a.
%
%   Output:
%   - a_prime: Output array representing the residual component of a after removing
%              the effect of b.
%
%   Example:
% Generate random data
% rng(1);  % Set the random seed for reproducibility
% n = 100;  % Number of samples
% a = randn(n, 1);  % Random normal distribution for variable a
% 
% % Create correlated variable b based on a
% correlation_coefficient = 0.7;  % Specify the desired correlation coefficient
% b = correlation_coefficient * a + sqrt(1 - correlation_coefficient^2) * randn(n, 1);
% 
% % Perform regression while preserving the mean of a
% a_prime = regress_out(a, b);
% 
% figure;
% subplot(1,2,1)
% scatter(a, b, 'b', 'filled');
% xlabel('value a');
% ylabel('value b');
% subplot(1,2,2)
% scatter(a_prime, b, 'b', 'filled');
% xlabel('a prime');
% ylabel('value b');
% 
% % Plot the original data and the regressed component
% figure;
% scatter(1:n, a, 'b', 'filled');
% hold on;
% scatter(1:n, a_prime, 'r', 'filled');
% xlabel('Sample');
% ylabel('Value');
% legend('Original Data (a)', 'Regressed Component (a\_prime)');
% title('Regression of Variable b from Variable a');
% 
% % Plot the distribution of the regressed component
% figure;
% histogram(a_prime, 'Normalization', 'pdf');
% xlabel('Value');
% ylabel('Density');
% title('Distribution of Regressed Component (a\_prime)');
%
%
%   Note:
%   - The function assumes that a and b are numeric arrays with the same size.
%
%   See also: pinv, mean
%
% Ryan Harvey 2023

a_mean = mean(a);
a = a - a_mean;
b = b - mean(b);
a_prime = a - b * pinv(b) * a;
a_prime = reshape(a_prime + a_mean, size(a));
end