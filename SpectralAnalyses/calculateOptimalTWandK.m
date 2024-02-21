function [TW, K, T] = calculateOptimalTWandK(data, Fs)
    % Calculate the duration of the data
    numTimePoints = size(data, 2); % Assuming columns represent time points
    T = numTimePoints / Fs; % Duration in seconds
    
    % Define your desired time-bandwidth product here
    % This is a starting point; you may need to adjust based on your specific needs
    TW = 3; % Time-bandwidth product, a common starting point
    
    % Calculate the number of tapers
    % K is often chosen as 2*TW - 1 for good frequency resolution and variance trade-off
    K = 2*TW - 1;
    
    % Display the results
    fprintf('Optimal Time-Bandwidth Product (TW): %f\n', TW);
    fprintf('Number of Tapers (K): %d\n', K);
    fprintf('Duration (T) in seconds: %f\n', T);
end

