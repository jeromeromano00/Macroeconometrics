function [c_hat, phi] = ar_estimation(p, data)
    % Function to estimate an autoregressive (AR) model of order p.

    % Get the length of the data
    n = length(data);

    % Initialize a matrix to store lagged data
    lagged_data = zeros(n, p);

    % Calculate the lagged values of the data
    for i = 1:p
        lagged_data(:, i) = lagmatrix(data, i);
    end

    % Remove rows with insufficient lagged data
    lagged_data = lagged_data(p+1:end, :);

    % Split the lagged data into predictor variables (X) and the response variable (y)
    X = lagged_data;
    y = data(p+1:end);

    % Fit a linear regression model using the lagged data
    model_automatic = fitlm(X, y);

    % Extract the estimated coefficients from the model
    c_hat = model_automatic.Coefficients.Estimate(1);  % Intercept
    phi = model_automatic.Coefficients.Estimate(2:end); % AR coefficients
end
