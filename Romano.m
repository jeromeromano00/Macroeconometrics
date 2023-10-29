% Macroeconometrics Problem set (1)


% Exercise 1
% From FREDII data base (http://research.stlouisfed.org/fred2/) 
% download the series GDPC1 (quarterly US real GDP). Transform the series
% in growth rates. Suppose that real GDP growth rates follow an AR(1)
                        % y_t = c + phi*y_{t-1} + e_t
                        % e ~ iid(o,sigma^2)

GDPC1 = readtable('GDPC1.csv')
%folderPath = 'C:\Users\Jerome Romano\OneDrive\Desktop\Macroeconometrics'
%% (a) Plot the series

% we import manually the data from FREDII

% Extract the GDP values
gdp_values = GDPC1.GDPC1 ;

% Calculate the growth rates
gdp_growth = diff(gdp_values) ./ gdp_values(1:end-1) * 100;

% Create an array of quarters for plotting
quarters = 1:numel(gdp_growth);

% Plot the GDP growth rates
plot(quarters, gdp_growth);

% Add labels and title to the plot
xlabel('Quarters');
ylabel('GDP Growth Rate (%)');
title('US Real GDP Growth Rate');

% Add a grid for better visualization
grid on;

%filename = fullfile(folderPath, 'first_image.png');  % Change the filename and format as needed

% Save the figure as an image file
%saveas(gcf, filename);

%% (b) Estimate the parameters c and φ with OLS.

% Lag the GDP growth data by one and two quarters - this is done in order
% to achieve harmonization in the length of the matrix
lag1_gdp_growth = lagmatrix(gdp_growth, 1);
lag2_gdp_growth = lagmatrix(gdp_growth, 2);

% Extract lagged values of GDP growth (excluding the first two rows)
y_tminus2 = gdp_growth(3:end);     % y(t-2)
y_tminus1 = lag1_gdp_growth(3:end); % y(t-1)
y_t = lag2_gdp_growth(3:end);      % y(t)

% Create the independent variable matrix 'x' using lagged values
x = [y_tminus1, y_tminus2];

% Fit a linear regression model using the lagged values
model = fitlm(x, y_t);

% Estimate an AR(1) model using the 'ar_estimation' function
[c, phi] = ar_estimation(1, gdp_growth);

disp('Value c:');
disp(c);
disp('Value phi:');
disp(phi);

%% (c) Estimate and plot the first 10 autocorrelations.

% Define the number of lags
maxLag = 10;

% Preallocate an array to store autocorrelation values
autocorrelations = zeros(1, maxLag + 1);

% Calculate the autocorrelation values based on the AR(1) model
autocorrelations(1) = 1; % Autocorrelation at lag 0 is always 1
for lag = 1:maxLag
    autocorrelations(lag + 1) = phi ^ lag;
end

% Plot the autocorrelations
figure;
stem(0:maxLag, autocorrelations);
xlabel('Lag');
ylabel('Autocorrelation');
ylim([min(-0.1, min(autocorrelations) - 0.05), max(0.2, max(autocorrelations) + 0.05)]);
title('Autocorrelations Based on AR(1) Model Coefficients');

% Display the autocorrelation values
fprintf('Autocorrelation values for lags 0 to %d:\n', maxLag);
for lag = 0:maxLag
    fprintf('Lag %d: %f\n', lag, autocorrelations(lag + 1));
end

%filename = fullfile(folderPath, 'second_image.png');
%saveas(gcf, filename);


%% D) Obtain and plot the coefficients of the Wold representation of yt

% Define the AR(1) process parameters (replace with your estimated values)
% phi: Autoregressive coefficient (estimated value)

% Define the maximum lag for the coefficients
maxLag = 10;

% Create a unit impulse at lag 0
impulse = [1; zeros(maxLag, 1)];

% Use the filter function to obtain the Wold coefficients
WoldCoefficients = filter(1, [1, -phi], impulse);

% Plot the Wold coefficients
figure;
stem(0:maxLag, WoldCoefficients);
xlabel('Lag');
ylabel('Wold Coefficient');
ylim([min(WoldCoefficients) - 0.1, max(WoldCoefficients) + 0.1]);
title('Wold Representation Coefficients for yt (AR(1) Process)');

% Display the Wold coefficients
fprintf('Wold Coefficients for AR(1) Process:\n');
for lag = 0:maxLag
    fprintf('Lag %d: %f\n', lag, WoldCoefficients(lag + 1));
end

%filename = fullfile(folderPath, 'third_image.png');
%saveas(gcf, filename);

%% 2.Now suppose that GDP admits an AR(2) representation
%         y_t = c + φ_{1} y_{t−1} + φ_{2} y_{t−1} + ε

%(a) Estimate the parameters c, φ_1, φ_2 with OLS.

[c,phi] = ar_estimation(2,gdp_growth);
fprintf('Estimated Mean: %f\n', c);
fprintf('Estimated Autoregressive Coefficient (phi1): %f\n', phi(1));
fprintf('Estimated Autoregressive Coefficient (phi2): %f\n', phi(2));
%% (b) Compute the roots of the AR polynomial

%see that the polynomial is given by \phi(z) = 1 - 1.2z - 0.10L^2 = 0 we can
%use the root function of matlab to compute the roots of that polynomial,

% Define the coefficients a, b, and c
a = -phi(2,1);
b = -phi(1,1);
c = 1;

% Solve the quadratic equation
roots_1 = roots([a, b, c]);
abs_roots = abs(roots_1);

abs_roots % since abs_roots are outside the unit circle, the process is stationary
%% (c) State the conditions under which the process is causal and stationary.
 % see pdf file

 %% (d) Obtain and plot the coefficients of the Wold representation of y_t

% Define the AR(2) process parameters (replace with your estimated values)

% Define the maximum lag for the coefficients
maxLag = 10;

% Create a unit impulse at lag 0
impulse = [1; zeros(maxLag, 1)];

% Use the filter function to obtain the Wold coefficients for the AR(2) process
AR2Coefficients = [1, -phi(1), -phi(2)];
WoldCoefficients = filter(1, AR2Coefficients, impulse);

% Plot the Wold coefficients
figure;
stem(0:maxLag, WoldCoefficients);
xlabel('Lag');
ylabel('Wold Coefficient');
ylim([min(-0.1), max(0.2)]);
title('Wold Representation Coefficients for yt (AR(2) Process)');

% Print the Wold coefficients
fprintf('Wold Representation Coefficients for AR(2) Process:\n');
for lag = 0:maxLag
    fprintf('Lag %d: %f\n', lag, WoldCoefficients(lag + 1));
end

%filename = fullfile(folderPath, 'fourth_image.png');
%saveas(gcf, filename);

