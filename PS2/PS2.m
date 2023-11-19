% Read data from the 'GS10.csv', 'FEDFUNDS.csv', and 'GDPC1.csv' files
ten_gov_bond = readtable('GS10.csv');
fed_fund_rate = readtable('FEDFUNDS.csv');
real_GDP = readtable('GDPC1.csv');

% Specify the key variable for the join operation
keyVariable = 'DATE';

% Perform the join operation on 'DATE' between fed_fund_rate and ten_gov_bond
merged_data = join(fed_fund_rate, ten_gov_bond, 'Keys', keyVariable);

% Calculate the 'spread' by subtracting 'FEDFUNDS' from 'GS10'
merged_data.spread = merged_data.GS10 - merged_data.FEDFUNDS;

% Create a new table with 'DATE' and 'spread'
spreadTable = table(merged_data.DATE, merged_data.spread, 'VariableNames', {'DATE', 'spread'});

% Display the result
disp(spreadTable);

% Find the earliest common date in both tables
minDate = max(min(spreadTable.DATE), min(real_GDP.DATE));

% Find the latest common date in both tables
maxDate = min(max(spreadTable.DATE), max(real_GDP.DATE));

% Select data between the earliest and latest common dates
spreadTable = spreadTable(spreadTable.DATE >= minDate & spreadTable.DATE <= maxDate, :);
real_GDP = real_GDP(real_GDP.DATE >= minDate & real_GDP.DATE <= maxDate, :);

% Convert the table to a timetable
spreadTable = table2timetable(spreadTable);

% Resample the data to quarterly frequency, calculating the mean
spreadTable = retime(spreadTable, 'quarterly', 'mean');

% Rename the 'spread' variable to 'quart_spread'
spreadTable = renamevars(spreadTable, 'spread', 'quart_spread');

% Display the resampled and renamed data
disp(spreadTable);

% Calculate the percentage growth of GDP using the 'diff' function
real_GDP.gdp_growth = [0; diff(real_GDP.GDPC1) ./ real_GDP.GDPC1(1:end-1) * 100];

%%


% estimate VAR(4)
data = [real_GDP.gdp_growth, spreadTable.quart_spread];
p = 4;
[B,Y,X]=VAREstim(data,p);
n=size(data,2);

%compute and plot the Wold representation

% 1) Companion form F
F=[B(2:end,:)'; eye(n*(p-1)) zeros(n*(p-1),n)];

% 2) Powers of F
for j=1:24
    FF=F^(j-1);
    C(:,:,j) = FF(1:n,1:n);
end

% 3) Take the appropriate submatrix from F to compute the imput response
% function

ii=0;
for j=1:n
    for i=1:n
        ii=ii+1;
        subplot(n,n,ii),plot(squeeze(C(j,i,:)));
        title ("Effect of error " + num2str(i) + " on variable " + num2str(j))
    end
end   
