% These are the run times, on my computer, for the peptide mass
% inversion problem. I ran these in such a way that they always
% returned only 1 answer.
function plotRunTimes

Data = [      8   0.092
              9   0.260
              10   0.723
              11   1.914
              12   5.063
              13  12.562
              14  28.597
              15  66.162
              16 139.579];

useNamedFigure('PeptideRunTimes'); clf;
plot(Data(:,1),Data(:,2));
xlabel('Peptide Length');
ylabel('Run Time');

Ratios = Data((2:end),2)./Data(1:(end-1),2)
title(sprintf('Run Time For Peptide Composition Search: Ratio-%.3f', ...
              mean(Ratios)));

prettyPlot;
print('-dpng','RunTimes.png');
