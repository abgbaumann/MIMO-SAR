function [aoi] = round_aoi(aoi, digits)
% A function to round a cell with matrices to a user-defined number of
% digits.
%
% 24th December 2020 / ABA
%

if ~exist('digits','var')
    digits = 3;
end

for aoi_i = 1:length(aoi)
    aoi{aoi_i} = round(aoi{aoi_i},digits);
end

end

