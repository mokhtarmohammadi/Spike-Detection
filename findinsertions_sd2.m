function insertions = findinsertions_sd2(actualpeaks, clusterspikes, varargin)
% find any output spikes which are outside of some minimum time from any of
% the ground truth data
% Minimum time default is 0.001 seconds. Alterable (mintime).
% Needs to use the structure returned by the spike generating system since
% it needs access to all the ground truth data
%
% form a single array from all the ground truth data , and sort it into
% time order
mintime = 0.0001 ;
for i=1:(nargin-2)
    if ischar(varargin{i})
        switch varargin{i}
            case 'MinTime'
                mintime = varargin{i+1};
        end
    end
end

% note that clusterspikes has time in ms, but groundtruth in seconds
% clusterspikes(:,2) = clusterspikes(:,2)/1000 ;

% numgt = length(groundtruthstruct) ;
gtall = [] ;
% for i = 1:numgt
%     gtall = [gtall groundtruthstruct(i).actualpeaks] ;
% end
gtall=actualpeaks;
gtall = sort(gtall) ;
%
glength = length(gtall) ;
cllength = size(clusterspikes,2) ;
insertions = 0 ;
%
gi = 1 ;
for ci = 1:cllength
    while gtall(gi) < (clusterspikes(ci) - mintime)
        gi = gi + 1 ;
        if gi > glength break; end
    end
    if gi > glength
        gi = glength ;
    end
    if abs(clusterspikes(ci) - gtall(gi)) <= mintime
        ;
    else
        insertions = insertions + 1 ;
        
    end
    % reset gi to something reasonable
    gi = max(1, gi-10) ;
end
