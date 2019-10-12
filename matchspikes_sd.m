function [matched missing] = matchspikes_sd(groundtruth, clusterspikes, varargin)
% match the ground truth (1 * N vector) with the clustered spikes (2 * M
% vector). Matched means that there was at least one spike within a certain
% time limit, default being 1ms (variable using varargin: 'MinTime'
%
% returns the matching and the number missing
%
glength = length(groundtruth) ;
 temp = size(clusterspikes) ;

mintime = 0.003 ; % 1ms default 
missing = 0 ;
matched = 0 ;
% for i = 1: temp
%     matched(i) = i-1 ;
% end

for i=1:(nargin-2)
    if ischar(varargin{i})
        switch varargin{i}
            case 'MinTime'
                mintime = varargin{i+1};
        end
    end
end

% note that clusterspikes has time in ms, but groundtruth in seconds
% clusterspikes = clusterspikes/1000 ;

ci = 1 ; % clusterindex

for gi=1:glength % ground  truth index

    pindex = 1 ;
    possibles = [] ;
    while clusterspikes(ci) < (groundtruth(gi) - mintime)
        ci = ci + 1 ;
        if (ci > temp) break; end
    end
    if (ci > temp) ci = temp ; end
    while abs(clusterspikes(ci) - groundtruth(gi)) <= mintime
        possibles(pindex) = clusterspikes(ci) ;
        pindex = pindex + 1 ;
        ci = ci + 1 ;
        if (ci > temp) break; end
    end
    if (isempty(possibles))
        missing = missing + 1 ;
    else
        % find element of possibles with minimum distance
%         for j = 1:pindex-1
%             best = 1 ;
%             if abs(possibles(j) - groundtruth(gi)) < abs(possibles(best) - groundtruth(gi))
%                 best = j ;
%             end
%         end
        matched= matched+ 1 ;
    end
    % reset ci to something reasonable
    ci = max(1, ci-10) ;
end
        

    