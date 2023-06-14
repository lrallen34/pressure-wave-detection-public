%%% Function for making currently non-existent directories
function makeneededdirs(varargin)
% input(s): 
%   directories (str): as many directories as desired to (1) check if they
%   exist and (2) if not, create them

for curdir = 1:length(varargin)
    if ~exist(varargin{curdir}, 'dir')
        mkdir(varargin{curdir})
    end
end
end
