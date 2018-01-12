%DISPARITY  High contrast colormap with subtle gradient discontinuities
%
% Examples:
%   map = disparity;
%   map = disparity(len);
%   B = disparity(A);
%   B = disparity(A, lims);
%
% A colormap designed for depth and disparity maps. It has subtle gradient
% discontinuities to bring out contours, and maximizes the range of colors
% used in order to improve contrast between intensity levels. It converts
% linearly to grayscale, such that black & white prints come out nicely.
%
% The function can additionally be used to convert a real valued array into
% a truecolor array using the colormap.
%
% IN:
%   len - Scalar length of the output colormap. If len == Inf the concise
%         table is returned. Default: len = size(get(gcf, 'Colormap'), 1);
%   A - Non-scalar numeric array of real values to be converted into
%       truecolor.
%   lims - 1x2 array of saturation limits to be used on A. Default:
%          [min(A(:)) max(A(:))].
%
% OUT:
%   map - (len)xJ colormap table. J = 3, except in the concise case, when
%         J = 4, map(1:end-1,4) giving the relative sizes of the 
%         inter-color bins.
%   B - size(A)x3 truecolor array.

% Copyright: Oliver Woodford, 2012

function map = disparity(varargin)
map = [            0                   0                   0   1.000000000000000; ...
                   0                   0   0.548245614035088                   0; ...
   0.105761794654448   0.000000000000013   0.270852836827302   1.000000000000000; ...
   0.059459459459459                   0   0.940540540540541                   0; ...
   0.162903903427862   0.002287878882682   0.657445157639783   1.000000000000000; ...
   0.366693676801250   0.059947169013601   0.374294758117918                   0; ...
   0.397297297297297                   0   0.602702702702703   1.000000000000000; ...
   0.735135135135135                   0   0.264864864864865                   0; ...
   0.619346855510299   0.046614971192608   0.328528965897891   1.000000000000000; ...
   1.000000000000000                   0   0.118421052631579                   0; ...
   0.878236661421783   0.026878650964249   0.299381316832219   1.000000000000000; ...
   0.879983391744831   0.020603479222390   0.875357224252214                   0; ...
   1.000000000000000                   0   0.666666666666667   1.000000000000000; ...
   0.859195402298851   0.140804597701150   0.859195402298851                   0; ...
   0.796752292660960   0.203005124169332   0.702693479008556   1.000000000000000; ...
   0.500000000000000   0.500000000000000   0.500000000000000                   0; ...
   0.591584728266931   0.512311475905112   0.196397630630585   1.000000000000000; ...
   0.140804597701149   0.859195402298851   0.140804597701149                   0; ...
   0.237154669316016   0.801555080604367   0.184894048769717   1.000000000000000; ...
                   0   1.000000000000000   0.333333333333334                   0; ...
   0.099802756002636   0.943534783080114   0.362316300764779   1.000000000000000; ...
                   0   1.000000000000000   0.881578947368421                   0; ...
   0.122068579355955   0.967299920607226   0.729793345404631   1.000000000000000; ...
   0.264864864864865   1.000000000000000   0.735135135135135                   0; ...
   0.371207986455179   0.946429224305391   0.732060152479271   1.000000000000000; ...
   0.602702702702703   1.000000000000000   0.397297297297297                   0; ...
   0.626666491156471   0.941072277957457   0.637870982308667   1.000000000000000; ...
   0.940540540540541   1.000000000000000   0.059459459459459                   0; ...
   0.842450533320871   0.992890454106401   0.353338543654404   1.000000000000000; ...
   0.894238205345540   0.999999999999999   0.729147163172666                   0; ...
   1.000000000000000   1.000000000000000   0.451754385964912   1.000000000000000; ...
   1.000000000000000   1.000000000000000   1.000000000000000                   0];
map = colormap_helper(map, varargin{:});