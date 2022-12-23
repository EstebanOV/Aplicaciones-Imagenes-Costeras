function insts = makeDJIInsts201510011529
%   insts = makeDJIInsts
%
% creates pixel instruments for DJI video.  Types can be line or matrix
% this one is empty - just test image products.
cnt = 1;

% vBar instruments
y = [-200 200];
x = [200: 25: 300];
z = 0;
for i = 1: length(x)
    insts(cnt).type = 'line';
    insts(cnt).xyz = [x(i) y(1) z; x(i) y(2) z];
    eval(['insts(cnt).name = ''vBar' num2str(x(i)) ''';']);
    eval(['insts(cnt).shortName = ''vBar' num2str(x(i)) ''';']);
    cnt = cnt+1;
end

% some runup lines
x = [150 350];
y = [-200:50:0];
z = 0;
for i = 1: length(y)
    insts(cnt).type = 'line';
    insts(cnt).xyz = [x(1) y(i) z; x(2) y(i) z];
    eval(['insts(cnt).name = ''runup' num2str(y(i)) ''';']);
    eval(['insts(cnt).shortName = ''runup' num2str(y(i)) ''';']);
    cnt = cnt+1;
end

% cBathy array
x = [150 5 350];   % determine sample region and spacing
y = [-250 5 200];    % format is [min del max]
z = 0;
insts(cnt).type = 'matrix';
insts(cnt).name = 'cBathyArray';
insts(cnt).shortName = 'mBW';
insts(cnt).x = x;
insts(cnt).y = y;
insts(cnt).z = z;
cnt = cnt+1;

% make some slices to check stability
insts(cnt).type = 'line';
insts(cnt).xyz = [300 140 7; 300 100 7];
insts(cnt).name = 'x = 300 pier transect';
insts(cnt).shortName = 'x300Slice';
cnt = cnt+1;

insts(cnt).type = 'line';
insts(cnt).xyz = [100 220 3; 115 220 3];
insts(cnt).name = 'y = 220 Piling x-transect';
insts(cnt).shortName = 'y220Slice';


%
%   Copyright (C) 2017  Coastal Imaging Research Network
%                       and Oregon State University

%    This program is free software: you can redistribute it and/or  
%    modify it under the terms of the GNU General Public License as 
%    published by the Free Software Foundation, version 3 of the 
%    License.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.

%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see
%                                <http://www.gnu.org/licenses/>.

% CIRN: https://coastal-imaging-research-network.github.io/
% CIL:  http://cil-www.coas.oregonstate.edu
%
%key UAVProcessingToolbox
%

