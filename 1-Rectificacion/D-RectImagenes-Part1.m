% sample an example Aerielle movie.  Set this up as generic for other
% movies 

% The important things are:
%   - the init material below.  Stored as a structure so I can preserve
%       recordsfor each UAV analysis.
%   - a gcp file from a current survey
%   - 

clear
close all

addpath(genpath('.')) % adds folder dependencies 

% required USER INPUT material.  Adapt to each collection
demoInputFile_12V1P1;      %Modifique de acuerdo al vuelo correspoendiente this file contains all of the required inputs.

% create cx directory and load instruments
if ~exist(inputs.pncx)
    mkdir(inputs.pncx)
end


[dname, fname, ~]=fileparts(inputs.instsFn);
if isempty(dname)
    insts = eval(fname);
else
    p = pwd;
    cd(dname);
    insts = eval(fname);
    cd(p);
    clear p;
end

for i = 1: length(insts)        % save inst info as part of stacks
    stack(i).inst = insts(i);
end
load(inputs.gcpFn);             % load gcps

% find the input folders (may be more than one)
e0 = matlab2Epoch(inputs.dn0);
dum=[inputs.pnIn filesep inputs.frameFn,'*'];
clipFns = dir(dum);
NClips = length(clipFns);
for i = 1: NClips
    fns{i} = dir([inputs.pnIn filesep clipFns(i).name filesep '*png']);
     if isempty(fns{i})
         fns{i} = dir([inputs.pnIn filesep clipFns(i).name filesep '*jpg']);
     end
    Nf(i) = length(fns{i});
end
nt = sum(Nf);
dn = datenum(inputs.dateVect) + ([1:nt]-1)*inputs.dt;

% read the first frame and display.  Do a manual geometry on it if needed.
I = imread([inputs.pnIn filesep clipFns(1).name filesep fns{1}(1).name]);
[NV, NU, NC] = size(I);
Ig = rgb2gray(I);           % for later sampling.
meta.showFoundRefPoints = inputs.showFoundRefPoints; % easy way to pass.
meta.showInputImages = inputs.showInputImages;

% Because nlinfit requires globals, we set up a variable globals (under 
% metadata) that contains the lcp as well as the flags and values for
% known geometry variables (eg if we know the camera roll).  To minimize
% the use of globals we pass this explicitly except for nlinfit uses when
% we declare it to be global within the calling subroutine.
meta.globals.lcp = makeLCPP3(inputs.stationStr,NU,NV);
meta.globals.knownFlags = inputs.knownFlags;
meta.globals.knowns = inputs.knowns;   
clear NU NV NC 

% When geometries are solved for each frame in turn, the results are saved
% in a metadata file in the cx folder for this day.  First search to see if
% such a folder already exists, in which we don't have to re-do all of the
% geometries.  Allow a few second slop in time (look for first 9 out of 10
% digits in epoch time.

bar = num2str(e0);
foo = dir([inputs.pncx bar(1:9) '*meta*']);
% report the use of an old metadata file, ask for ok
if ~isempty(foo)
    yn = input(['Found old metadata file ' foo(1).name '\n' ...
        'Hit return to continue or anything else to start over...'], ...
        's');
    if ~isempty(yn)
        foo = [];
    end
end
if ~isempty(foo)
    oldGeoms = 1;           % flag that old geoms exist and load
    load([inputs.pncx filesep foo(1).name]); 
    betas = meta.betas;
else
    % if no metafile is found, I initialize one
    oldGeoms = 0;
    betas = nan(nt,6);              % we will save all the betas
    [betas(1,:),meta.refPoints] =initUAVAnalysis(I, gcp, inputs, meta);
end


%%

save('./Outputs/OutputD-Part1')






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
