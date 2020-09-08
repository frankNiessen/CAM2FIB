function CAM2FIB
%% CAM2FIB.m
% Interpretation and interpolation of g-code to coordinate list + creation
% of a stream file for patterning with a FIB/SEM instrument
% (Optimized for use with FEI Helios Nanolab)
% *************************************************************************
% Copyright (c) 2019, Frank Niessen, EMC, AIIM, University of Wollongong
% All rights reserved
% If you find this program useful and apply it for your research I would
% appreciate a citation to the associated research paper
% [F. Niessen, M.J.B. Nancarrow, Computer-aided manufacturing and focused
% ion beam technology enable machining of complex micro- and nano-structures,
% Nanotechnology, 2019, https://doi.org/10.1088/1361-6528/ab329d]
% *************************************************************************
close all; clc
fprintf('*****************************************************************\n')
fprintf('                         CAM2FIB.m \n');
fprintf('*****************************************************************\n')
%% Initialization
checkCurrentFolder;
fprintf('\n*** Initializing');
% *** File Input settings *************************************************
[f.Name,f.Path] = uigetfile([fileparts(mfilename('fullpath')),...
                            '\Input_Gfiles\*.*'],'Select a G-code file');  % Select G-Code file
fNameFull = [f.Path,'\',f.Name];                                           % Assemble full filename
% *** File settings **********************************************
f.unitLabel = 'nm';                                                        % Output length unit - 'm' 'mm' '�m' 'nm' or '�'
f.unitOut = convertUnit(f.unitLabel);                                      % Conversion of Output unit to m
f.overlap = 50;                                                            % Relative overlap of beam path [%]
f.dcrit = 0;                                                               % Minimum length for G0/beam-blank repositioning
f.binFac = 1;                                                              % Interpolation binning factor of entire beam path: >1: Not active; >=1: Reinterpolation with f.binFac -> Recommended: 1
% *** Visual Output Options ***********************************************
Out.scrPrint = 0;                                                          % Flag: Screen print
Out.plot = 1;                                                              % Flag: Plotting output
Out.scatter = 1;                                                           % Flag: Scatter points instead of lines
Out.plotCnt = 100;                                                         % Nr of updates for live plot
% *** Stream file output settings *****************************************
str.DAC = 's16';                                                           % Digital Analog Converter (DAC) type ['s16': DAC is 16 bits]
str.nrIter = 1;                                                            % Nr of iterations for pattern
str.tMach = [50 100];                                                      % Pattern machining time [s] / Multiple times possible [t1 t2 t3 t4 ...]
str.fType = 'str';                                                         % Output file extension
% *** Scaling parameters - Patterning-Imaging-Acquisition (PIA) ***********
PIA.rot = 0;                                                               % Clockwise pattern rotation [�]
PIA.flip = [0 0];                                                          % Flag: Flip [x y]
PIA.scalFac = 1;                                                           % Scaling factor (linear scaling of beam path - default [1])
PIA.cal = 0.31618*1e6;                                                     % Bits per m per Mag.
PIA.max = [65536 56576];                                                   % Adressable max x and y range of Patterning-Imaging-Acquisition
PIA.offs = [0 PIA.max(1)-PIA.max(2)];                                      % Offset of x and y range of PIA (Leave as [0 0] by default)
PIA.AR = PIA.max(2)/PIA.max(1);                                            % PIA aspect ratio
PIA.relScale = 0.5;                                                        % Relative range of points in 'PIA.max' used for patterning [0 to 1]
PIA.centrePattern = 1;                                                     % Flag: Centre pattern in Pattering range
PIA.availMags = [.1e3 .5e3  1e3  2e3  5e3  10e3  20e3  30e3  40e3 ...
                50e3 60e3];                                                % Available magnifications
%% Processing
[pos,bb.raw,opt,h] = gCode2beamPath(fNameFull,Out,f);                      % Interpolate beampath from gCode
[pos.BPinterp,bb.interp] = interpBeamPath(pos.BP(:,1:2),bb.raw,f,opt);     % Filter beamPath
[pos,str] = scaleBeamPath(pos,PIA,f,str);                                  % Scale beamPath
%% Plotting
% *** Create machining strategy plot **************************************
if Out.plot
    plotMaps(h,pos,bb,f);
    tileFigs();
    drawnow;
end
%% Writing Stream files
fprintf('\n*** Output of Stream files');
for i=1:length(str.tMach)
    fprintf('\n   -> Writing streamfile %.0f of %.0f: Mag.: %.0fx t_Mill: %.0f s',i,length(str.tMach),str.patMag,str.tMach(i));
    strFilePath = writeStreamFile(pos,bb,str,f,str.tMach(i));      % Write stream pattern file
end
fprintf('\nAll done!\n');
end
%% Read in gCode file
function out = readGfile(filename, startRow, endRow)
% *** Initialize variables ************************************************
delimiter = ' ';
if nargin<=2
    startRow = 1;
    endRow = inf;
end
% *** Open the text file **************************************************
formatSpec = '%s%s%s%s%s%s%s%[^\n\r]';
fileID = fopen(filename,'r');
% *** Read columns of data according to the format ************************
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'EmptyValue' ,NaN,'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'EmptyValue' ,NaN,'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end
% *** Close the text file and create output variable **********************
fclose(fileID);
out = [dataArray{1:end-1}];
end
%% Interpolate beamPath from gCode  commands
function [pos,bb,opt,h] = gCode2beamPath(fName,Out,f)
%function beamPath = gCode2beamPath(fName,dRes,angRes,Out,f)
% Read in of g-code file and output of interpolated beam path cooredinates
% for a given distance and angular resolution
% fName:    File name, either relative to current folder or including
%           absolutepath
% dRes:     Distance resolution in input unit
% scrPrint: Screen output
% beamPath: x, y, and z coordinate list of interpolated beampath
% ************************************************************************
% This code curently supports commands G0, G1, G2 and G3 for rapid
% movement, rapid cordinated movement, clockwise arc move and counter
% clockwise arc move
% ************************************************************************
% Adapted from an initial version by Tom Williamson (18/06/2018)
% [https://mathworks.com/matlabcentral/fileexchange/67767-g-code-reader]
% Original copyrright:
% Copyright (c) 2018, Tom Williamson
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% * Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
%
% * Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% *** Ini *****************************************************************
opt.absIJ = 0;                                                             %Absolute [1] or relative [0] coordinates I and J for arc center of G2 and G3
opt.preAllocFac = 10;                                                      %Preallocation factor for output array
opt.dcrit = f.dcrit;                                                       %Minimum distance for 'G0' beam blanked movements (leave empty '[]' if no threshould should be applied)
%opt.dRes = f.dRes;                                                        %Linear resolution in length unit
opt.filtFac = f.binFac;                                                   %Filtering factor
% Initialize variables
k = 1;                                                                     %beamPath counter
sw.G17 = 1;                                                                %G17 plane selected
curMode = NaN;                                                             %Current mode
pos.new = [0,0,0];                                                         %newPosition
pos.interp = [];                                                           %Interpolated positions
pos.dAcc = 0;                                                              %Accumulated distance
% ******************************** READ OUT *******************************
fprintf('\n   -> Opening gCode file ''%s''',f.Name);                          %Screen Output
rawFile = readGfile(fName);                                                %Read in g-code file
f.rawFlength = length(rawFile);                                            %Length of raw file
pos.BP = NaN(opt.preAllocFac*f.rawFlength,3);                              %Allocate memory for beamPath
bb = ones(opt.preAllocFac*f.rawFlength,1);                                 %Initialize beam blanking flag array
Out.plotCnt = ceil(f.rawFlength/Out.plotCnt);                              %Plot Counter
% ********** DETERMINE LINEAR RESOLUTION FROM beam DIAMETER ***************
tmp.dLine = min(find(any(contains(rawFile,'D='),2)));                      %Find line with beam diameter
tmp.l2 = [rawFile{tmp.dLine,:}];                                           %Get lind with beam diameter
tmp.Index = strfind(tmp.l2, 'D=');                                         %Find position of 'D='
if isempty(tmp.Index)
    tmp.dbeam = str2double(inputdlg('Enter beam width [nm]:','Input',1));
else
    tmp.dbeam = sscanf(tmp.l2(tmp.Index(1) + length('D='):end), '%g', 1);  %Extract diameter
end
opt.dRes = (100-f.overlap)/100*tmp.dbeam;                                  %Get linear resolution
assert(isa(opt.dRes,'numeric')&& opt.dRes>0,...
       'Invalid interpolation resolution!');                               %Error checking
fprintf('\n   -> The tool diameter was identified as %.1f [input unit]',tmp.dbeam);     %Screen Output
fprintf('\n   -> The input unit is interpreted as ''%s''',f.unitLabel);       %Screen Output
fprintf('\n   -> The interpolation resolution is %.1f %s',opt.dRes,f.unitLabel);%Screen Output
clear tmp
% ********** Live Output **************************************************
if Out.plot; [h] = plotMaps([],[],[],f,'IniLive',Out); end                 %Initialize Plot
fprintf('\n*** Interpretation and Interpolation\n');
textprogressbar('   -> Interpreting and interpolating gCodes: ','Ini');    %Initialize text progress bar
prgInd = round(f.rawFlength/200);                                          %Progress index
% ******************************** Interpret g-code ***********************
for row = 1:f.rawFlength
    if ~mod(row,prgInd); textprogressbar(row/f.rawFlength*100,'Update');end%Update progress bar
    fullLine = rawFile(row,:);                                             %Get line
    fullLine = fullLine(~cellfun(@isempty,fullLine));                      %Remove empty line fragments
    %Remove potential line numbers
    if strcmp(fullLine{1}(1),'N')
         fullLine = fullLine(2:end);
    end
    %Filter out G18 G19 sections
    if any(contains(fullLine,{'G18','G19'}))
        sw.G17 = 0;
    elseif any(contains(fullLine,'G17'))
        sw.G17 = 1;
    end
    %Ignore command G28 - Return to origin
    if any(contains(fullLine,'G28'))
        continue
    end
    pos.dArc = [0,0,0];                                                    %Reset arcOffsets
    for i = 1:length(fullLine)                                             %Loop over line fragments
        if Out.scrPrint; disp(fullLine{i}); end                            %Screen print
        switch fullLine{i}                                                 %Check for commands G0 - G3
            case 'G0' %Rapid Positioning
                if Out.scrPrint; disp('Rapid positioning'); end
                curMode = 'G0';
            case 'G1' %Linear Interpolation
                if Out.scrPrint; disp('Linear interpolation'); end
                curMode = 'G1';
            case 'G2' %Controlled Arc Move, clockwise
                if Out.scrPrint; disp('Controlled Arc Move, clockwise'); end
                curMode = 'G2';
                if ~sw.G17; curMode = ''; end
            case 'G3' %Controlled Arc Move, counterclockwise
                if Out.scrPrint; disp('Controlled Arc Move, counterclockwise'); end
                curMode = 'G3';
                if ~sw.G17; curMode = ''; end
            otherwise
                switch fullLine{i}(1)                                      %Check for coordinates X, Y, Z, I, J
                    case 'X'
                        pos.new(1) = str2double(fullLine{i}(2:end));
                    case 'Y'
                        pos.new(2) = str2double(fullLine{i}(2:end));
                    case 'Z'
                        pos.new(3) = str2double(fullLine{i}(2:end));
                    case 'I'
                        pos.dArc(1) = str2double(fullLine{i}(2:end));
                    case 'J'
                        pos.dArc(2) = str2double(fullLine{i}(2:end));
                end
         end
    end
    % ******************************** Interpolate coordinates ****************
    if any(strcmp(curMode,{'G0','G1'})) %Lines
        [pos,bb] = interpLine(curMode,pos,opt,bb,k);
    elseif any(strcmp(curMode,{'G2','G3'})) %Arcs
        [pos,bb] = interpArc(curMode,pos,opt,bb,k);
    end
    % Update variables and plots ******************************************
    pos.BP(k:k+size(pos.interp,1)-1,:) = pos.interp;                       %Update pos.BP
    k = k + size(pos.interp,1);                                            %Update counter
    if (Out.plot && ~mod(row,Out.plotCnt))
        set(h.plt,'XData',pos.BP(find(~isnan(pos.BP(:,1))),1),'YData',pos.BP(find(~isnan(pos.BP(:,2))),2));
        drawnow;
    end
    % Eventually Reallocate memory ****************************************
    if k > 0.5*length(pos.BP)
       pos.BP = [pos.BP;NaN(ceil(0.5*opt.preAllocFac*f.rawFlength),3)];    %Allocate memory for beamPath
       bb = [bb;ones(ceil(0.5*opt.preAllocFac*f.rawFlength),1)];           %Initialize beam blanking flag array
    end
end
% *** Output **************************************************************
set(h.plt,'XData',pos.BP(find(~isnan(pos.BP(:,1))),1),'YData',pos.BP(find(~isnan(pos.BP(:,2))),2));
drawnow;
textprogressbar(100,'Update'); textprogressbar('','Deini');                %Finalizing progress bar
aNind = ~all(isnan(pos.BP),2);                                             %Find NaN rows
pos.BP = pos.BP(aNind,:);                                                  %Delete all NaN rows
bb = bb(aNind);                                                            %Delete all NaN rows
end
%% Interpolate Line G-codes
function [pos,bb] = interpLine(Mode,pos,opt,bb,k)
if k > 1
    d = norm((pos.new(1:2) - pos.BP(k-1,1:2)));                                   %Determine distance of linear movement
else
    d = nan;
end
%G0 - Rapid positioning ***************************************************
if strcmp(Mode,'G0')
    if isempty(opt.dcrit) || d > opt.dcrit ||  k ==1                       %G0 - Beam-blanked movement
        pos.interp = pos.new;                                              %Update position
        bb(k) = 0;                                                    %Blank beam
    else                                                                   %G1 - Do standard linear interpolation
        Mode = 'G1';
    end
end
%G0 - Linear interpolation ************************************************
if strcmp(Mode,'G1')
    if d > 0 % Check non-zero distance
        direction = (pos.new - pos.BP(k-1,:))/d;                           %Compute direction
        pos.interp = pos.BP(k-1,:) + direction.*(0:opt.dRes:d)';           %Interpolate points
        pos.interp = [pos.interp;pos.new];                                 %Append final position
    end
end
end
%% Interpolate Arc G-codes
function [pos,bb] = interpArc(Mode,pos,opt,bb,k)
% Get Arc center positions
if opt.absIJ
    cntrPos = pos.dArc;                                                    %Absolute center position of Arc I J
else
    cntrPos = pos.BP(k-1,:) + pos.dArc;                                    %Relative center position of Arc I J
end
% Get vectors
v(1,:) = (pos.BP(k-1,1:2)-cntrPos(1:2));                                   %First vector
v(2,:) = (pos.new(1:2)-cntrPos(1:2));                                      %Second vector
% Get angles
ang(1) = atan2d(v(1,2),v(1,1));                                            %Start Angle
ang(2) = atan2d(v(2,2),v(2,1));                                            %End Angle
% G2 - Clockwise circle interpolation *************************************
if strcmp(Mode,'G2')
    r = norm(pos.new(1:2)-cntrPos(1:2));
    if ang(2) > ang(1)
       ang(2) = ang(2)-360;
    end
end
% G3 - Counterclockwise circle interpolation ******************************
if strcmp(Mode,'G3')
    r = norm(v(2,:));                                                      %Radius
    ang(find([norm(v(1,:)),norm(v(2,:))] <0.1)) = 0;                       %Set angle to 0 for short vectors
    if ang(2) < ang(1)
        ang(2) = ang(2)+360;
    end
end
d = abs((ang(2)-ang(1))*pi/180*r);                                         %Arc length
%Interpolate positions
nrPts = round(d/opt.dRes)-1;                                               %Number of points on arc
if nrPts > 0
    angles = linspace(ang(1),ang(2),nrPts+2);                              %Vector with angles
    angles = angles(2:end-1);
    pos.interp = [cntrPos(1:2) + [cosd(angles)',sind(angles)']*r,...
                  linspace(cntrPos(3),pos.new(3),length(angles))'];        %Interpolate points
    pos.interp = [pos.interp;pos.new];                                     %Append end position
else
    pos.interp = pos.new;                                                  %Only take over end position
end
end
%% Interpolate beamPath
function [BPinterp,BBinterp] = interpBeamPath(BP,bb,f,opt)
%Additional interpolation of beamPath
%Even the beamPath is already interpolated from G-code, local densification of
%points can occur when the interpolated distance of individual G-codes is
%shorter than the interpolation distance (always at least 1 start and end point
%are determined)
%The present function reinterpolates the beampath to ensure a homogeneous
%point distribution. Further, by specifying f.binFac > 1, binning of the beamPath
%can be achieved.
%**************************************************************************
fprintf('\n*** Interpolation of entire beam path');
% Check binning factor
if f.binFac < 1
    BPinterp = BP;
    BBinterp = bb;
    fprintf('\n   -> No interpolation of beam path applied ***');
    return
else
    fprintf('\n   -> Interpolating beam path with %.1f %s',...
             opt.dRes*f.binFac,f.unitLabel);                               %Screen Output
end
% Initialize
bbInd = find(bb==0);                                                       %Indeces of beamPath for beam blanking
BP = [BP;BP(end,:)];                                                       %Duplicate last position
bbInd = [bbInd;size(BP,1)];                                                %Add last position to beam blanking list
nrbb = size(bbInd,1);                                                      %Nr of beam blanking operations
ind(1) = 1;                                                                %Initialization of index variable
BPinterp = [];                                                             %Interpolated beamPath
BBinterp = [];                                                             %Associated beam blanking array
nrPts.in = size(BP,1);                                                     %Nr of input points
% *** Stepwise linear interpolation ***************************************
for i = 1:nrbb %Loop over block separated by beam blanking
   if bbInd(i) == 1 %Check beam blanking at first coordinate
      ind(i) = 1;
      dBP{1} = BP(1,:);
      rng{i} = 1;
      continue
   end
   ind(i) = bbInd(i);                                                      %Beam blanking indice
   rng{i}= ind(i-1)+1:ind(i)-1;                                            %Range of points between two beam blanking operations
   if size(rng{i},2) <= 1                                                  %No interpolation for single point
       dBP{i} = BP(ind,:);                                                 %Set beam path
       continue
   end
   BPtmp = BP(rng{i},:);                                                   %Choose section of beamPath
   if ~any(any(diff(BPtmp)))                                               %Check whether distance is > 0
       dBP{i} = BP(ind,:);                                                 %Set beam path
       continue
   end
   dChord = sqrt(sum(diff(BPtmp,[],1).^2,2));                              %Chord lengths
   dCumArcAbs = [0;cumsum(dChord)];                                        %Cummulated absolute arc length
   dChord = dChord/sum(dChord);                                            %Normalized chordlength
   dCumArc = [0;cumsum(dChord)];                                           %Cummulated arc length
   nrPts.out = ceil(max(dCumArcAbs)/opt.dRes/f.binFac);                    %Nr of output points
   if nrPts.out > size(rng{i},2)                                           %Intervene if Nrbins > NrPoints
       nrPts.out = size(rng{i},2);                                         %Correct NrOfBins
       warning(['Number of points lower than number of bins.',...
                'Number of bins was reduced']);                            %Issue warning
   end
   dCumArcBin = linspace(0,1,nrPts.out)';                                  %Normalized vector of output points
   [~,binInd] = histc(dCumArcBin,dCumArc);                                 %Bin data
   binInd((binInd <= 0) | (dCumArcBin <= 0)) = 1;                          %Fetch errors
   binInd((dCumArcBin >= 1)) = size(rng{i},2) - 1;                         %Fetch errors
   binInd(binInd >= size(rng{i},2)) = size(rng{i},2) - 1;                  %Fetch errors
   s = (dCumArcBin - dCumArc(binInd))./dChord(binInd);                     %Interpolate
   if isnan(s(1));s(1) = 0; end
   dBP{i} = BPtmp(binInd,:) + (BPtmp(binInd+1,:) - ...
            BPtmp(binInd,:)).*repmat(s,1,size(BP,2));                      %Save incremental beamPath
end
% *** Assembly of interpolated beamPath and BB arrays *********************
for i = 1:nrbb
   BPinterp = [BPinterp;dBP{i}];                                           %Write beamPath section
   BPinterp(end+1,:) = BP(bbInd(i),:);                                     %Add beam blanking position
   BBinterp = [BBinterp,ones(1,size(dBP{i},1))];                           %Write beam blaking array
   BBinterp(end+1) = 0;                                                    %Add beam blanking position
end
BBinterp = BBinterp';                                                      %Flip array
% *** Output
nrPts.rmvd = 100*(1-size(BPinterp,1)/nrPts.in);                            %Fraction of removed points [%]
fprintf('\n   -> Applying a binning factor of %0.1f, %0.0f%% of %0.0f points were removed',...
          f.binFac,nrPts.rmvd,nrPts.in);                                   %Screen output
end
%% Scale beamPath
function [pos,str] = scaleBeamPath(pos,PIA,f,str)
fprintf('\n*** Transforming and Scaling beam path');
% *** Scaling
pos.BPinterp = pos.BPinterp.*PIA.scalFac;                                  %Apply linear scaling
patternSz = (max(pos.BPinterp)-min(pos.BPinterp));                         %Get pattern size in output unit
maxMag = min(min(PIA.max)./(patternSz*f.unitOut*PIA.cal));                 %Find maximum magnification for fitting pattern on screen
str.patMag = PIA.availMags(max(find(maxMag*PIA.relScale > PIA.availMags)));%Find suitable magnification for pattering
pos.BPDAC = round(pos.BPinterp.*f.unitOut.*PIA.cal.*str.patMag);           %Convert beamPath to bits
% *** Rotation
PIA.rot = 360 - PIA.rot;                                                   %Apply clockwise rotation convention
R=[cosd(PIA.rot) -sind(PIA.rot); sind(PIA.rot) cosd(PIA.rot)];             %Create rotation matrix
pos.BPDAC = pos.BPDAC*R';                                                  %Rotate beamPath
if PIA.rot ~= 360
    fprintf('\n   -> The pattern was rotated clockwise by %0.0f�',...
            360 - PIA.rot);                                                %Screen output
end
% *** Mirroring
PIA.flip(1) = abs(PIA.flip(1)-1);                                          %Adapting convention to ensure no mirroring for PIA.flip = [0 0]
if PIA.flip(1) %X
    pos.BPDAC(:,2) = -1*(pos.BPDAC(:,2)-max(pos.BPDAC(:,2)));
else
    fprintf('\n   -> The pattern was mirrored at its horizontal axis');    %Screen output
end
if PIA.flip(2) %Y
   pos.BPDAC(:,1) = -1*(pos.BPDAC(:,1)-max(pos.BPDAC(:,1)));
   fprintf('\n   -> The pattern was mirrored at its vertical axis');       %Screen output
end
% *** Translation
pos.BPDAC = pos.BPDAC - min(pos.BPDAC);                                    %Anquer at 0
if PIA.centrePattern
    pos.BPDAC = pos.BPDAC + 0.5*(PIA.max-max(pos.BPDAC));                  %Centre
end
pos.BPDAC = pos.BPDAC + PIA.offs;                                          %Add offset to PIA scanning range
pos.DACsz = max(pos.BPDAC)-min(pos.BPDAC);                                 %Scanning range
fprintf(['\n   -> The pattern size of %.1f x %.1f %s^2 was transformed to '...
         '%.0f x %.0f patterning points'],patternSz(1),patternSz(2),...
         f.unitLabel,pos.DACsz(1),pos.DACsz(2));                           %Screen output
end
%% Plot Maps
function h = plotMaps(h,pos,bb,f,varargin)
% *** Ini Live Plot
if ~isempty(varargin) && strcmp(varargin{1},'IniLive')
    Out = varargin{2};
    h.fig = figure('units','normalized','outerposition',[0.1 0.1 0.4 0.4]);%Create figure
    h.ax(1) = axes;                                                        %Create axes
    title(h.ax(1),'Live Interpolation');
    set(h.ax(1),'Parent',h.fig,'Color','k');                               %Create axes
    hold on;
    if Out.scatter
        h.plt = scatter(h.ax(1),NaN,NaN,1,'w');
    else
        h.plt = plot(h.ax(1),NaN,NaN,'w');
    end
    daspect([1 1 1]);
    xlabel(h.ax(1),['x [',f.unitLabel,']']);
    ylabel(h.ax(1),['y [',f.unitLabel,']']);
    return
end
% *** Plot PostProcessing
h.titles = {'LivePlot','MachiningStrategy','DensityMap-BeforeInterpolation','DensityMap-AfterInterpolation','ScaledAndTransformedPattern'};
nrPlts = length(h.titles);
% Initialize figures
for i = 2:nrPlts
    h.fig(i) = figure;                                                     %Create figure
    h.ax(i) = axes;                                                        %Create axes
    daspect(h.ax(i),[1 1 1]);                                              %1:1:1 Aspect ratio
    hold(h.ax(i),'on');
    title(h.ax(i),h.titles{i});                                            %Set title
    set(h.ax(i),'Color','k');                                              %Set color and y axis
    axis(h.ax(i),'tight');
end
%Scanning strategy
k = size(pos.BPinterp,1);
Xgr = [pos.BPinterp(:,1),pos.BPinterp(:,1)];
Ygr = [pos.BPinterp(:,2),pos.BPinterp(:,2)];
Zgr = [zeros(k,1),zeros(k,1)];
C = [1:k;1:k]';
Xgr(bb.interp==0,:) = nan; Ygr(bb.interp==0,:) = nan; Zgr(bb.interp==0) = nan;
surface(h.ax(2),Xgr,Ygr,Zgr,C,'facecol','no','edgecol','interp','linew',2);
colormap(h.ax(2),'jet');
h.c(1) = colorbar(h.ax(2));
ylabel(h.c(1),'Pattering order [1]');
xlabel(h.ax(2),['x [',f.unitLabel,']']);
ylabel(h.ax(2),['y [',f.unitLabel,']']);
bbInd.u = find(bb.interp==0);
bbInd.l = bbInd.u-1;
bbInd.l(bbInd.l==0) = 1;
for i = 1:length(bbInd.u)
    plot(h.ax(2),[pos.BPinterp(bbInd.l(i,:),1),pos.BPinterp(bbInd.u(i,:),1)],...
                 [pos.BPinterp(bbInd.l(i,:),2),pos.BPinterp(bbInd.u(i,:),2)],...
                 ':w','linewidth',2);
end
if exist('hist3')
    %DensityPlot uncleaned
    tmpFig = figure('visible','off');
    xBin = ceil(sqrt(ceil(size(pos.BPinterp(:,1),1)))/2);
    %xBin = 100;
    yBin = xBin;
    n = hist3([pos.BP(:,1),pos.BP(:,2)],[xBin yBin],'CDataMode','auto','EdgeColor','none');
    close(tmpFig);
    n1 = n';
    n1(size(n,1)+1,size(n,2)+1)=0;
    xb = linspace(min(pos.BP(:,1)),max(pos.BP(:,1)),size(n,1)+1);
    yb = linspace(min(pos.BP(:,2)),max(pos.BP(:,2)),size(n,1)+1);
    h.pcolor = pcolor(h.ax(3),xb,yb,n1);
    set(h.pcolor,'edgecolor','none');
    colormap(h.ax(3),'jet');
    h.c(2) = colorbar(h.ax(3));
    ylabel(h.c(2),'Nr. of Points [1]');
    xlabel(h.ax(3),['x [',f.unitLabel,']']); ylabel(h.ax(3),['y [',f.unitLabel,']']);
    %DensityPlot cleaned
    tmpFig = figure('visible','off');
    n = hist3([pos.BPinterp(:,1),pos.BPinterp(:,2)],[xBin yBin],'CDataMode','auto','EdgeColor','none');
    close(tmpFig);
    n1 = n';
    n1(size(n,1)+1,size(n,2)+1)=0;
    xb = linspace(min(pos.BPinterp(:,1)),max(pos.BPinterp(:,1)),size(n,1)+1);
    yb = linspace(min(pos.BPinterp(:,2)),max(pos.BPinterp(:,2)),size(n,1)+1);
    h.pcolor = pcolor(h.ax(4),xb,yb,n1);
    set(h.pcolor,'edgecolor','none');
    colormap(h.ax(4),'jet');
    h.c(3) = colorbar(h.ax(4));
    ylabel(h.c(3),'Nr. of Points [1]');
    caxis(h.ax(4),h.c(2).Limits);
    xlabel(h.ax(4),['x [',f.unitLabel,']']); ylabel(h.ax(4),['y [',f.unitLabel,']']);
else
    warning('\nCould not plot density maps with command ''hist3'' - No Statistics toolbox installed.\n');
end
%Scaled and Transformed Pattern
plot(h.ax(5),pos.BPDAC(:,1),pos.BPDAC(:,2),'w');
set(h.ax(5),'Ydir','reverse');
xlabel('DAC-pts X [1]'); ylabel('DAC-pts Y [1]');
%Save images
saveImgs(h,f);
end
%% Save images
function saveImgs(h,f)
%function saveImgs(h)
if ~isdir([fileparts(mfilename('fullpath')),'\Output_StreamFiles\',f.Name])
    mkdir([fileparts(mfilename('fullpath')),'\Output_StreamFiles\',f.Name]);
end
for i = 1:length(h.fig)
    set(h.fig(i),'renderer','opengl','invertHardcopy','off',...
                     'units','inch','outerposition',[1,1,6,6],'color','w');
    title(h.ax(i),'');
    %set(hax_new, 'Position', get(0, 'DefaultAxesPosition'));
    %title(hax_new,'');
    print(h.fig(i),[fileparts(mfilename('fullpath')),'\Output_StreamFiles\',...
                    f.Name,'\FIG',num2str(i),'_',h.titles{i},'.tiff'],'-dtiff','-r300');
    title(h.ax(i),h.titles{i});
end
end
%% ConvertUnit - Determine conversion factor
function unitOut = convertUnit(unitLabel)
%function unitOut = convertUnit(unitLabel)
%Conversion factor from inputunit to m
switch unitLabel
    case 'm'
        unitOut = 1;
    case 'mm'
        unitOut = 1e-3;
    case '�m'
        unitOut = 1e-6;
    case 'nm'
        unitOut = 1e-9;
    case '�'
        unitOut = 1e-10;
end
end
%% Write stream file
function fNameFull = writeStreamFile(pos,bb,settings,fInfo,t)
%function strFilePath = writeStreamFile(pos.BP,settings,fileInfo)
%*** Prepare matrix format ************************************************
nrLines = size(pos.BPDAC,1);                                               %Nr of coordinate lines
dwTime = round(t/nrLines*1e7,-1);                                          %Calculate dwell time
pos.BPDAC = [repmat(dwTime,nrLines,1),pos.BPDAC,bb.interp];                %Add dwell time and BeamBlank flag
pos.BPDAC(end+1,:) = [pos.BPDAC(end,1:end-1),0];                           %Add beam blank line to the end of file
%*** Write header *********************************************************
if ~isdir([fileparts(mfilename('fullpath')),'\Output_StreamFiles\',fInfo.Name])
    mkdir([fileparts(mfilename('fullpath')),'\Output_StreamFiles\',fInfo.Name]);
end
fNameFull = [fileparts(mfilename('fullpath')),'\Output_StreamFiles\',fInfo.Name,'\',...
              num2str(settings.patMag),'xMag_',num2str(t),'s.',settings.fType];
fID = fopen(fNameFull,'w');
fprintf(fID, '%s\n%s\n%s\n', settings.DAC, num2str(settings.nrIter),...
        num2str(nrLines+1));                                               % Write header
fclose(fID);
dlmwrite(fNameFull,pos.BPDAC,'delimiter','\t','precision','%1.0f',...
         'newline','unix','-append');
end
%% Text Progress bar
function textprogressbar(c,mode)
% This function creates a text progress bar. It should be called with a
% STRING argument to initialize and terminate. Otherwise the number correspoding
% to progress in % should be supplied.
% INPUTS:   C   Either: Text string to initialize or terminate
%                       Percentage number to show progress
% OUBPUTS:  N/A
% Example:  Please refer to demo_textprogressbar.m

% Author: Paul Proteus (e-mail: proteus.paul (at) yahoo (dot) com)
% Version: 1.0
% Changes tracker:  29.06.2010  - First version

% Inspired by: http://blogs.mathworks.com/loren/2007/08/01/monitoring-progress-of-a-calculation/

% Initialization
persistent strCR;           %   Carriage return pesistent variable

% Vizualization parameters
strPercentageLength = 10;   %   Length of percentage string (must be >5)
strDotsMaximum      = 10;   %   The total number of dots in a progress bar

% Main
if strcmp(mode,'Ini')
    % Progress bar - initialization
    fprintf('%s',c);
    strCR = -1;
elseif strcmp(mode,'Deini')
    % Progress bar  - termination
    strCR = [];
    fprintf([c '\n']);
elseif strcmp(mode,'Update')
    % Progress bar - normal progress
    c = floor(c);
    percentageOut = [num2str(c) '%%'];
    percentageOut = [percentageOut repmat(' ',1,strPercentageLength-length(percentageOut)-1)];
    nDots = floor(c/100*strDotsMaximum);
    dotOut = ['[' repmat('.',1,nDots) repmat(' ',1,strDotsMaximum-nDots) ']'];
    strOut = [percentageOut dotOut];

    % Print it on the screen
    if strCR == -1,
        % Don't do carriage return during first run
        fprintf(strOut);
    else
        % Do it during all the other runs
        fprintf([strCR strOut]);
    end

    % Update carriage return
    strCR = repmat('\b',1,length(strOut)-1);

else
    % Any other unexpected input
    error('Unsupported argument type');
end
end
%% tileFigs - Tile figures accross screen
function tileFigs()
%function tileFigs()
%Tile all figures evenly spread accros the screen
%% Initialization
mon = 1;                                                                   %Choose monitor number
offset.l = 70; offset.r = 0; offset.b = 70; offset.t = 0;                  %Offsets left right botton top (possible taskbars)
grid = [2 2 2 2 2 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4;
3 3 3 3 3 3 3 3 4 4 4 5 5 5 5 5 5 5 5 6]';                                 %Define figure grid
%% Find figures and screen dimension
h.figs = flip(findobj('type','figure'));                                   %Get figure handles
h.figs = h.figs(find(strcmp(get(h.figs,'visible'),'on')));                 %Only treat visible figures
set(h.figs,'unit','pixels');                                               %Set figure units to [pxs]
nFigs = size(h.figs,1);                                                    %Get number of visible figures
scr.Sz = get(0,'MonitorPositions');                                        %Get screen size
scr.h = scr.Sz(mon,4)-offset.t;                                            %Get screen height
scr.w = scr.Sz(mon,3)-offset.l-offset.r;                                   %Get screen width
scr.orX = scr.Sz(mon,1)+offset.l;                                          %Get screen origin X
scr.orY = scr.Sz(mon,2);                                                   %Get screen origin Y
%% Check limits
if ~nFigs; error('figures are not found'); return; end                     %Stop for no figures
if nFigs > 20; error('too many figures(maximum = 20)'); return; end        %Check for limit of 20 figures
%% Define grid according to screen aspect ratio
if scr.w > scr.h %Widescreen
    n.h = grid(nFigs,1);                                                   %Define number of figures in height
    n.w = grid(nFigs,2);                                                   %Define number of figures in width
else
    n.h = grid(nFigs,2);                                                   %Define number of figures in height
    n.w = grid(nFigs,1);                                                   %Define number of figures in width
end
%% Determine height and width for each figure
fig.h = (scr.h-offset.b)/n.h;                                              %Figure height
fig.w =  scr.w/n.w;                                                        %Figure width
%% Resize figures
k = 1;                                                                     %Initialize figure counter
for i =1:n.h %Loop over height
    for j = 1:n.w  %Loop over width
        if k > nFigs; return; end                                          %Stop when all figures have been resized
        fig_pos = [scr.orX + fig.w*(j-1) scr.h-fig.h*i fig.w fig.h];   %Compute new figure position
        set(h.figs(k),'OuterPosition',fig_pos);                            %Set new figure position
        k = k + 1;                                                         %Increase figure counter
    end
end
end
