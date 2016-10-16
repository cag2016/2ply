function CalciumSignalTrace3 

%% Data to work with
% 1. csv file from fiji of Channel 2 (GFP) with Frame #, ROI1, ROI2, Background, 

% 2. csv file from fiji of Channel 3 (Voltage signal for LED from Teensy)
%   with Frame # and Average of each frame = this is used to coordinate the
%   timing of LED step changes with the imaging data.

% 3.From the metadata in the header of each Tif of interest, extract
%   Frame rate and trigger start time 

% 4. Step values  and Step duration from Step log file for each acquisition.

%% Flow chart 

% Import metadata of tif of interest into matlab (rename metadata)
% Extract as variable: Trigger start time from metadata
% Extract as varaible: Frame rate
metadata = scim_openTif()
frameRate = metadata.acq.frameRate;
triggerStartTime = metadata.internal.triggerTimeString;
PreStimtime=60;

% Based on Trigger start time, search and read in  matching log file of steps.
triggerDateVec = datevec( triggerStartTime, 'mm/dd/yyyy HH:MM:SS.FFF' );
logFile = datestr( triggerDateVec, 'yyyy-mm-dd_HH-MM-SS.txt' );
fprintf ('opening logfile = %s\n', logFile);
fidLog = fopen( logFile, 'r' );
patternSteps = 'DutyCycleSteps = [';
patternPrestimTime = 'PrestimTime = ';
patternStepduration = 'Stepduration = ';

while 1
    tLine = fgetl( fidLog );
    if( ~ischar(tLine) )
        break
    end

    if( strfind( tLine, patternSteps ))
        eval( [ tLine ';' ]);
    end
    
     if( strfind( tLine, patternPrestimTime ))
        eval( [ tLine ';' ]);
     end
    
      if( strfind( tLine, patternStepduration ))
        eval( [ tLine ';' ]);
    end
end
fclose( fidLog );

% Load the Channel 3 csv.
% Take note of frames that differ significantly from the two neighboring
%   frames in mean pixel value. These are bad frames which catch the Step changes. 


[ filename, pathname ] = uigetfile( '*', 'Pick Channel 3 CSV' );
Channel3Table = dlmread( [ pathname filename ], '\t', 1, 0 );
Voltagemean= Channel3Table (:,2);

startFrame = ceil(frameRate * 60);

StepTransitionFrame = startFrame:4:length(Voltagemean);

%perfect transition frames manually put in
indx = find( StepTransitionFrame >= 377 );
StepTransitionFrame(indx) = StepTransitionFrame(indx) - 1;

indx = find( StepTransitionFrame >= 1732 );
StepTransitionFrame(indx) = StepTransitionFrame(indx) - 1;



StepValuePerFrame = zeros( size( Voltagemean ));
indxStep = 1;
for i = 1 : 1 : length( StepTransitionFrame )

    if( mod( i, 100 ) == 0 )
        fprintf('%i\n', i );
    end
    
    cIndex = StepTransitionFrame(i);
    if( i == 1 )
        StepValuePerFrame(1:cIndex) = 0;
    else
        cIndexLast = StepTransitionFrame(i-1);
        
        if( indxStep > length( DutyCycleSteps ))
            StepValuePerFrame(cIndexLast:cIndex) = 0;
        else
            StepValuePerFrame(cIndexLast:cIndex) = DutyCycleSteps(indxStep);
            indxStep = indxStep+1;
        end
    end
    
end



    


% StepTransitionFrame = zeros(size(Voltagemean));
% 
%  for i = 2 : 1 : length(Voltagemean)-1
%      if( Voltagemean(i)~=Voltagemean(i+1) && Voltagemean(i)~=Voltagemean(i-1) )
%         StepTransitionFrame(i)= 1;
%   
%      end
%    
%  end
 
 % Incorporate step values and associate every frame with a step value (bad frames not discarded yet). 
     
% StepValueperFrame=zeros(size(Voltagemean));
% indx = 0;
% 
% for i=1:1:length(StepValueperFrame)
%  
%     if (StepTransitionFrame==1)
%         indx= indx+1;
%         
%     end
%     
%     if (indx>0)
%         StepValueperFrame(i)= DutyCycleSteps(indx)
%     end
%     
%     
%     
% end

% Load Channel 2 csv. 
% Separate ROI1, ROI2, Background 
% Calculate, ROI1-Background and ROI2-Background
% Pick Prestim fluorescence Frames for ROI1 and ROI2 
% Calculate F0 Average Prestim fluorescence for ROI1 and ROI2
% Calculate Normalized Calcium signal for ROI1 and ROI2

[ filename, pathname ] = uigetfile( '*', 'Pick Channel 2 CSV' );
Channel2Table = dlmread( [ pathname filename ], '\t', 1, 0 );

Frameno = Channel2Table(:,1);
ROI1=Channel2Table(:,2);
%ROI2= ResultsTable(:, 3);
BackgroundROI = Channel2Table(:, 3);
PrestimFramesUsed = round(PreStimtime * frameRate)-15;


% Substract background from ROIs
ROI1minusBackground = ROI1-BackgroundROI;
% ROI2minusBackground= ROI2-BackgroundROI;

%Pick Prestim fluorescence Frames for ROI1 and ROI2
PrestimFluoroscenceROI1 = ROI1minusBackground(1:PrestimFramesUsed);
% PrestimFluoroscenceROI2 = ROI2minusBackground(1:PrestimFramesUsed);

%F0 Average Prestim fluorescence for ROI1 and ROI2
F0ROI1 = mean (PrestimFluoroscenceROI1);
% F0ROI2 = mean (PrestimFluoroscenceROI2);

%Normalized Calcium signal for ROI1 and ROI2
NormalizedROI1signal = ROI1minusBackground/F0ROI1;
% NormalizedROI2signal = ROI2minusBackground/F0ROI2;



% Remove bad frames noted from Channel 3 csv. 
% Have 4 columns of final data from which plots will be made.
%    1. Frame number (excluding bad frames)(
%    2. Step value associated with each Frame number.(StepValuePerFrame)
%    3. Normalized calcium signal for ROI1 (NormalizedROI1signal)
%    4. Normalized calcium signal for ROI2

CleanFrames = setdiff( 1:1:length(StepValuePerFrame), StepTransitionFrame );
Time= CleanFrames/frameRate;
CleanNormalizedROI1signal= NormalizedROI1signal(CleanFrames);
CleanStepValuePerFrame=StepValuePerFrame(CleanFrames);

%Plot normalized calcium signals for ROI1 and ROI2
% Coumn 3 and 4 are the values for the y axis (normalized calcium signal for ROI1 and ROI2)
% The Time axis will computed from the Step value (column 2)s and the 
    %Step Duration (variable extracted from log file)
    %For example, if the step duration is 1 second, then a consecutive series of identical step
   %values in columm 2 will be equal to 1 second and so on. 
   figure;
   
   plot( Time, CleanNormalizedROI1signal, 'bo' );
   
%find peaks
minPeakProminence = 1;
minPeakDistance = 5;

[eventValues, eventIndices, eventWidth, eventProminence] = ... 
    findpeaks(CleanNormalizedROI1signal, 'MinPeakProminence', minPeakProminence, 'MinPeakDistance', minPeakDistance);
    TimeEvent=Time (eventIndices);

 figure; findpeaks(CleanNormalizedROI1signal, 'MinPeakProminence', minPeakProminence, 'MinPeakDistance', minPeakDistance);
 figure; plot( Time, CleanNormalizedROI1signal);
 hold on;
 plot(TimeEvent, eventValues,'ro');
 
 
 %peak triggered stepvalue
 
 nSamplesForAvg = 71;
 avgEventStepValue= zeros( nSamplesForAvg, 1 );
 avgEventSignal = zeros( nSamplesForAvg, 1 );
 avgEventSamples = zeros( nSamplesForAvg, 1 );
 
 
 EventStepValues = zeros(nSamplesForAvg, length(TimeEvent));

 for i = 1:length(TimeEvent)
     
     framesToAvg = ( -((nSamplesForAvg-1)/2) : 1 : ((nSamplesForAvg-1)/2) ) + eventIndices(i);
     
     for iFrame = 1 : 1 : length( framesToAvg )
        
         cFrame = framesToAvg(iFrame);
         
         % If the current frame to be added is one of the transition
         % frames, don't add it
%          if( ismember( cFrame, StepTransitionFrame ))
             
%          else % Otherwise, if it's a clean frame, add it
             avgEventStepValue(iFrame) = avgEventStepValue(iFrame) + StepValuePerFrame(cFrame);
             avgEventSignal(iFrame) = avgEventSignal(iFrame) + NormalizedROI1signal(cFrame);
             avgEventSamples(iFrame) = avgEventSamples(iFrame) + 1;
             
             EventStepValues(iFrame, i) = StepValuePerFrame(cFrame);
%          end
         
     end
     
     
     figure;
     plot( (-((nSamplesForAvg-1)/2) : 1 : ((nSamplesForAvg-1)/2))/frameRate, StepValuePerFrame(framesToAvg) );
     hold on;
     plot( (-((nSamplesForAvg-1)/2) : 1 : ((nSamplesForAvg-1)/2))/frameRate, StepValuePerFrame(framesToAvg), 'b*' );
 end
     
%plot calcium triggered average
figure;stairs( (-((nSamplesForAvg-1)/2) : 1 : ((nSamplesForAvg-1)/2))/frameRate,avgEventStepValue./avgEventSamples)
     
 %errorbar

 figure;errorbar((-((nSamplesForAvg-1)/2) : 1 : ((nSamplesForAvg-1)/2))/frameRate, (avgEventStepValue./avgEventSamples),std(EventStepValues'))

     
     
     
     
%      cEvent= int64(TimeEvent(i));
%      cIndex = (cEvent==StepTime);
%      subplot( 2, 1, 1 );
%      plot(cEvent,eventValues(i), 'bo');
%      subplot( 2, 1, 2 );
% 
%     if( ~isempty( cIndex ))
% 
%      if( cIndex > 2 )
% %          plotyy( cEvent, eventValues(i), StepTime(cIndex-2:cIndex), Steps( cIndex-2:cIndex ));
%          plot(StepTime(cIndex-2:cIndex),Steps(cIndex-2:cIndex), 'r*' );
%      else
% %          plotyy( cEvent, eventValues(i), StepTime(cIndex), Steps(cIndex ));
%          plot(StepTime(cIndex),Steps(cIndex), 'r*' );
%      end
%      %
% %      plotyy( cEvent, eventValues(i), StepTime(cIndex-2:cIndex), Steps(cIndex-2:cIndex) );
% %      keyboard;
%     end
%  end    
 
keyboard;
end

function thisIsOldStuff
%% Variables

Stepoffset=60;
dataStart = 59.9752;
dataEnd = 557.35; % start of last pulse of intense triplet near end
ledStart = 60;
ledEnd = 458; % start of last pulse of intense triplet near end
scaleFactor = ( dataEnd - dataStart ) / ( ledEnd - ledStart );

minPeakProminence = 1;
minPeakDistance = 5;

%% Pick CSV Results table of measured ROI means

[ strFilename, strPathname ] = uigetfile( '*', 'Pick a csv ROI results file' );
strFullPath = fullfile( strPathname, strFilename );
ResultsTable = dlmread(strFilename,'\t',1,0);
load StepsX010.mat;
StepTime = zeros( size( StepsX010 ));
for i = 1 : 1 : length( StepTime )
    StepTime(i) = i-1;
end
StepTime = StepTime + Stepoffset;
%StepTime = (( StepTime - ledStart ) * scaleFactor ) + dataStart;


% Steps = [ zeros( Stepoffset, 1 ); Steps ];
% StepTime = zeros( size( Steps ));
% for i = 1 : 1 : length( StepTime )
%     StepTime(i) = (1/4.06901041666667) + (i-1);
% end

% Stepdelta= [0; diff(Steps)];

indx = find( StepTime > Stepoffset, 1, 'first' );
Steps = StepsX010(indx:end);
StepTime = StepTime(indx:end);

% Steps = StepsX010(1:indx);
% StepTime = StepTime(1:indx);

% Stepdelta = Stepdelta(1:indx);

% Pick ROI1, ROI2 and Background columns
Frameno = ResultsTable(:,1);
Time=Frameno/4.06901041666667;
ROI1=ResultsTable(:,2);
ROI2= ResultsTable(:, 3);
BackgroundROI = ResultsTable(:, 4);
Backdif = abs(ResultsTable(:,end-1)-ResultsTable(:,end));

Goodindex=find(Backdif<6);

Time=Time(Goodindex);
ROI1=ROI1(Goodindex);
ROI2=ROI2(Goodindex);
BackgroundROI= BackgroundROI(Goodindex);

% Substract background from ROIs
ROI1minusBackground = ROI1-BackgroundROI;
ROI2minusBackground= ROI2-BackgroundROI;

%Pick Prestim fluorescence Frames for ROI1 and ROI2
PrestimFluoroscenceROI1 = ROI1minusBackground(1:244);
PrestimFluoroscenceROI2 = ROI2minusBackground(1:244);

%F0 Average Prestim fluorescence for ROI1 and ROI2
F0ROI1 = mean (PrestimFluoroscenceROI1);
F0ROI2 = mean (PrestimFluoroscenceROI2);

%Normalized Calcium signal for ROI1 and ROI2
NormalizedROI1signal = ROI1minusBackground/F0ROI1;
NormalizedROI2signal = ROI2minusBackground/F0ROI2;

%Plot normalized calcium signals for ROI1 and ROI2
figure;
subplot( 2, 1, 1 );
plot( Time, NormalizedROI1signal);
hold on;
%plot( Time, NormalizedROI2signal);
tmp = axis;
axis([ 0 560 tmp(3) tmp(4) ]);

subplot( 2, 1, 2 );
stairs( StepTime, Steps );
% plot( StepTime, Steps, 'b*' );
%plot( Steps );
tmp = axis;
axis([ 0 560 tmp(3) tmp(4) ])

% subplot( 3, 1, 3 );
%  stairs( StepTime, Stepdelta );
%  tmp = axis;
%  axis([ 100 tmp(2) tmp(3) tmp(4) ])

% figure;
% plotyy( Time, NormalizedROI1signal, StepTime, Steps );
% hold on;
% plot( Time, NormalizedROI2signal, 'g' );
% 
% figure;plot(Time, NormalizedROI1signal, 'b--o')
% hold on
% plot(Time, NormalizedROI2signal, 'r--o')
% stairs( StepTime, Steps );

% Finding peaks/events
[eventValues, eventIndices, eventWidth, eventProminence] = ... 
    findpeaks(NormalizedROI1signal, 'MinPeakProminence', minPeakProminence, 'MinPeakDistance', minPeakDistance);
TimeEvent=Time (eventIndices);

 figure; findpeaks(NormalizedROI1signal, 'MinPeakProminence', minPeakProminence, 'MinPeakDistance', minPeakDistance);
 figure; plot( Time, NormalizedROI1signal);
 hold on;
 plot(TimeEvent, eventValues,'ro');
 
figure;
subplot( 2, 1, 1 );
hold on;
subplot( 2, 1, 2 );
hold on;

 
 for i = 1:length(TimeEvent)
     cEvent= int64(TimeEvent(i));
     cIndex = find(cEvent==StepTime);
     subplot( 2, 1, 1 );
     plot(cEvent,eventValues(i), 'bo');
     subplot( 2, 1, 2 );

    if( ~isempty( cIndex ))

     if( cIndex > 2 )
%          plotyy( cEvent, eventValues(i), StepTime(cIndex-2:cIndex), Steps( cIndex-2:cIndex ));
         plot(StepTime(cIndex-2:cIndex),Steps(cIndex-2:cIndex), 'r*' );
     else
%          plotyy( cEvent, eventValues(i), StepTime(cIndex), Steps(cIndex ));
         plot(StepTime(cIndex),Steps(cIndex), 'r*' );
     end
     %
%      plotyy( cEvent, eventValues(i), StepTime(cIndex-2:cIndex), Steps(cIndex-2:cIndex) );
%      keyboard;
    end
 end    

keyboard;
end
