function [mapTarget_,mapError_] = SATDefinitions()
%SATDEFINITIONS Definitions of columns in the following variables from RH's
%               SAT data:
%               Target_
%               Error_
%   Detailed explanation goes here

    %% Target_ variable as column names and descriptions as map:
    % From teba/Users/Rich/Mat_Code/Import/EventTranslator.m Lines:123-186
    %  Comment for each column Index is used as column names for that column
    % 
    colNames = {};
    colDesc = {};
    % Column #1
    colDesc = [colDesc; 'Target align time'];
    colNames =  [colNames;'TargetAlignTime'];      
    % Column #2
    colDesc = [colDesc; 'Target location (255 = catch trial)'];
    colNames =  [colNames;'TargetLocation_255_IsCatchTrial'];      
    % Column #3
    colDesc = [colDesc; 'TARGET COLOR'];
    colNames =  [colNames;'TargetColor'];  
    % Column #4
    colDesc = [colDesc; 'Min Target Color (CLUT value; higher values = BRIGHTER [more white])'];
    colNames =  [colNames;'MinTargetColorCLUT_HighIsBright'];  
    % Column #5
    colDesc = [colDesc; 'Set Size'];
    colNames =  [colNames;'SetSize'];  
    % Column #6
    colDesc = [colDesc; 'SET SIZE CONDITION (0 = 2; 1 = 4; 2 = 8; 3 = random)'];
    colNames =  [colNames;'SetSizeCondition'];  
    % Column #7
    colDesc = [colDesc; 'Easy hard  Target type; must index to LAST ''3007'' to find correct value.  I(RH) do notunderstand this coding scheme...'];
    colNames =  [colNames;'EasyTargetType']; 
    % Column #8
    colDesc = [colDesc; 'Easy hard  Target type; must index to LAST ''3007'' to find correct value.  I(RH) do notunderstand this coding scheme...'];
    colNames =  [colNames;'HardTargetType']; 
    % Column #9
    colDesc = [colDesc; 'Task Type  (0 = fixation; 1 = detection; 2 = search/MG??'];
    colNames =  [colNames;'TaskType']; 
    % Column #10
    colDesc = [colDesc; 'Hold Time (use to determine search or memory guided)'];
    colNames =  [colNames;'HoldTime']; 
    % Column #11
    colDesc = [colDesc; 'Homogeneous (0)/Non-homogeneous (1)'];
    colNames =  [colNames;'HomogeneousOrHeterohenous']; 
    % Column #12
    colDesc = [colDesc; 'Eccentricity'];
    colNames =  [colNames;'Eccentricity']; 

    mapTarget_.IndexToName = containers.Map((1:12)',colNames);
    mapTarget_.IndexToDesc = containers.Map((1:12)',colDesc);


%% Error_ variable as column names and descriptions as map:
    colNames = {};
    colDesc = {};
    % From teba/Users/Rich/Mat_Code/Import/EventTranslator.m Lines:87-95
    %  Comment for each column Index is used as column names for that column
    %find Error Codes
    %Types of Errors
    %Column:
    %1 = CatchError
    colDesc = [colDesc; 'CatchError'];
    colNames =  [colNames;'CatchError'];
    %2 = HoldError
    colDesc = [colDesc; 'HoldError'];
    colNames =  [colNames;'HoldError'];
    %3 = Latency Error
    colDesc = [colDesc; 'Latency Error'];
    colNames =  [colNames;'LatencyError'];
    %4 = Target Hold Error
    colDesc = [colDesc; 'Target Hold Error'];
    colNames =  [colNames;'TargetHoldError'];
    %5 = Saccade Direction Error
    colDesc = [colDesc; 'Saccade Direction Error'];
    colNames =  [colNames;'SaccadeDirectionError'];
       
    mapError_.IndexToName = containers.Map((1:5)',colNames);
    mapError_.IndexToDesc = containers.Map((1:5)',colDesc);
end


