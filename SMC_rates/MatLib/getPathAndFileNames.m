% This script is for retrieving path and file names from multiple
% directories

function [pathAndFileNames,FileName,PathName] = getPathAndFileNames(fileSelectionOptions,multiselect,suggestedFilePath,dialogue)
if(~exist('multiselect','var'))
    multiselect = 1;
end
if(~exist('suggestedFilePath','var'))
    suggestedFilePath = '';
end
if(~exist('fileSelectionOptions','var'))
    fileSelectionOptions = {'*.*','All Files '};
end
if(~exist('dialogue','var'))
    dialogue = 'Select a file...';
end
response = 'yes';
FileName = cell(0);
PathName= cell(0);
while (strcmp(response,'No')==0 )
    [FileNameTemp,PathNameTemp,FilterIndex] = uigetfile(fileSelectionOptions,'multiselect','on',dialogue,suggestedFilePath);   
    if FilterIndex
        if iscell(FileNameTemp)
            for f=1:size(FileNameTemp,2)
                FileName{end+1,1} = FileNameTemp{f}; %#ok<*AGROW>
                PathName{end+1,1} = PathNameTemp;
            end
        else
            FileName{end+1,1} = FileNameTemp;
            PathName{end+1,1} = PathNameTemp;
        end
    end
    if(multiselect)
        response = questdlg('Keep choosing files?','Select a file','Yes','No','Yes');
        if ~isempty(PathName); 
            suggestedFilePath = PathName{end}; 
        end
    else
        response = 'No';
    end
end
pathAndFileNames = cell(size(FileName,1),1);
for f=1:size(FileName,1)
    pathAndFileNames{f} = strcat(PathName{f},FileName{f});
end
return;