function checkCurrentFolder
%% Change the current folder to the folder of this m-file.
if(~isdeployed)
    folder = cd(fileparts(which(mfilename)));
    % Add that folder plus all subfolders to the path.
    addpath(genpath(folder));
end