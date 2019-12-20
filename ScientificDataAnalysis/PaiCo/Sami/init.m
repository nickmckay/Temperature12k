function init
% Adds all the subdirectories to path
filename = mfilename('fullpath');
[pathstr, name, ext] = fileparts(filename);
addpath(genpath(pathstr));
javaaddpath('common/ParforProgressMonitor.jar');