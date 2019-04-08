% add folders to path
addpath('data','functions');
% load data
load('DATA_APC.mat');
load('Elektromodellflug');
load('axi_motor_db.mat');


% create the properties of the classes
[aero,prop,motor,esc,bat,env,mission] = createProperties( ...
    'settings_default.m', DATA_APC,axi_motor_db,Elektromodellflug);


disp('successfully loaded settings');