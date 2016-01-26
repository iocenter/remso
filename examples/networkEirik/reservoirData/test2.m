% mrstPath = '/usr/local/MRST/mrst-2013a';
% run(fullfile(mrstPath,'startup.m'))

% change controls
i1 = [300, 252.736, 300, 237.045, 138.37]*meter^3/day;
i2 = [152.222, -2.50224e-09, 73.7928, 13.3054, 86.2628]*meter^3/day;
p1 = [226.08, 232.332, 255.718, 276.489, 286.352]*barsa;
p2 = [227.019, 241.283, 269.806, 350, 350]*barsa;
p3 = [0, 0, 0, 0, 0]*meter^3/day;
p4 = [226.256, 233.6, 255.611, 272.09, 282.53]*barsa;
p5 = [0, 0, 0, 0, 0]*meter^3/day;
u  = [i1; i2; p1; p2; p3; p4; p5];

% try
 %write control file:
 fid = fopen('SIMPLE10x5x10_CONTROLS.TXT', 'w');
 fprintf(fid, '%+12.6e %+12.6e %+12.6e %+12.6e %+12.6e %+12.6e %+12.6e\n',  [i1; i2; p1; p2; p3; p4; p5]);
 fclose(fid);

 runSim2;
% catch err
%  warning('something went wrong with MRST');
%  exit(3)
% end
