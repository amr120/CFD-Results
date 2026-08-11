name = 'SMALLDF2_avg';
g = ts_read_hdf5([dr.ts name '.hdf5']);

%stator - i_in = 2; i_out = 0;
%rotor - i_in = 9; i_out = 7;


Rc = 0.08; %load in propeller tip radius

val1 = [];

for i = 1:length(g{1}.r)
    % Assuming val1mean is a function that works with a slice of g{1}.r up to index i
    val1 = [val1, mean(g{1}.r(:,i))]; 
end

% Calculate the minimum difference between val1 and Rc
[~, idx] = min(abs(val1 - Rc));

% Take cuts
i_out = 5; i_in = 5; j_max = idx;

%Define Block IDs of inlet and exit to the stator and Rotor
statorin = 2;
statorout = 0;
rotorin = 9;
rotorout=7;



%--------------CALCULATE LOSSES FOR ROTOR--------------

inlet = ts_structured_cut(g,rotorin,i_in,i_in,1,j_max,1,'en');
outlet = ts_structured_cut(g,rotorout,i_out,i_out,1,j_max,1,'en');


% Calculate inlet conditions
V1_rel = ts_mass_average(inlet,'V_rel',3);
P01 = ts_mass_average(inlet,'Po_rel',2);
[s1,mass_1] = ts_mass_average(inlet,'s',2);
p1 = ts_area_average(inlet,'P',2);

% Exit conditions
T2 = ts_mass_average(outlet,'T',3);
[s2,mass_2] = ts_mass_average(outlet,'s',2);
P02 = ts_mass_average(outlet,'Po_rel',2);


% Cumilative mass flow fraction
mass_2 = [0 cumsum(mass_2)];
mass_1 = [0 cumsum(mass_1)];

% Non-dimensionalise
mass_1 = mass_1 / mass_1(end);
mass_2 = mass_2 / mass_2(end);

% Interpolate inlet entropy to the same mass flow fraction as outlet
s1_interp = interp1(mass_1,s1,mass_2);
P01_interp = interp1(mass_1,P01,mass_2);
P02_interp = interp1(mass_1,P02,mass_2);
p1_interp = interp1(mass_1,p1,mass_2);


% Loss coefficient
%Zeta = (T2 * (s2 - s1)) / (0.5 * V1_rel^2);
Zeta = (T2 * (s2 - s1_interp)) / (0.5 * V1_rel^2);
Yp = (P01_interp - P02_interp)./(P01_interp - p1_interp);

figure(); hold on; grid on; box on;
xlabel('Yp'); ylabel('Cumilative mass flow fraction');
plot(Yp,mass_2);
title('Rotor')

writematrix(Yp,[dr.ts name 'ROTORTURBOSTREAMYp.csv'])
writematrix(mass_2,[dr.ts name 'ROTORTURBOSTREAMMass.csv'])

ts_export_paraview(g,[dr.pv name '_rotors.hdf5'],'HighSpeed',0:6,1);



%-----------------------STATORS------------------------------------------------------


inlet = ts_structured_cut(g,statorin,i_in,i_in,1,j_max,1,'en');
outlet = ts_structured_cut(g,statorout,i_out,i_out,1,j_max,1,'en');


% Calculate inlet conditions
V1 = ts_mass_average(inlet,'V',3);
P01 = ts_mass_average(inlet,'Po',2);
[s1,mass_1] = ts_mass_average(inlet,'s',2);
p1 = ts_mass_average(inlet,'P',2);


% Exit conditions
T2 = ts_mass_average(outlet,'T',3);
[s2,mass_2] = ts_mass_average(outlet,'s',2);
P02 = ts_mass_average(outlet,'Po',2);

p2 = ts_mass_average(outlet,'P',2);

% Cumilative mass flow fraction
mass_2 = [0 cumsum(mass_2)];
mass_1 = [0 cumsum(mass_1)];

% Non-dimensionalise
mass_1 = mass_1 / mass_1(end);
mass_2 = mass_2 / mass_2(end);

% Interpolate inlet entropy to the same mass flow fraction as outlet
s1_interp = interp1(mass_1,s1,mass_2);
P01_interp = interp1(mass_1,P01,mass_2);
P02_interp = interp1(mass_1,P02,mass_2);
p2_interp = interp1(mass_1,p2,mass_2);
p1_interp = interp1(mass_1,p1,mass_2);

% Loss coefficient
Zeta = (T2 * (s2 - s1_interp)) / (0.5 * V1^2);
Yp = (P01_interp - P02_interp)./(P01_interp - p1_interp);


figure(); hold on; grid on; box on;
xlabel('Yp'); ylabel('Cumilative mass flow fraction');
plot(Yp,mass_2);
title('Stator')

writematrix(Yp,[dr.ts name 'STATORTURBOSTREAMYp.csv'])
writematrix(mass_2,[dr.ts name 'STATORTURBOSTREAMMass.csv'])

ts_export_paraview(g,[dr.pv name '_stators.hdf5'],'HighSpeed',7:13,1);

%ts_export_paraview(g,[dr.pv 'Feb7Prop' 'HS.hdf5'],'HighSpeed',[],1);