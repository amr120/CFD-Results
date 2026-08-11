des_name = '22042026Prop';
name = des_name;
%g = ts_read_hdf5([dr.ts des_name 'restart2_avg.hdf5']);
 g = ts_read_hdf5([dr.ts des_name 'restart2_O2_avg.hdf5']);

s = importdata([dr.geom des_name '-mission.json'],'r');
mymission = jsondecode(s{1});

Rc = mymission.Rc; %load in propeller tip radius

val1 = [];

for i = 1:length(g{1}.r)
    % Assuming val1mean is a function that works with a slice of g{1}.r up to index i
    val1 = [val1, mean(g{1}.r(:,i))]; 
end

% Calculate the minimum difference between val1 and Rc
[~, idx] = min(abs(val1 - Rc));

% Take cuts
i_out = 5; i_in = 5; j_max = idx;

rotorin = 2;
rotorout=0;

inlet = ts_structured_cut(g,rotorin,i_in,i_in,1,j_max,1,'en');
outlet = ts_structured_cut(g,rotorout,i_out,i_out,1,j_max,1,'en');

% Calculate inlet conditions
%V1_rel = ts_mass_average(inlet,'V_rel',3);
%s1 = ts_mass_average(inlet,'s',3);

% Exit conditions
%T2 = ts_mass_average(outlet,'T',3);
%[s2,mass] = ts_mass_average(outlet,'s',2);

% Loss coefficient
%Z%eta = (T2 * (s2 - s1)) / (0.5 * V1_rel^2);

% Mass flow fraction
%mass_cumi = [0 cumsum(mass)];

% Non-dimensionalise
%mass_cumi = mass_cumi / mass_cumi(end);

% Plot loss
%figure(); hold on; grid on; box on;
%xlabel('Entropy loss coefficient'); ylabel('Cumilative mass flow fraction');
%plot(Zeta,mass_cumi);

% Export paraview for streamlines
%ts_export_paraview(g,[dr.pv 'test.hdf5'],'HighSpeed',[],1);




% Calculate inlet conditions
V1_rel = ts_mass_average(inlet,'V_rel',3);
P01 = ts_mass_average(inlet,'Po_rel',2);
[s1,mass_1] = ts_mass_average(inlet,'s',2);
p1 = ts_area_average(inlet,'P',2);

% Exit conditions
T2 = ts_mass_average(outlet,'T',3);
[s2,mass_2] = ts_mass_average(outlet,'s',2);
P02 = ts_mass_average(outlet,'Po_rel',2);
p2 = ts_area_average(outlet,'P',2);

M2 = mass_2;
M1 = mass_1;



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
p2_interp = interp1(mass_1,p2,mass_2);


% Loss coefficient
%Zeta = (T2 * (s2 - s1)) / (0.5 * V1_rel^2);
Zeta = (T2 * (s2 - s1_interp)) / (0.5 * V1_rel^2);
Yp = (P01_interp - P02_interp)./(P01_interp - p1_interp);



% Plot loss
figure(); hold on; grid on; box on;
xlabel('DeltaP'); ylabel('Cumilative mass flow fraction');
plot(p2_interp ,mass_2);
title('Rotor')

figure(); hold on; grid on; box on;
xlabel('Yp'); ylabel('Cumilative mass flow fraction');
plot(Yp,mass_2);
plot(mymission.LOSS, (mymission.R - mymission.Rh)./(mymission.Rc - mymission.Rh))
ylim([0,1])
xlim([-0.01, 0.05])
title('Rotor')


fmt = ['The rotor Zeta is: [', repmat('%g, ', 1, numel(Zeta)-1), '%g]\n'];
%fprintf(fmt, Zeta)
fmt = ['The rotor vector YP is: [', repmat('%g, ', 1, numel(Yp)-1), '%g]\n'];
fprintf(fmt, Yp)
fmt = ['The rotor mass-flow-ratio is: [', repmat('%g, ', 1, numel(mass_2)-1), '%g]\n'];
fprintf(fmt, mass_2)

writematrix(Yp,[dr.ts name 'ROTORTURBOSTREAMYp.csv'])
writematrix(mass_2,[dr.ts name 'ROTORTURBOSTREAMMass.csv'])

ts_export_paraview(g,[dr.pv name 'HS.hdf5'],'HighSpeed',[],1);

p1 = ts_area_average(inlet,'P',3);
p2 = ts_area_average(outlet,'P',3);
PressureRatio = p2/p1

%rho1 = ts_mass_average(inlet,'Ro',2);

%mass = ts_area_average(inlet,'Mdot',3);
%ts_export_paraview(g,[dr.pv 'Feb7Prop' 'HS.hdf5'],'HighSpeed',[],1);

leg = {'1' '10' '20' '30'};
plot_thrust = 1;
%% Plot thrust from nacelle cases and compare against target values
if plot_thrust == 1
    
    % Block ids and indices to assemble control volume
    bid_far = 18; bid_jet = 17; 
%     i_in = 7; i_out = 7;
    i_in = 41; i_out = 72;
    
    % Loop over all operating points
    for o = 1:length(leg)
    
        % Read the case
        %g = ts_read_hdf5([dr.ts des_name 'E_O' num2str(o) '_avg.hdf5']);
%         g = ts_read_hdf5([dr.ts des_name 'E_O' num2str(o) '.hdf5']);
    
        % Take cuts around propulsor
        j = g{bid_far+1}.attribute.nj-i_out;
        i = g{bid_jet+1}.attribute.ni-i_out;
        C{1} = ts_structured_cut(g,bid_far,i_in,i_in,1,j,1,'en');
        C{2} = ts_structured_cut(g,bid_far,i_in+1,g{bid_far+1}.attribute.ni-1,j,j,1,'en');
        C{3} = ts_structured_cut(g,bid_jet,i,i,1,'en',1,'en');

        % Join cuts together
        varnames = fieldnames(C{1}); flips = [3 3 1];
        for v = 1:length(varnames)
            if isnumeric(C{1}.(varnames{v})) == 1
                c.(varnames{v}) = [];
                for n = 1:length(C)
                    c.(varnames{v}) = [c.(varnames{v}) ; flip(C{n}.(varnames{v}),flips(n))];
                end
            end
        end
        
        % Record other variables
        c.av = C{1}.av; c.bv = C{1}.bv; c.attribute = C{1}.attribute;
        
        % Calculate diagonal lengths
        c = ts_secondary(c); c.t = c.rt ./ c.r;
        t_av = 0.25 * (c.t(2:end,2:end) + c.t(1:end-1,1:end-1) + ...
            c.t(1:end-1,2:end) + c.t(2:end,1:end-1));
        x1 = c.x(2:end,2:end) - c.x(1:end-1,1:end-1); x2 = c.x(1:end-1,2:end) - c.x(2:end,1:end-1);
        r1 = c.r(2:end,2:end) - c.r(1:end-1,1:end-1); r2 = c.r(1:end-1,2:end) - c.r(2:end,1:end-1);
        rt1 = c.r(2:end,2:end) .* (c.t(2:end,2:end) - t_av)...
            - c.r(1:end-1,1:end-1) .* (c.t(1:end-1,1:end-1) - t_av);
        rt2 = c.r(1:end-1,2:end) .* (c.t(1:end-1,2:end) - t_av)...
            - c.r(2:end,1:end-1) .* (c.t(2:end,1:end-1) - t_av);

        % Correct velocities on walls
        q = c.mwall == 0; c.rovx(q) = 0; c.rovr(q) = 0;

        % Area components
        Ax = 0.5 * (r1.*rt2 - r2.*rt1);
        Ar = 0.5 * (x2.*rt1 - x1.*rt2);

        % Calculate pressure force components
        P = 0.25 * (c.P(1:end-1,1:end-1) + c.P(2:end,1:end-1) + c.P(1:end-1,2:end) + c.P(2:end,2:end));
        F = P .* Ax;

        % Calculate mass flux
        rovx_av = 0.25 * (c.rovx(1:end-1,1:end-1) + c.rovx(2:end,1:end-1) + ...
            c.rovx(1:end-1,2:end) + c.rovx(2:end,2:end));
        rovr_av = 0.25 * (c.rovr(1:end-1,1:end-1) + c.rovr(2:end,1:end-1) + ...
            c.rovr(1:end-1,2:end) + c.rovr(2:end,2:end));
        Vx_av = 0.25 * (c.Vx(1:end-1,1:end-1) + c.Vx(2:end,1:end-1) + ...
            c.Vx(1:end-1,2:end) + c.Vx(2:end,2:end));
        m = Ax.*rovx_av + Ar.*rovr_av;

        % Calculate axial momentum flux
        p = m .* Vx_av;

        % Calculate stagnation enthalpy flux
        To_av = 0.25 * (c.To(1:end-1,1:end-1) + c.To(2:end,1:end-1) + ...
            c.To(1:end-1,2:end) + c.To(2:end,2:end));    
        Ho = c.av.cp * m .* To_av;

        % Sum forces, momentum and enthalpy
        F_total = sum(F(:)); p_total = sum(p(:)); Ho_total = sum(Ho(:));

        % Calculate thrust and power
        T(o) = double(c.bv.nblade) * (F_total + p_total); 
        Wx(o) = double(c.bv.nblade) * Ho_total;
        
        % Calculate fan power from rotor row alone
        inlet = ts_structured_cut(g,2,1,1,1,'en',1,'en');
        outlet = ts_structured_cut(g,7,'en','en',1,'en',1,'en');
        [ho_in,m_in] = ts_mass_average(inlet,'ho',3);
        ho_out = ts_mass_average(outlet,'ho',3);
        Wx_fan(o) = m_in * double(inlet.bv.nblade) * (ho_out - ho_in);
        rpm(o) = inlet.bv.rpm;
        
    end

end