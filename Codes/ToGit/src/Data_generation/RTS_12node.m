% Input data file

% Wind data

load('AV_AEMO');

% Scaling the system
Scale_Factor = 1;

% Electricity Network Data
% From       To     X       CAP 
ElNetwork=[
    101     102   0.0146    175; %l1
    101     103   0.2253    175; %l2
    101     105   0.0907    400; %l3
    102     104   0.1356    175; %l4
    102     106    0.205    175; %l5
    103     109   0.1271    400; %l6
    103     124    0.084    200; %l7
    104     109    0.111    175; %l8
    105     110    0.094    400; %l9
    106     110   0.0642    400; %l10
    107     108   0.0652    600; %l11
    108     109   0.1762    175; %l12
    108     110   0.1762    175; %l13
    109     111    0.084    200; %l14
    109     112    0.084    200; %l15
    110     111    0.084    200; %l16
    110     112    0.084    200; %l17
    111     113   0.0488    500; %l18
    111     114   0.0426    500; %l19
    112     113   0.0488    500; %l20
    112     123   0.0985    500; %l21
    113     123   0.0884    500; %l22
    114     116   0.0594    1000; %l23
    115     116   0.0172    500; %l24
    115     121   0.0249    1000; %l25
    115     124   0.0529    500; %l26
    116     117   0.0263    500; %l27
    116     119   0.0234    500; %l28
    117     118   0.0143    500; %l29
    117     122   0.1069    500; %l30
    118     121   0.0132    1000; %l31
    119     120   0.0203    1000; %l32
    120     123   0.0112    1000; %l33
    121     122   0.0692    500; %l34
    ];

% Calculation of PTDF matrix
PTDF_calc;


GenDATA=[
    0    152    40      40        1        12       12.65   12; %1
    0    152    40      40        2        13       0       0; %2
    0	 300    70      70        7        11       0       0; %3
    0	 591    60      60        13       17       0       0; %4
    0	  60    30      30        15       18       11.12   10; %5
    0	 155    30      30        15       14       0       0; %6
    0	 155    30      30        16       15       14.88   7; %7
    0	 400    50      50        18       5        0       0; %8
    0	 400    50      50        21       7        0       0; %9
    0	 300    50      50        22       20       0       0; %10
    0	 310    60      60        23       10.52    16.80   6; %11
    %0	 310    60      60        23       10.52    0       0; %11
    0	 350    40      40        23       10.89    0       0; %12
    ];

CDATA=[
    0; %1 
    44.55; %2
    65.61; %3
    30.82; %4
    0    ; %5
    33.22; %6
    0; %7
    20.84; %8
    26.90; %9
    33.25; %10
    0; %11
    %40.32; %11
    32.22; %12
    ];



%  Min      Max     Node    Cost  
GasWellData=[
    0      6000     1       2      ;
    0      6000     3       2.4    ;
    0     15000     11       2.8    ;
    ];
GasWellData(1,2) = 3500;
GasWellData(2,2) = 2500;
GasWellData(3,2) = 7500;

GasNetwork=[
    101     102    1      6000        1        1           121     40000      2; %l1
    102     104    1       6000        2        1.2           0         0      2; %l2
    103     105    1      6433        1        1           150     50000      2; %l3
    104     105    1      4479        1        1           186     60000      2; %l4
    105     106    1      4523        1        1           189     55000      2; %l5    
    104     107    1      5634        1        1           184     55000      2; %l6
    106     108    1      7513        1        1           150     40000      2; %l7
    107     108    1      2546        1        1           179     45000      1; %l8
    108     109    1       5040        2        1.3           0         0      2; %l9
    109     110    1      5040        1        1           148     40000      2; %l10
    110     111    1      4038        1        1           150     40000      1; %l11
    111     112    1      7323        1        1           130     30000      2; %l12
    ];
GasNetwork(9,4)=1500;
GasNetwork(2,4)=2000;

% Calculation of PTDF matrix
PTDF_calc_Gas;

% Demand       El_Bus  Cshed           
DemandDATA=[
    0.038      1     2000; %d1
    0.034      2     2000; %d2
    0.063      3     2000; %d3
    0.026      4     2000; %d4
    0.025      5     2000; %d5
    0.048      6     2000; %d6
    0.044      7     2000; %d7
     0.06      8     2000; %d8
    0.061      9     2000; %d9
    0.068      10    2000; %d10
    0.093      13    2000; %d11
    0.068      14    2000; %d12
    0.111      15    2000; %d13
    0.035      16    2000; %d14
    0.117      18    2000; %d15
    0.064      19    2000; %d16
    0.045      20    2000; %d17
    ];

DemandDATA(:,3) = 1000;

Total_Demand = 2650;

%   Level      Node
G_DemandDATA=[
    0.25        5; %d1
    0.25        7; %d2
    0.35        6; %d3
    0.15        12; %d4
    ];

Total_Demand_Gas = 5000;

% Wmin   Wmax   El_Bus
WindDATA=[
    0     200     1; %1 
    0     200     2; %2 
    0     200     3; %3 
    0     200     4; %4 
    0     200     5; %5 
    0     200     6; %6 
    ];

WindDATA(:,2) = 250+25;
WindDATA(:,3) = [1,2,11,12,12,16];


% Fill system info structure
system_info = [];

system_info.PTDF = round([PTDF_nrf(:,1:ref_node-1), zeros(size(PTDF_nrf,1),1), PTDF_nrf(:,ref_node:size(PTDF_nrf,2))],2); % The final PTDF matrix
system_info.PTDF_gas = round([PTDF_nrf_Gas(:,1:ref_node_gas-1), zeros(size(PTDF_nrf_Gas,1),1), PTDF_nrf_Gas(:,ref_node_gas:size(PTDF_nrf_Gas,2))],2); % The final PTDF matrix
system_info.F = ElNetwork(:,4)/Scale_Factor;
system_info.D = Total_Demand * DemandDATA(:,1)/Scale_Factor;
system_info.Pmax = GenDATA(:,2)/Scale_Factor;
system_info.Pmin = GenDATA(:,1)/Scale_Factor;
system_info.R = (system_info.Pmax + system_info.Pmin)/2;
system_info.ResCap = GenDATA(:,2)*0.40;

system_info.Gmax = GasWellData(:,2);
system_info.Gmin = GasWellData(:,1);
system_info.Dg = Total_Demand_Gas * G_DemandDATA(:,1)/Scale_Factor;
system_info.FG = GasNetwork(:,4)/Scale_Factor;

system_info.Wmax = WindDATA(:,2)/Scale_Factor;
system_info.DiagWmax = diag(system_info.Wmax);
system_info.Wmin = WindDATA(:,1)/Scale_Factor;

%system_info.C = 0.5*round(GenDATA(:,6)*Scale_Factor,3);
system_info.Clsh = DemandDATA(:,3)*Scale_Factor;
system_info.C = CDATA(:,1);%[35;40;30;55;60;45;50;10;15;65;20;25]/2;
system_info.A = diag(0.0001*system_info.C);
% system_info.Cr = [60;70;50;100;110;80;90;10;20;120;30;40];
system_info.Cr1 = 0.2*system_info.C;
system_info.Cr2 = 0.2*system_info.C;

system_info.Cg = GasWellData(:,4);
system_info.Ag = diag(0.0001*system_info.Cg);

N_Gas_nodes = 12;
system_info.NGN = N_Gas_nodes;

% Mapping on the network
system_info.AG = zeros(N_El_nodes, size(GenDATA,1));
system_info.AW = zeros(N_El_nodes, size(WindDATA,1));
system_info.AD = zeros(N_El_nodes, size(DemandDATA,1));
system_info.PG = zeros(N_Gas_nodes, size(GenDATA,1));
system_info.GG = zeros(N_Gas_nodes, size(GasWellData,1));
system_info.DG = zeros(N_Gas_nodes, size(G_DemandDATA,1));
system_info.PLG = zeros(N_Gas_nodes, size(GasNetwork,1)); %system_info.PTDF_gas;%




for n=1:N_El_nodes
    for gg=1:size(GenDATA,1)
        if GenDATA(gg,5) == n
        system_info.AG(n,gg) = 1;
        end
    end
    for ww=1:size(WindDATA,1)
        if WindDATA(ww,3) == n
        system_info.AW(n,ww) = 1;
        end
    end
    for dd=1:size(DemandDATA,1)
        if DemandDATA(dd,2) == n
        system_info.AD(n,dd) = 1;
        end
    end
end



for m = 1:N_Gas_nodes
    for ggg = 1:size(GenDATA,1)
        if GenDATA(ggg,8) == m
            system_info.PG(m,ggg) = 1;
        end
    end
    for gggg = 1:size(GasWellData,1)
        if GasWellData(gggg,3) == m
            system_info.GG(m,gggg) = 1;
        end
    end
    for ddd=1:size(G_DemandDATA,1)
        if G_DemandDATA(ddd,2) == m
        system_info.DG(m,ddd) = 1;
        end
    end    
    for ppp=1:size(GasNetwork,1)
        if GasNetwork(ppp,1)-100 == m 
        system_info.PLG(m,ppp) = -1;
        end
        if GasNetwork(ppp,2)-100 == m
        system_info.PLG(m,ppp) = 1;        
        end
    end       
end
% system_info.PLG = system_info.PLG;
system_info.phi = GenDATA(:,7)' * 0.5;
system_info.phi2 = [1;5;7;11];


system_info.IG = zeros(size(system_info.Pmax,1),nnz(system_info.phi));

for i = 1:size(GenDATA,1)
    for j = 1:nnz(system_info.phi)
        if system_info.phi2(j,1) == i
            system_info.IG(i,j) = 1;
        end
    end
end


% PTDF and mapping matrices combination

system_info.Qg = system_info.PTDF*system_info.AG;
system_info.Qw = system_info.PTDF*system_info.AW;
system_info.Qd = system_info.PTDF*system_info.AD;

