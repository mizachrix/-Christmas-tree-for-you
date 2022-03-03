%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Energy and Reserves Dispatch with\\ Distributionally Robust Joint Chance Constraints
% Christos ORDOUDIS, Viet Anh NGUYEN, Daniel KUHN, Pierre PINSON
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Node data
N_Gas_nodes = 12;
ref_node_gas = 4;

% ADMITANCE MATRIX
B_N_gas = zeros(N_Gas_nodes, N_Gas_nodes);

% Off-diagonal elements B-matrix
for l=1:size(GasNetwork,1) % Number of Lines
    B_N_gas(GasNetwork(l,1)-100, GasNetwork(l,2)-100) = -1/GasNetwork(l,3);
    B_N_gas(GasNetwork(l,2)-100, GasNetwork(l,1)-100) = -1/GasNetwork(l,3);
end

% Diagonal elements B-matrix
for k=1:N_Gas_nodes
    B_N_gas(k,k) = -sum(B_N_gas(k,:));
end

B_L_Gas = zeros(size(GasNetwork,1), N_Gas_nodes);

% Off-diagonal elements B-matrix
for l=1:size(GasNetwork,1) % Number of Lines
    B_L_Gas(l, GasNetwork(l,1)-100) = 1/GasNetwork(l,3);
    B_L_Gas(l, GasNetwork(l,2)-100) = -1/GasNetwork(l,3);
end

% Remove ref node
B_NN_Gas = B_N_gas;
B_NN_Gas(:,ref_node_gas)=[];
B_NN_Gas(ref_node_gas,:)=[];

B_LL_Gas = B_L_Gas;
B_LL_Gas(:,ref_node_gas) = [];

PTDF_nrf_Gas = B_LL_Gas * (B_NN_Gas ^ (-1));