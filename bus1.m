% busdata: [bus#, type(1=slack),Pg,Qg Pd, Qd, Vm_spec, Va_spec]
function busdata=bus1()
busdata = [
    1 1  0   0  0     0    1.06  0;
    2 2  0.4 0  0.2  0.1   1    0;
    3 2  0   0  0.45 0.15  1    0;
    4 2  0   0  0.4  0.05  1    0;
    5 2  0   0  0.6  0.1   1    0;
];
