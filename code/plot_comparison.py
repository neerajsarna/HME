import matplotlib.pyplot as plt
from numpy import*

x_DG_order1 = loadtxt("result_high_resolution/locationKn0p5order1",unpack=True);
solution_DG_order1 = loadtxt("result_high_resolution/solutionKn0p5order1",unpack=True);

x_DG_order2 = loadtxt("result_high_resolution/locationKn0p5order2",unpack=True);
solution_DG_order2 = loadtxt("result_high_resolution/solutionKn0p5order2",unpack=True);

data_PRICE = loadtxt("results_Julian/HMEKn0p5PRICELinearPath3GQ.csv",unpack=True);
data_WP2 = loadtxt("results_Julian/HMEKn0p5WP2Min.csv",unpack=True);

plt.plot(x_DG_order1,solution_DG_order1[0,:],label='density_DG_p0');
plt.plot(x_DG_order2,solution_DG_order2[0,:],label='density_DG_p1');
plt.plot(data_PRICE[0,:],data_PRICE[1,:],label='density_PRICE');
plt.plot(data_WP2[0,:],data_WP2[1,:],label='density_WP2');
plt.legend(loc='best');
plt.xlabel("x");
plt.ylabel("density");
plt.xlim([-1,1.5]);
plt.title("variation of density");
plt.grid(True);
plt.show();


plt.plot(x_DG_order1,multiply(solution_DG_order1[0,:],solution_DG_order1[2,:]),label='pressure_p0');
plt.plot(x_DG_order2,multiply(solution_DG_order2[0,:],solution_DG_order2[2,:]),label='pressure_p1');
plt.plot(data_PRICE[0,:],multiply(data_PRICE[1,:],data_PRICE[3,:]),label='pressure_PRICE');
plt.plot(data_WP2[0,:],multiply(data_WP2[1,:],data_WP2[3,:]),label='pressure_WP2');
plt.legend(loc='best');
plt.xlabel("x");
plt.ylabel("pressure");
plt.xlim([-1,1.5]);
plt.title("variation of pressure");
plt.grid(True);
plt.show();

plt.plot(x_DG_order1,solution_DG_order1[1,:],label='velocity_DG_p0');
plt.plot(x_DG_order2,solution_DG_order2[1,:],label='velocity_DG_p1');
plt.plot(data_PRICE[0,:],data_PRICE[2,:],label='velocity_PRICE');
plt.plot(data_WP2[0,:],data_WP2[2,:],label='velocity_WP2');
plt.legend(loc='best');
plt.xlabel("x");
plt.ylabel("velocity");
plt.xlim([-1,1.5]);
plt.title("variation of velocity");
plt.grid(True);
plt.show();	

plt.plot(x_DG_order1,solution_DG_order1[3,:],label='f3_DG_p0');
plt.plot(x_DG_order2,solution_DG_order2[3,:],label='f3_DG_p1');
plt.plot(data_PRICE[0,:],data_PRICE[4,:],label='f3_PRICE');
plt.plot(data_WP2[0,:],data_WP2[4,:],label='f3_WP2');
plt.legend(loc='best');
plt.xlabel("x");
plt.ylabel("f3");
plt.xlim([-1,1.5]);
plt.title("variation of f3");
plt.grid(True);
plt.show();	

plt.plot(x_DG_order1,solution_DG_order1[4,:],label='f4_DG_p0');
plt.plot(x_DG_order2,solution_DG_order2[4,:],label='f4_DG_p1');
plt.plot(data_PRICE[0,:],data_PRICE[5,:],label='f4_PRICE');
plt.plot(data_WP2[0,:],data_WP2[5,:],label='f4_WP2');
plt.legend(loc='best');
plt.xlabel("x");
plt.ylabel("f4");
plt.xlim([-1,1.5]);
plt.title("variation of f4");
plt.grid(True);
plt.show();
