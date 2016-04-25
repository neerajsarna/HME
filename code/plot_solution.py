import matplotlib.pyplot as plt
from numpy import*

x = loadtxt("result_high_resolution/locationKn0p5order2",unpack=True);
solution = loadtxt("result_high_resolution/solutionKn0p5order2",unpack=True);

plt.plot(x,solution[0,:],label='density');
plt.plot(x,multiply(solution[0,:],solution[2,:]),label='pressure');

plt.legend(loc='best');
plt.xlabel("x");
plt.ylabel("u");
plt.xlim([-1,1.5]);
plt.title("variation of density and pressure");
plt.grid(True);
plt.show();


plt.plot(x,solution[2,:],label='theta');
plt.legend(loc='best');
plt.xlabel("x");
plt.ylabel("theta");
plt.title("variation of theta");
plt.grid(True);
plt.show();

plt.plot(x,solution[1,:],label='velocity');
plt.legend(loc='best');
plt.xlabel("x");
plt.ylabel("velocity");
plt.xlim([-1,1.5]);
plt.title("variation of velocity");
plt.grid(True);
plt.show();

plt.plot(x,solution[3,:],label='f3');
plt.legend(loc='best');
plt.xlabel("x");
plt.ylabel("f3");
plt.xlim([-1,1.5]);
plt.title("variation of f3");
plt.grid(True);
plt.show();

plt.plot(x,solution[4,:],label='f4');
plt.legend(loc='best');
plt.xlabel("x");
plt.ylabel("f4");
plt.xlim([-1,1.5]);
plt.title("variation of f4");
plt.grid(True);
plt.show();	
