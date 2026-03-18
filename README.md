![image](https://github.com/zhangjm-geo/RSG_FSM/blob/main/Figures/RSGFSM.jpg)

# Description
'RSGFSM' is a novel fast sweeping method based on rotated staggered grids is proposed for solving the eikonal equation, offering both efficiency and accuracy in traveltime computation.
The conventional fast sweeping method (FSM) is based on finite-difference on standard rectangular grid, suffers severely from source singularity near the source. Consequently, errors originating at the source propagate throughout the entire computational domain during the sweeping process in off-axis positions, leading to a significant degradation in traveltime accuracy. To address this problem, a rotated staggered grid fast sweeping method (RSGFSM) is proposed. By rotating coordinates and placing traveltime and slowness on a staggered grid, their spatial separation is reduced, therefore, improving the numerical accuracy of traveltime computation. Meanwhile, additional wave propagation directions are considered in the proposed RSGFSM, thereby effectively mitigating the source singularity problem.

# Requirements
C/C++ compiler

C Make

Matlab

Linux SeisUnix (SU) (Optional)

# Usage
## compile
You can compile it directly using 'make' 
<pre><code>
make
</code></pre>

## run
and then run the executable code 'eikonal_fsm' with a parameter file 'parameter_***.txt'.
<pre><code>
./eikonal_fsm parameters_homo.txt
</code></pre>

## parameter file
A reasonable set of input parameters are as follows:
<pre><code>
--x grid number: nx 101
--z grid number: nz 101
--x grid interval (m): dx 10
--z grid interval (m): dz 10   
--x source grid position: x0 51  
--z source grid position: z0 51
--velocity model file name: velmodel ./model/homo_velocity.bin
--method(1_FSM;2_2ndFSM;3_SGFSM): method 3
--analy(1_homo;2_gradient_slow;3_no): analy 1
--velocity model: velocity 1000
--constant gradient: gradient 0.1
</code></pre>

You can choose different model parameter cards to calculate travel times for different models,
such as the homogeneous velocity model (parameters_homo.txt), 
constant gradient velocity model (parameters_gradient.txt), 
and Marmousi velocity model (parameters_mar.txt).  

For the last three parameters:

**analy**: Specifies the calculation of analytical solutions.
1: Compute the analytical solution for the homogeneous model.
2: Compute the analytical solution for the constant gradient model.
3: Do not compute analytical solutions for other models.

**velocity**: Defines the velocity value. For the homogeneous model, this is the constant velocity; for the constant gradient model, this is the velocity at the starting point.

**gradient**: Specifies the velocity gradient for the constant gradient model.

## Plot
Please use MATLAB software to open the 'MATLAB_Plot' file path and run the .m file in that folder to generate the plot.


# Example
 We test on the constant velocity model with a velocity of 1000 m/s on a 101 × 101 mesh, using 10 m spatial intervals. The source location xs is set at (500 m, 500 m). 
 <pre><code>
make
</code></pre>
 <pre><code>
./eikonal_fsm parameters_homo.txt
</code></pre>
<pre><code>
run MATLAB_Plot/homo_err1.m
run MATLAB_Plot/homo_err2.m
run MATLAB_Plot/homo_err3.m
</code></pre>
![image](https://github.com/zhangjm-geo/RSG_FSM/blob/main/Figures/homo_FSM.jpg)
![image](https://github.com/zhangjm-geo/RSG_FSM/blob/main/Figures/homo_2ndFSM.jpg)
![image](https://github.com/zhangjm-geo/RSG_FSM/blob/main/Figures/homo_RSGFSM.jpg)

Figure: The left column from top to bottom shows the comparison between the numerical solution calculated by FSM, second-order FSM, RSGFSM and the analytical solution; the middle column presents the absolute errors of the three methods relative to the analytical solution; and the right column shows the relative errors contours obtained using the three methods.


If you encounter any issues, please feel free to contact me via email at  15756461017@163.com or zhangjianming@ouc.edu.cn.
