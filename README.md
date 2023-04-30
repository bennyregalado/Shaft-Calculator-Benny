# Shaft-Calculator-Benny
# What it does 
  There are several assumptions that are made in this script immediately after the inputs are defined. At the interactions between the Gears, it is assumed that the tangential components on each Gear are equal due to Newton’s third law. Additionally, it is assumed that there is a Pressure angle between Gear interactions of 20 degrees clockwise from the positive x direction. The main nomenclature of this code is based on a hand drawn Free body diagram. Locations of the Forces are at the midpoints of the Bearings and Gears. Gear forces are represented with the letters C, F, G and J. Bearing Forces are represented by A, B, D, E, H, and I. The code is written based on the input being the right side (Shaft 1), the Output being the left side (Shaft 3) and the intermediate shaft in the middle (Shaft 2). Please refer to the picture of the FBD for more clarification.
  
  <a href="https://ibb.co/pn3vj5q"><img src="https://i.ibb.co/b7rg5Dt/Screenshot-2023-04-23-232940.jpg" alt="Screenshot-2023-04-23-232940" border="0"></a>
  
  Based on these assumptions, the code will begin conducting Static Analyses by determining the unknown Forces through a series of Sum of Forces and Moments calculations on each the xy and yz plane for each shaft. The Torques were also calculated as they go through each shaft based on the Gear ratios. The Lengths of each shaft defined by the numbers in the FBD. These numbers were calculated and referenced figure 18-2 in Shigley’s. 
  
  After all the Forces and Torques have been calculated. A series of Vectors for Forces, Moments and Torques have been initialized and defined. The size of these vectors are made so that each index represents 0.01 inches along the length of each shaft. Using these defined Force vectors, the Moment vectors were defined by performing a Trapezoidal Sum integration within a for loop. This same procedure was done for the other two Shafts. 
  
  Once all the Forces, Moments and Torques have been determined along the length of each Shaft, Fatigue Analysis will be performed. A separate script of code has been formatted to operate as a function call to determine the smaller and larger diameter of each shaft. The inputs for the function call are; the largest Moment, Torque and the Material properties Yield Strength and Ultimate Tensile Strength. This will output a smaller diameter, larger Diameter, radius of the notch and yielding Factor of Safety.
  
  Assuming a fatigue factor of safety of 2 and an initial diameter of 0.5 in, The fatigue code starts with the following:
 
 Calculate Ka using the a and b values for Cold drawn/machined materials
  
  Kb using the initial diameter
  
  Kc and kf is assumed to be one
  
  Kd is calculated by assuming a maximum temperature of 482 F
  
  Ke is 0.814 for 99% Reliability
    
  Endurance limit is calculated for the Material’s Ultimate Tensile Strength and multiplied by all the modification factors
Afterwards, it assumes that the Static Notch Stress concentration for Bending is 1.7 and equal to the Fatigue Stress Concentration (i.e. kt = kf = 1.7). The Shear Notch Stress Concentration for Torsion is assumed to be 1.5 and equal to the Fatigue Stress Concentration as an assumption as well (i.e. kts = kfs = 1.5). Using and reorganizing the Distortion Energy Goodman Equation, the new outer diameter and radius was calculated. Then, The Neuber constant (sqrt_a) was calculated for both Bending and shear using the Ultimate tensile Strength. 
	  
  The While loop is used to calculate the larger diameter, smaller diameter and radius of the Notch within the constraints until the difference between diameters between iterations are less than 0.001 inches. The constraints are that the larger diameter is 1.1 times the smaller diameter and the notch radius is 5% of the larger diameter. This while loop uses multiple complex equations to automate and approximate the Static Stress Concentration for Bending and Shear Stress. It then calculates the notch sensitivity factors and then the Fatigue Stress Concentration Factors. Adjusting the Size factor, the Endurance limit was recalculated and then the Distortion Energy Goodman equation was used again to calculate the new diameter. This loop will continue until the error criteria is met. This was done for each shaft so there will be a total of 6 diameters calculated (2 different ones on each shaft).
    
  The Deflection Calculations were done after Fatigue Analysis. It starts by averaging the Moment vectors on each plane to create a maximum Moment vector for each shaft. Then using an input of different diameters and location of each diameter change, Inertia vectors were calculated. Then a series of for loops were created to manually calculate Inertia, Theta and deflection vectors with respect to the size of each Shaft and with each index representing 0.01 inches of the Shaft. A series of Checks were performed at the edges of the Gears and Bearings to determine whether the current design passes the Maximum deflection criteria for Gears and Maximum Slope criteria for Bearings. The results are just printed onto the Command Window. Locations are again based on the figure 18-2 in Shigley’s. 
	  
  There are 3 final sections to this main script which is a print out of all the Force, Moment, and Torque diagrams for each plane and on each Shaft. This should be read from INPUT SHAFT TO OUTPUT SHAFT. Ex, At one inch on the Force Diagram for shaft one, the Force at that point is 1 inch away from where the Loads are applied.

# Inputs:
  Material Properties of material of choice 
  
  sut: Ultimate Tensile Strength [ksi] (line 10)
  
  sy: Yielding Strength [ksi] (line 11)
  
  E: Elastic Modulus [psi] (line 12)
  
  Geometry
  
  Dc: diameter of Gear 2 (N2) (line 13)
  
  Df: diameter of Pinion 3 (N3) (line 13)
  
  Dg: diameter of Gear 4 (N4) (line 13)
  
  Dj: diameter of Pinion 5 (N5) (line 13)
  
  GR_CF: gear ratio of first set (line 14)
  
  GR_JG: gear ratio of second set (line 15)

*Note: The lengths of the shaft are fixed within the script. If the user wants to change the lengths, throughout the script the lengths of the script, there are lengths indicated on each respective line that the lengths are utilized for each section.

# Outputs
  d1: The diameter of each section of the input shaft found in the workspace denoted as vector
  
  d2: The diameter of each section of the intermediate shaft found in the workspace denoted as vector
  
  d3: The diameter of each section of the output shaft found in the workspace denoted as vector
  
  Yielding and Fatigue Factor of safety found in each section of the shaft found in the command window
  
  Filet radius [in] for each section of the shaft found in the command window
  
  Deflection and slope checks for each shaft based on a maximum deflection of 0.003 [in] and slope of 0.0005 [rad]. This will be displayed in the command window.
  
  Rotational speeds of each shaft are displayed in the command window.
  Shear, moment and torque diagrams for each respective shaft in each plane (YX and ZX)

# Summary of Fatigue Analysis
The critical points on a shaft include the points of contact that gears make with each other and their respective loadings. Another critical point on a shaft is where the stress concentration may occur, the connections between varying diameters or sections of the shaft. For this problem, the shaft calculation code will consider the inputs given by the user and assumptions for the user to not have to input. The calculations made within the code will iterate a smaller diameter for the user to find all of the following outputs. This code was created in response to the assignment given by Professor Yi ren from Arizona State University for his MEE 342: Principles of Mechanical Design course. The outputs and inputs are based on the project given to design a 2-speed gearbox/transmission. Below are the inputs required for the calculation:

# Input:
M: the maximum moment on the shaft [lbf*in]

T: the maximum torque on the shaft [lbf*in]

Sut: The ultimate tensile strength of the material chosen [psi]


Assumptions of Calculator:

a: factor for calculating ka, for machined/cold-drawn material 

b: factor for calculating ka, for machined/cold-drawn material 

n: factor of safety is 2

D/d = 1.1, the bigger and smaller diameter ratio used to minimize material and weight.

r/d = .05, the radius and smaller diameter ratio used to minimize material and weight.

kd: temperature modification factor, assuming that the maximum temperature is experienced for shaft (~482 F or 250 C)

ke: reliability factor, assuming 99% reliability 

kf: miscellaneous-effects modification factor, assumption of 1, is only added in special cases. It comes from experience as a design engineer.


# Output:
D: Bigger diameter of the segment
d: smaller diameter of the segment
r: filet radius 
kt: stress concentration factor, using Table 7-1 an assumption is made
kts: stress concentration factor for torsion, using Table 7-1 an assumption is made
kf: stress concentration due to fatigue from Moment
kfs: stress concentration due to fatigue from Torsion

# Example: 
  For this example, we are to use AISI 4340 Alloy Steel, Sut = 108 [ksi] and Sy = 68.2 [ksi]. The shaft is to have a safety of factor of 2, n = 2. A reliability of 99% and 10 year continuous operation. The maximum torque and moment that the shaft undergoes is 1260 [lbf*in] and 1100 [lbf*in] respectfully (numbers pulled from Shigley’s Example 7-1). Using the assumption that the smaller and bigger diameter of the shaft has a ratio of D/d =1.1, the smaller diameter is to be determined. For the first iteration the filet radius ratio is to be determined as r/d =.05. The stress concentration due to static and fatigue failure are to be determined as well for the rest of the calculations needed for the gearbox design. The key calculations will follow the shaft calculations to determine the best design.  
https://www.azom.com/article.aspx?ArticleID=6772 


  First the endurance limit is determined for the shaft using the assumptions made due to lack of information and/or the givens. From Shigley’s, the endurance limit is calculated using the following.

<a href="https://imgbb.com/"><img src="https://i.ibb.co/9tghhxL/Screenshot-2023-04-23-233051.jpg" alt="Screenshot-2023-04-23-233051" border="0"></a>

  Typically, the shaft is made of a machined and/or cold-drawn material. For this example we are to use the assumption that it is a cold-drawn material and the first iteration assumptions. 

First, the surface condition modification factor is to be determined using the following table. 

ka=a(Sutb)

a= 2.70 [kpsi]

b= -0.265

ka=2.70(108-.265)=.78074

The size modification factor is to be determined using the following. Typically, gearbox or transmission shafts contain shafts that have a diameter of 0.11  d  2 [in]. For our first guess or iteration, we are to assume that the first diameter is 1 [in]. 
<a href="https://imgbb.com/"><img src="https://i.ibb.co/Zhn0yNf/Screenshot-2023-04-23-233201.jpg" alt="Screenshot-2023-04-23-233201" border="0"></a>

kb=0.879d-0.107

kb=0.897(1)-0.107=1

The load modification factor we are going to use the value for bending using the following table. 

kc=1

The temperature modification factor can be found using the following. 

kd=0.975+0.432(10-3)*Tf-0.115*(10-5)Tf2+0.104(10-8)Tf3-0.595(1012)Tf4

Using an average maximum temperature that a shaft undergoes, 482 [K], we can find the following. This is also not specified, it can be assumed that it is 1. It is also found using this temperature specified, the temperature modification factor is the following.
kd=1

Now we are to determine the reliability factor. We are to use the specified value of 99%. 

<a href="https://ibb.co/TTfVLZD"><img src="https://i.ibb.co/sPf0Hhp/Screenshot-2023-04-23-233649.jpg" alt="Screenshot-2023-04-23-233649" border="0"></a>

ke=0.814

Since there is no information that pertains to the miscellaneous effects on the shaft, the respective factor can be specified as the following. 

kf=1

Now, the rotary-beam test specimen endurance limit is to be determined using the following. 

Se'=Sut2

It is safe to say that for our first assumption, our fatigue stress concentration due to bending and torsion are equal to the static equivalents. For our first iteration we are going to use the following numbers based on Table provided in Shigley’s 

Kf=Kt

Kfs=Kts

Kf=Kt=1.7

Kfs=Kts=1.5

<a href="https://ibb.co/pZ2BsMM"><img src="https://i.ibb.co/Qn93qTT/Screenshot-2023-04-23-233746.jpg" alt="Screenshot-2023-04-23-233746" border="0"></a>

The axial loads are not considered for the rest of the calculations from here due to the bearings supporting the shaft due to the reaction forces. The torsional and moment loads are primarily considered. 

1/n=a'Se+m'Sut

<a href="https://imgbb.com/"><img src="https://i.ibb.co/swPTrDW/Screenshot-2023-04-23-233858.jpg" alt="Screenshot-2023-04-23-233858" border="0"></a>

For the initial iteration we are going to assume that the shaft diameter ratio will be 1.1. The radius ratio will be 0.1.

Dd=1.1

rd=0.05

The assumption of Kts = Kfs and Kt = Kf, should not be considered to get the most accurate results possible. The iteration will compare the calculation of the factor of safety and the smaller diameter and see if they correlate with each other. Now the fatigue stress concentration values are found using the following:

<a href="https://imgbb.com/"><img src="https://i.ibb.co/Lt152bW/Screenshot-2023-04-23-234008.jpg" alt="Screenshot-2023-04-23-234008" border="0"></a>

Now, for the following iterations, it is required to recalculate the stress concentration using the smaller diameter found in each respective iteration. The static stress concentration factors must be iterated using the table, in order to complete this an external resource was used outside of Shigley’s. The case study on how the calculator fits with the tables provided in the Shigley’s is provided below the output summary. 

kb=0.879d-0.107

<a href="https://ibb.co/z7JMV9y"><img src="https://i.ibb.co/nnCGMdH/Screenshot-2023-04-23-234054.jpg" alt="Screenshot-2023-04-23-234054" border="0"></a>

Now, using the fatigue stress concentration factors, the design factor is going to be verified for each. The diameter is iterated over a while loop to determine the final diameter that will be calculated 

<a href="https://imgbb.com/"><img src="https://i.ibb.co/0BLQ7d2/Screenshot-2023-04-23-234151.jpg" alt="Screenshot-2023-04-23-234151" border="0"></a>

sigma_a'=32KfMd3

sigma_m'=3KfsTd3

Here is the final use case via MATLAB:

<a href="https://imgbb.com/"><img src="https://i.ibb.co/x1FjZDt/Screenshot-2023-04-23-234234.jpg" alt="Screenshot-2023-04-23-234234" border="0"></a>

Output:

![image](https://user-images.githubusercontent.com/122718640/233919004-0a349b1c-9368-404a-986e-7c59b2815555.png)

Dn, is the smaller diameter. D is the bigger diameter. Rn is the filet radius. Kt and Kts are the static stress concentration factors due to bending and torsion respectively. Similarly, Kf and Kfs are the fatigue stress concentration factors due to bending and torsion respectively. 

Case Study in regards to the table fitting:

** Below is a conversation between the professor and the creator of the code**

I found that the following site: https://amesweb.info/stress-concentration-factor-calculator/formulas.aspx , I found that the error between Table A-15-8/9 and the equations provided have a near negligible difference between the two. After taking four points, I found that the difference of the theoretical stress-concentration factor due to bending (K_t) is ~1.21%. I found that the difference between the theoretical stress-concentration factor due to torsion is ~7.16%. Taking the average of the two, it is a 4.185% difference. 

I am going to accept this difference because of the fact that we are "eye-balling" a graph, which propagates precision error. This graph that is displayed in Shigley's doesn't give exact values for every point as well. Thus, I will proceed with my automation code using the website's formulas. 

Below are the relations, derivations and calculations that we conducted together. 

<a href="https://ibb.co/cCbjM5w"><img src="https://i.ibb.co/5nYPtDk/Screenshot-2023-04-23-234435.jpg" alt="Screenshot-2023-04-23-234435" border="0"></a>

<a href="https://ibb.co/Rg3VNLT"><img src="https://i.ibb.co/j8y1z7k/Screenshot-2023-04-23-234520.jpg" alt="Screenshot-2023-04-23-234520" border="0"></a>

Links to the desmos calculations: 
https://www.desmos.com/calculator/psmktvbova 
https://www.desmos.com/calculator/7pqkk7286j

end of email

