# Shaft-Calculator-Benny
The critical points on a shaft include the points of contact that gears make with each other and their respective loadings. Another critical point on a shaft is where the stress concentration may occur, the connections between varying diameters or sections of the shaft. For this problem, the shaft calculation code will consider the inputs given by the user and assumptions for the user to not have to input. The calculations made within the code will iterate a smaller diameter for the user to find all of the following outputs. This code was created in response to the assignment given by Professor Yi ren from Arizona State University for his MEE 342: Principles of Mechanical Design course. The outputs and inputs are based on the project given to design a 2-speed gearbox/transmission. Below are the inputs required for the calculation:

Input:
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


Output:
D: Bigger diameter of the segment
d: smaller diameter of the segment
r: filet radius 
kt: stress concentration factor, using Table 7-1 an assumption is made
kts: stress concentration factor for torsion, using Table 7-1 an assumption is made
kf: stress concentration due to fatigue from Moment
kfs: stress concentration due to fatigue from Torsion

Example: 
For this example, we are to use AISI 4340 Alloy Steel, Sut = 108 [ksi] and Sy = 68.2 [ksi]. The shaft is to have a safety of factor of 2, n = 2. A reliability of 99% and 10 year continuous operation. The maximum torque and moment that the shaft undergoes is 1260 [lbf*in] and 1100 [lbf*in] respectfully (numbers pulled from Shigley’s Example 7-1). Using the assumption that the smaller and bigger diameter of the shaft has a ratio of D/d =1.1, the smaller diameter is to be determined. For the first iteration the filet radius ratio is to be determined as r/d =.05. The stress concentration due to static and fatigue failure are to be determined as well for the rest of the calculations needed for the gearbox design. The key calculations will follow the shaft calculations to determine the best design.  
https://www.azom.com/article.aspx?ArticleID=6772 


First the endurance limit is determined for the shaft using the assumptions made due to lack of information and/or the givens. From Shigley’s, the endurance limit is calculated using the following. 

Typically, the shaft is made of a machined and/or cold-drawn material. For this example we are to use the assumption that it is a cold-drawn material and the first iteration assumptions. 

First, the surface condition modification factor is to be determined using the following table. 

ka=a(Sutb)
a= 2.70 [kpsi]
b= -0.265
ka=2.70(108-.265)=.78074

The size modification factor is to be determined using the following. Typically, gearbox or transmission shafts contain shafts that have a diameter of 0.11  d  2 [in]. For our first guess or iteration, we are to assume that the first diameter is 1 [in]. 

kb=0.879d-0.107
kb=0.897(1)-0.107=1

The load modification factor we are going to use the value for bending using the following table. 

kc=1

The temperature modification factor can be found using the following. 

kd=0.975+0.432(10-3)*Tf-0.115*(10-5)Tf2+0.104(10-8)Tf3-0.595(1012)Tf4

Using an average maximum temperature that a shaft undergoes, 482 [K], we can find the following. This is also not specified, it can be assumed that it is 1. It is also found using this temperature specified, the temperature modification factor is the following.
kd=1

Now we are to determine the reliability factor. We are to use the specified value of 99%. 

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

The axial loads are not considered for the rest of the calculations from here due to the bearings supporting the shaft due to the reaction forces. The torsional and moment loads are primarily considered. 
1n=a'Se+m'Sut

For the initial iteration we are going to assume that the shaft diameter ratio will be 1.1. The radius ratio will be 0.1.
Dd=1.1
rd=0.05

The assumption of Kts = Kfs and Kt = Kf, should not be considered to get the most accurate results possible. The iteration will compare the calculation of the factor of safety and the smaller diameter and see if they correlate with each other. Now the fatigue stress concentration values are found using the following:




Now, for the following iterations, it is required to recalculate the stress concentration using the smaller diameter found in each respective iteration. The static stress concentration factors must be iterated using the table, in order to complete this an external resource was used outside of Shigley’s. The case study on how the calculator fits with the tables provided in the Shigley’s is provided below the output summary. 
kb=0.879d-0.107

Now, using the fatigue stress concentration factors, the design factor is going to be verified for each. The diameter is iterated over a while loop to determine the final diameter that will be calculated 

a'=32KfMd3
m'=3KfsTd3


Here is the final use case via MATLAB:

Output:

Dn, is the smaller diameter. D is the bigger diameter. Rn is the filet radius. Kt and Kts are the static stress concentration factors due to bending and torsion respectively. Similarly, Kf and Kfs are the fatigue stress concentration factors due to bending and torsion respectively. 

Case Study in regards to the table fitting:
** Below is a conversation between the professor and the creator of the code**
I found that the following site: https://amesweb.info/stress-concentration-factor-calculator/formulas.aspx , I found that the error between Table A-15-8/9 and the equations provided have a near negligible difference between the two. After taking four points, I found that the difference of the theoretical stress-concentration factor due to bending (K_t) is ~1.21%. I found that the difference between the theoretical stress-concentration factor due to torsion is ~7.16%. Taking the average of the two, it is a 4.185% difference. 

I am going to accept this difference because of the fact that we are "eye-balling" a graph, which propagates precision error. This graph that is displayed in Shigley's doesn't give exact values for every point as well. Thus, I will proceed with my automation code using the website's formulas. 

Below are the relations, derivations and calculations that we conducted together. 

Links to the desmos calculations: 
https://www.desmos.com/calculator/psmktvbova 
https://www.desmos.com/calculator/7pqkk7286j


