#include <iostream>
#include <vector>
#include <cmath>
#include <iterator>
#include <algorithm>
#include "TdS_Math.h"
#include <fstream>
#include <sstream>
#include <ctime>
#include <random>

#include "poly34_rev.h"

using namespace std;

//const double EPS = 1e-14;
#define pi 3.14159265358979323846264338327950288
#define EPS (1e-14)
#define MaxNp 100
#define SAVE_SIMULATION false
#define SAVE_TABLE true
#define CORRELATED_NOISE false
#define VertDiff 0


char fout_name[100];

int P3P_MinkovSawada(double angleA, double angleB, double angleAB, double angleBC, double angleCA, double(&interpretations)[10][3]);
int SrereoEyeRecovery_Optimization(double *eyeL, double *eyeR, int numPoints, double *px, double *py, double *pz,  double noise, double ang_condition);
double SrereoEyeRecovery_Cost(double degA, double degB, double(&visual_angles)[3][MaxNp][2], int numPoints,	// Step 1: Recover 3D scene & Compute cost
	double(&distance_veridical)[MaxNp][2], double(&vergence_veridical)[MaxNp], double *diff_dist, double *diff_verg);	// Step 2: Compare recovered 3D scene with veridical 3D scene

int ComputeAnotherProjectionLine(double *vq, double *eye, double *pa, double *pb, double *pc, double anga, double angb, double angc, bool raddeg);
int SolvePolyEquation(double *x, double coef4, double coef3, double coef2, double coef1, double coef0);

int main() {

	int random_seed = (unsigned int)time(0);
	srand(40);
	srand(random_seed);

	clock_t time0 = clock();
	
	//2 eyes recovery
	double eyeL[3] = { -3.3, 0, 0 };
	double eyeR[3] = { +3.3, 0, 0 };


	/* Simulation 1 * 
	float tempf;
	FILE *setting_experiment;
	fopen_s(&setting_experiment, "Setting_experiment.txt", "r");
	fscanf_s(setting_experiment, "%d", &num_repeat);
	fscanf_s(setting_experiment, "%d", &numPoints);
	if (numPoints > MaxNp)
	{
		printf("Too many points\n\n");
		getchar();
		return 0;
	}
	fscanf_s(setting_experiment, "%f", &tempf); base_noise = tempf;
	fscanf_s(setting_experiment, "%f", &tempf); base_minZ = tempf;
	fclose(setting_experiment);
	/**/


	/* Simulation 2 *
	double px[MaxNp], py[MaxNp], pz[MaxNp];
	px[0] = -20;	py[0] = -20;	pz[0] = 57;
	px[1] = -20;	py[1] = +20;	pz[1] = 57;
	px[2] = +20;	py[2] = -20;	pz[2] = 57;
	px[3] = +20;	py[3] = +20;	pz[3] = 57;
	px[4] = 0;		py[4] = 0;	    pz[4] = 57;
	int numPoints = 5;
	/**/

	/**
	int error_recovery = SrereoEyeRecovery_Optimization(eyeL, eyeR, numPoints, px, py, pz, 0, 0);
	if (error_recovery != 0)
		printf("Error type: %d\n", error_recovery);
	else
		printf("Recovered\n");
	/**/


	/* Simulation 3 */
	int num_repeat = 100; //temp
	int numPoints = 5; //temp
	int base_minZ = 40;
	int minZ = base_minZ;

	double base_noise = 0;
	double noise = base_noise;
	double base_ang = 0;
	double ang_condition = base_ang;
	
	double px[MaxNp], py[MaxNp], pz[MaxNp];
	ofstream fout;

	//10-20 (r=10); 20-40 (r=20); 40-80 (r=40); 80-160 (r=80); 160-320 (r=160)
	//Np = 5,6,7,8,...
	for (int sessionA = 0; sessionA < 1; sessionA++)
	{
		minZ = base_minZ;
		for (int sessionB = 0; sessionB < 1; sessionB++) // minZ //6
		{
			if (SAVE_SIMULATION)
			{
				//sprintf_s(fout_name, "Output_data_%3d_%3d.txt", numPoints, minZ);
				sprintf_s(fout_name, 100, "Output_data_nn%3d_%3d.txt", numPoints, minZ);
				//fout.open("Output_data.txt", ios::out | ios::app);
				fout.open(fout_name, ios::out | ios::app);
				fout << "Repeat\t" << num_repeat << endl;
				fout << "VertDiff\t" << VertDiff << endl;
				fout << "LR correlated noise\t" << CORRELATED_NOISE << endl;
				fout.close();
			}

			noise = base_noise;
			ang_condition = base_ang;
			for (int sessionC = 0; sessionC < 1; sessionC++) // noise or ang_condition
			{
				ofstream fout;
				if (SAVE_SIMULATION)
				{
					//fout.open("Output_data.txt", ios::out | ios::app);
					fout.open(fout_name, ios::out | ios::app);
					fout << endl;
					fout << "Points\t" << numPoints << endl;
					fout << "Noise\t" << noise << endl;
					fout << "MinZ\t" << minZ << endl;
					fout << "Ang\t" << ang_condition << endl;
					fout.close();
				}

				for (int trial = 0; trial < num_repeat; trial++)
				{
					int error_recovery;

					int rseed = (unsigned int)time(0);
					srand(rseed);

					do {
						for (int p = 0; p < numPoints; p++)
						{
							pz[p] = TdS_RandomNumber(minZ, minZ * 2);

							double eccentricity = TdS_RandomNumber(0, 45) * pi / 180.0;
							double ori_xy = TdS_RandomNumber(0, 360) * pi / 180.0;
							px[p] = pz[p] * tan(eccentricity) * cos(ori_xy);
							py[p] = pz[p] * tan(eccentricity) * sin(ori_xy);
						}
						error_recovery = SrereoEyeRecovery_Optimization(eyeL, eyeR, numPoints, px, py, pz, noise, ang_condition);
					} while (error_recovery == -2 || error_recovery == -3);

					if (error_recovery != 0)
						printf("Error type: %d\n", error_recovery);
					else
						printf("Trial %d\n", trial);
				} // trial
				//noise /= 2.0;
				//ang_condition *= 2.0;
			} //sessionC
			minZ *= 4;

		} //sessionB
		numPoints *= 2;

	} //sessionA
	/**/

	clock_t time1 = clock();
	cout << double(time1 - time0)/CLOCKS_PER_SEC << endl;
	
	printf("Program end\n\n");
	getchar();

	return 0;
}

int P3P_MinkovSawada(double angleA, double angleB, double angleAB, double angleBC, double angleCA, double(&interpretations)[10][3])
{
	interpretations[0][0] = 0; interpretations[0][1] = 0; interpretations[0][2] = 0;
	interpretations[1][0] = 0; interpretations[1][1] = 0; interpretations[1][2] = 0;
	interpretations[2][0] = 0; interpretations[2][1] = 0; interpretations[2][2] = 0;
	interpretations[3][0] = 0; interpretations[3][1] = 0; interpretations[3][2] = 0;

    // Valid tetrahedron
    if ((angleAB <= 0) || (180 <= angleAB) || (angleBC <= 0) || (180 <= angleBC) || (angleCA <= 0) || (180 <= angleCA))
		return -1;
    if ((angleAB + angleBC < angleCA) || (angleBC + angleCA < angleAB) || (angleCA + angleAB < angleBC))
		return -2;
    if ((angleA <= 0) || (180 <= angleA) || (angleB <= 0) || (180 <= angleB))
		return -3;

    // Valid triangle
    if (180 <= angleA + angleB)
		return -4;
    if (360 <= angleAB + angleBC + angleCA)
		return -5;

    double angleC = 180 - (angleA + angleB);

	// Valid triangle and valid image but invalid tetrahedron
	//  #interpretations = 0
	if ((180 - angleAB) + (180 - angleBC) < angleB)
		return 0;
	if ((180 - angleBC) + (180 - angleCA) < angleC)
		return 0;
	if ((180 - angleCA) + (180 - angleAB) < angleA)
		return 0;





	// New Jun/28/2018
	double angAB, angBC, angCA, angA, angB, angC;
	double tempA[3] = { angleA, angleB, angleC };
	double tempB[3] = { angleB, angleC, angleA };
	double tempC[3] = { angleC, angleA, angleB };
	double tempAB[3] = { angleAB, angleBC, angleCA };
	double tempBC[3] = { angleBC, angleCA, angleAB };
	double tempCA[3] = { angleCA, angleAB, angleBC };
	double interpretationsX3[3][8][3] = {};
	int numInterpretationsX3[3] = {};

	for (int ang = 0; ang < 3; ang++)
	{
		angA = tempA[ang];  angB = tempB[ang];  angC = tempC[ang];
		angAB = tempAB[ang]; angBC = tempBC[ang]; angCA = tempCA[ang];

		double cosAB = cos(pi*angAB / 180.0);
		double cosBC = cos(pi*angBC / 180.0);
		double cosCA = cos(pi*angCA / 180.0);
		double cosA = cos(pi*angA / 180.0);
		double sinA = sin(pi*angA / 180.0);
		double cosB = cos(pi*angB / 180.0);
		double sinB = sin(pi*angB / 180.0);

		double Rab = 100;
		double Rbc = Rab / (cosB + sinB * cosA / sinA);
		double Rca = Rab / (cosA + sinA * cosB / sinB);

		//double tanABx, tanBCx, tanCAx;
		//double ratio_ang_R[3], ratio_max;
		//int i_max;
		//tanABx = (angAB<=84 ? tan(pi*angAB / 180.0) : 10);
		//tanBCx = (angBC<=84 ? tan(pi*angBC / 180.0) : 10);
		//tanCAx = (angCA<=84 ? tan(pi*angCA / 180.0) : 10);
		//ratio_ang_R[0] = Rab / tanABx;
		//ratio_ang_R[1] = Rbc / tanBCx;
		//ratio_ang_R[2] = Rca / tanCAx;

		//ratio_max = TdS_Max(3, ratio_ang_R, &i_max);

		//Rab /= ratio_max/100.0;
		//Rbc /= ratio_max/100.0;
		//Rca /= ratio_max/100.0;
		
		double K1 = ((pow(Rbc, 2)) / (pow(Rca, 2)));
		double K2 = ((pow(Rbc, 2) / pow(Rab, 2)));
		double K1K2 = K1*K2;
		double cos2CA = cosCA*cosCA;
		double cos2BC = cosBC*cosBC;

		double G4 = pow(K1K2 - K1 - K2, 2) - 4 * K1K2*cos2BC;
		double G3 = 4 * (K1K2 - K1 - K2)*(K2 - K1K2)*cosAB
			+ 4 * K1*cosBC*((K1K2 + K2 - K1)*cosCA + 2 * K2*cosAB*cosBC);
		double G2 = pow(2 * (K2 - K1K2)*cosAB, 2)
			+ 2 * (K1K2 + K1 - K2)*(K1K2 - K1 - K2)
			+ 4 * K1*((K1 - K2)*cos2BC + (K1 - K1K2)*cos2CA - 2 * (K2 + K1K2)*cosAB*cosCA*cosBC);
		double G1 = 4 * (K1K2 + K1 - K2)*(K2 - K1K2)*cosAB
			+ 4 * K1*((K1K2 - K1 + K2)*cosCA*cosBC + 2 * K1K2*cosAB*cos2CA);
		double G0 = pow(K1K2 + K1 - K2, 2) - 4 * K1*K1K2*cos2CA;

		double x[4] = {};
		int numPolySolutions;
		numPolySolutions = SolvePolyEquation(x, G4, G3, G2, G1, G0);
		
		if (numPolySolutions <= 0)
			continue;
		
		for (int i = 0; i < numPolySolutions; i++)
		{
			if (x[i] <= 0) continue;

			if (pow(x[i], 2) - 2 * x[i] * cosAB + 1 <= 0) continue;

			double a = Rab / sqrt(pow(x[i], 2) - 2 * x[i] * cosAB + 1);  // A24
			double b = x[i] * a;
			// a > 0 and b > 0 because Rab > 0

			if (angA == angBC && a < 1e-10) continue;
			if (angB == angCA && b < 1e-10) continue;

			double y0 = cosCA*cosCA + (Rca*Rca - a*a) / (a*a); // Equation A?
			if (y0 < 0) continue;


			double m1 = 1 - K1;
			double p1 = 2 * (K1 * cosCA - x[i] * cosBC);
			double q1 = (pow(x[i], 2) - K1);
			double m2 = 1;
			double p2 = 2 * (-x[i] * cosBC);
			double q2 = pow(x[i], 2) * (1 - K2) + 2 * x[i] * K2 * cosAB - K2;
			double mqmq = (m1 * q2 - m2 * q1);
			
			y0 += 0;
			if (false) // (fabs(mqmq) > 1e-10) // August/2018
			{
				double c0 = a * (p2*q1 - p1*q2) / mqmq;

				if (c0 > 0)
				{
					if (angC == angAB && c0 < 1e-10) continue;

					interpretationsX3[ang][numInterpretationsX3[ang]][0] = a;
					interpretationsX3[ang][numInterpretationsX3[ang]][1] = b;
					interpretationsX3[ang][numInterpretationsX3[ang]][2] = c0;
					numInterpretationsX3[ang]++;
				}
			}
			else
			{
				double c1 = a * (cosCA + sqrt(y0)); // A27
				double c2 = a * (cosCA - sqrt(y0)); // A27
				double criterion1 = fabs(b*b + c1*c1 - 2 * b*c1*cosBC - Rbc*Rbc); // A3
				double criterion2 = fabs(b*b + c2*c2 - 2 * b*c2*cosBC - Rbc*Rbc); // A3

				if (c1 > 0 && (angC == angAB && c1 < 1e-10)==false && criterion1 <= 1e-4*Rab)
				{
					//if (angC == angAB && c1 < 1e-10) continue;
					//if (fabs(b*b + c1*c1 - 2 * b*c1*cosBC - Rbc*Rbc) > 1e-10) continue; // A3

					interpretationsX3[ang][numInterpretationsX3[ang]][0] = a;
					interpretationsX3[ang][numInterpretationsX3[ang]][1] = b;
					interpretationsX3[ang][numInterpretationsX3[ang]][2] = c1;
					numInterpretationsX3[ang]++;
				}

				if (c2 > 0 && (angC == angAB && c2 < 1e-10)==false && criterion2 <= 1e-4*Rab)
				{
					//if (angC == angAB && c2 < 1e-10) continue;
					//if (fabs(b*b + c2*c2 - 2 * b*c2*cosBC - Rbc*Rbc) > 1e-10) continue; // A3

					interpretationsX3[ang][numInterpretationsX3[ang]][0] = a;
					interpretationsX3[ang][numInterpretationsX3[ang]][1] = b;
					interpretationsX3[ang][numInterpretationsX3[ang]][2] = c2;
					numInterpretationsX3[ang]++;
				}
			}
		} // for (int i = 0; i < numPolySolutions; i++)

		// Check redundancy
		int check_redundancy[8] = {};
		double diff;
		for (int i = 0; i < numInterpretationsX3[ang] - 1; i++)
		{
			if (check_redundancy[i] == 1) continue;
			for (int j = i+1; j < numInterpretationsX3[ang]; j++)
			{
				diff  = pow(interpretationsX3[ang][i][0] - interpretationsX3[ang][j][0], 2.0);
				diff += pow(interpretationsX3[ang][i][1] - interpretationsX3[ang][j][1], 2.0);
				diff += pow(interpretationsX3[ang][i][2] - interpretationsX3[ang][j][2], 2.0);

				if (diff < 1e-4*Rab)
					check_redundancy[j] = 1;
			}
		}
		int k = 0;
		for (int i = 0; i < numInterpretationsX3[ang]; i++)
		{
			if (check_redundancy[i] == 1)
				continue;
			else
			{
				interpretationsX3[ang][k][0] = interpretationsX3[ang][i][0];
				interpretationsX3[ang][k][1] = interpretationsX3[ang][i][1];
				interpretationsX3[ang][k][2] = interpretationsX3[ang][i][2];
				k++;
			}
		}
		numInterpretationsX3[ang] = k;

		// Examine results
		double la, lb, lc;
		double reangV, reangT1, reangT2;
		int validity[8] = {};
		double check3D1a, check3D1b, check3D2;
		for (int i = 0; i < numInterpretationsX3[ang]; i++)
		{
			la = interpretationsX3[ang][i][0];
			lb = interpretationsX3[ang][i][1];
			lc = interpretationsX3[ang][i][2];

			// (a,b)
			reangV  = acos( (la*la + lb*lb - Rab*Rab)/(2*la*lb) )  * 180.0/pi;
			reangT1 = acos( (Rab*Rab + la*la - lb*lb)/(2*Rab*la) ) * 180.0/pi;
			reangT2 = acos( (Rab*Rab + lb*lb - la*la)/(2*Rab*lb) ) * 180.0/pi;
			check3D1a = fabs(reangV - angAB);
			check3D1b = angAB;
			check3D2 = 180 - reangV - reangT1 - reangT2;
			if (fabs(reangV-angAB) > 0.01 && fabs(reangV - angAB) > angAB*0.01)
			{
				//											printf("Err1: %f, %f, %f\n", check3D1a, check3D1b, check3D2);
				validity[i] = 1;
			}
			//if (fabs(180 - reangV-reangT1-reangT2) > 0.01)
			//{
			//	printf("Err2: %f, %f, %f\n", check3D1a, check3D1b, check3D2);
			//	validity[i] = 1;
			//}
			
			// (b,c)
			reangV  = acos( (lb*lb + lc*lc - Rbc*Rbc)/(2*lb*lc) )  * 180.0/pi;
			reangT1 = acos( (Rbc*Rbc + lb*lb - lc*lc)/(2*Rbc*lb) ) * 180.0/pi;
			reangT2 = acos( (Rbc*Rbc + lc*lc - lb*lb)/(2*Rbc*lc) ) * 180.0/pi;
			check3D1a = fabs(reangV - angBC);
			check3D1b = angBC;
			check3D2 = 180 - reangV - reangT1 - reangT2;
			if (fabs(reangV-angBC) > 0.01 && fabs(reangV - angBC) > angBC*0.01)
			{
				//											printf("Err1: %f, %f, %f\n", check3D1a, check3D1b, check3D2);
				validity[i] = 1;
			}
			//if (fabs(180 - reangV - reangT1 - reangT2) > 0.01)
			//{
			//	printf("Err2: %f, %f, %f\n", check3D1a, check3D1b, check3D2);
			//	validity[i] = 1;
			//}
			
			// (c,a)
			reangV  = acos( (lc*lc + la*la - Rca*Rca)/(2*lc*la) )  * 180.0/pi;
			reangT1 = acos( (Rca*Rca + lc*lc - la*la)/(2*Rca*lc) ) * 180.0/pi;
			reangT2 = acos( (Rca*Rca + la*la - lc*lc)/(2*Rca*la) ) * 180.0/pi;
			check3D1a = fabs(reangV - angCA);
			check3D1b = angCA;
			check3D2 = 180 - reangV - reangT1 - reangT2;
			if (fabs(reangV-angCA) > 0.01 && fabs(reangV - angCA) > angCA*0.01)
			{
				//											printf("Err1: %f, %f, %f\n", check3D1a, check3D1b, check3D2);
				validity[i] = 1;
			}
			//if (fabs(180 - reangV-reangT1-reangT2) > 0.01)
			//{
			//	printf("Err2: %f, %f, %f\n", check3D1a, check3D1b, check3D2);
			//	validity[i] = 1;
			//}
		}
		k = 0;
		for (int i = 0; i < numInterpretationsX3[ang]; i++)
		{
			if (validity[i] == 1)
				continue;
			interpretationsX3[ang][k][0] = interpretationsX3[ang][i][0];
			interpretationsX3[ang][k][1] = interpretationsX3[ang][i][1];
			interpretationsX3[ang][k][2] = interpretationsX3[ang][i][2];
			k++;
		}
		numInterpretationsX3[ang] = k;
		if (4 < numInterpretationsX3[ang])
		{
			printf("Err3 (4<num): %f, %f, %f\n", check3D1a, check3D1b, check3D2);
			numInterpretationsX3[ang] = 4;
		}


		// Sort results
		double vj0, vj1, v_swap[3];
		for (int i = 0; i < numInterpretationsX3[ang]-1; i++)
			for (int j = numInterpretationsX3[ang]-1; j>i; j--)
			{
				vj0 = interpretationsX3[ang][j  ][0]*100 + interpretationsX3[ang][j  ][1]*10 + interpretationsX3[ang][j  ][2];
				vj1 = interpretationsX3[ang][j-1][0]*100 + interpretationsX3[ang][j-1][1]*10 + interpretationsX3[ang][j-1][2];

				if (vj0<vj1)
				{
					for (int k = 0; k < 3; k++) v_swap[k] = interpretationsX3[ang][j][k];
					for (int k = 0; k < 3; k++) interpretationsX3[ang][j][k] = interpretationsX3[ang][j-1][k];
					for (int k = 0; k < 3; k++) interpretationsX3[ang][j-1][k] = v_swap[k];
				}
			}

		for (int i = 0; i < numInterpretationsX3[ang]; i++)
		{
			double temp_interpretations[3];
			if (ang == 0)
			{
				temp_interpretations[0] = interpretationsX3[ang][i][0] / Rab;
				temp_interpretations[1] = interpretationsX3[ang][i][1] / Rab;
				temp_interpretations[2] = interpretationsX3[ang][i][2] / Rab;
			}
			else if(ang == 1)
			{
				temp_interpretations[0] = interpretationsX3[ang][i][2] / Rca;
				temp_interpretations[1] = interpretationsX3[ang][i][0] / Rca;
				temp_interpretations[2] = interpretationsX3[ang][i][1] / Rca;
			}
			else //ang == 2
			{
				temp_interpretations[0] = interpretationsX3[ang][i][1] / Rbc;
				temp_interpretations[1] = interpretationsX3[ang][i][2] / Rbc;
				temp_interpretations[2] = interpretationsX3[ang][i][0] / Rbc;
			}
			for (int k = 0; k < 3; k++)
				interpretationsX3[ang][i][k] = temp_interpretations[k];
		}
	} // for (int ang = 0; ang < 3; ang++)


	// Copy data
	int numInterpretations = 0;
	int check_consistency[3][4] = {};
	for (int ang = 0; ang < 3; ang++) for (int i = 0; i < numInterpretationsX3[ang]; i++) check_consistency[ang][i] = 1;
	//ang = 0
	for (int i = 0; i<numInterpretationsX3[0]; i++)
	{
		for (int abc = 0; abc < 3; abc++) interpretations[i][abc] = interpretationsX3[0][i][abc];
		check_consistency[0][i] = 0;
		numInterpretations++;

		//ang = 0 vs ang = 1+2
		for (int ang = 1; ang < 3; ang++) for (int j = 0; j<numInterpretationsX3[ang]; j++)
		{
			double diff;
			diff  = pow(interpretationsX3[0][i][0] - interpretationsX3[ang][j][0], 2.0);
			diff += pow(interpretationsX3[0][i][1] - interpretationsX3[ang][j][1], 2.0);
			diff += pow(interpretationsX3[0][i][2] - interpretationsX3[ang][j][2], 2.0);
			if (diff < 1e-2)
				check_consistency[ang][j] = 0;
		}
	}
	//ang = 1
	for (int j = 0; j<numInterpretationsX3[1]; j++)
	{
		if (check_consistency[1][j] == 0)
			continue;
		for (int abc = 0; abc < 3; abc++) interpretations[numInterpretations][abc] = interpretationsX3[1][j][abc];
		check_consistency[1][j] = 0;
		numInterpretations++;

		//ang = 1 vs ang = 2
		for (int k = 0; k<numInterpretationsX3[2]; k++)
		{
			double diff;
			diff  = pow(interpretationsX3[1][j][0] - interpretationsX3[2][k][0], 2.0);
			diff += pow(interpretationsX3[1][j][1] - interpretationsX3[2][k][1], 2.0);
			diff += pow(interpretationsX3[1][j][2] - interpretationsX3[2][k][2], 2.0);
			if (diff < 1e-2)
				check_consistency[2][k] = 0;
		}
	}
	for (int k = 0; k < numInterpretationsX3[2]; k++)
	{
		if (check_consistency[2][k] == 0)
			continue;
		for (int abc = 0; abc < 3; abc++) interpretations[numInterpretations][abc] = interpretationsX3[2][k][abc];
		numInterpretations++;
	}

	int total_consistency = 0;
	for (int ang = 0; ang < 3; ang++) for (int i = 0; i < 4; i++) total_consistency += check_consistency[ang][i];
	if (total_consistency > 0)
		total_consistency = 0;

	return numInterpretations;
}

int SrereoEyeRecovery_Optimization(double *eyeL, double *eyeR, int numPoints, double *px, double *py, double *pz, double noise, double ang_condition)
{

#define MaxA 900 //1800
#define MaxB 900 //1800

	if (numPoints > MaxNp) return -1;

	char x = 0, y = 1, z = 2;
	char a = 0, b = 1, c = 2;
	char l = 0, r = 1;

	double p[MaxNp][3];
	for (int i = 0; i < numPoints; i++)
	{
		p[i][x] = px[i];
		p[i][y] = py[i];
		p[i][z] = pz[i];
	}

	//Avoid a case seeing different faces of the triangle from two eyes for simplicity
	if(true)
	{
		double temp_vl_a[3], temp_vr_a[3];
		double temp_vl_b[3], temp_vr_b[3];
		double temp_vl_c[3], temp_vr_c[3];
		TdS_Vector3D_Minus(temp_vl_a, eyeL, &(p[a][x]), 1);
		TdS_Vector3D_Minus(temp_vr_a, eyeR, &(p[a][x]), 1);
		TdS_Vector3D_Minus(temp_vl_b, eyeL, &(p[b][x]), 1);
		TdS_Vector3D_Minus(temp_vr_b, eyeR, &(p[b][x]), 1);
		TdS_Vector3D_Minus(temp_vl_c, eyeL, &(p[c][x]), 1);
		TdS_Vector3D_Minus(temp_vr_c, eyeR, &(p[c][x]), 1);

		double vector_ab[3], vector_bc[3], vector_ca[3];
		TdS_Vector3D_Minus(vector_ab, &(p[b][x]), &(p[a][x]), 0); // b-a; a->b
		TdS_Vector3D_Minus(vector_bc, &(p[c][x]), &(p[b][x]), 0); // c-b; b->c
		TdS_Vector3D_Minus(vector_ca, &(p[a][x]), &(p[c][x]), 0); // a-c; c->a

		double cross_a[3], cross_b[3], cross_c[3];
		TdS_Vector3D_CrossProduct(cross_a, vector_ca, vector_ab);
		TdS_Vector3D_CrossProduct(cross_b, vector_ab, vector_bc);
		TdS_Vector3D_CrossProduct(cross_c, vector_bc, vector_ca);

		if (TdS_Vector3D_DotProduct(cross_a, temp_vl_a)*TdS_Vector3D_DotProduct(cross_a, temp_vr_a) < 0) return -2;
		if (TdS_Vector3D_DotProduct(cross_b, temp_vl_b)*TdS_Vector3D_DotProduct(cross_b, temp_vr_b) < 0) return -2;
		if (TdS_Vector3D_DotProduct(cross_c, temp_vl_c)*TdS_Vector3D_DotProduct(cross_c, temp_vr_c) < 0) return -2;
	}


	double ang[3][MaxNp][2] = {};
	double interoculardistance = TdS_Point3D_Distance(eyeL, eyeR);

	double distance_veridical[MaxNp][2] = {};
	double vergence_veridical[MaxNp] = {};

	for (int i = 0; i < numPoints; i++)
	{
		ang[a][i][l] = TdS_Point3D_Angle(eyeL, &(p[a][x]), &(p[i][x]), NULL, 1); //deg
		ang[a][i][r] = TdS_Point3D_Angle(eyeR, &(p[a][x]), &(p[i][x]), NULL, 1); //deg

		ang[b][i][l] = TdS_Point3D_Angle(eyeL, &(p[b][x]), &(p[i][x]), NULL, 1); //deg
		ang[b][i][r] = TdS_Point3D_Angle(eyeR, &(p[b][x]), &(p[i][x]), NULL, 1); //deg

		ang[c][i][l] = TdS_Point3D_Angle(eyeL, &(p[c][x]), &(p[i][x]), NULL, 1); //deg
		ang[c][i][r] = TdS_Point3D_Angle(eyeR, &(p[c][x]), &(p[i][x]), NULL, 1); //deg

		distance_veridical[i][l] = TdS_Point3D_Distance(eyeL, &(p[i][x])) / interoculardistance;
		distance_veridical[i][r] = TdS_Point3D_Distance(eyeR, &(p[i][x])) / interoculardistance;

		double temp_vl[3], temp_vr[3];
		TdS_Vector3D_Minus(temp_vl, eyeL, &(p[i][x]), 1);
		TdS_Vector3D_Minus(temp_vr, eyeR, &(p[i][x]), 1);
		vergence_veridical[i] = TdS_Vector3D_Angle(temp_vl, temp_vr, NULL, 1);
	}

	double angMax = 0;
	for (int abc = 0; abc < 3; abc++) for (int i = 0; i < numPoints; i++)
	{
		if (abc == i) continue;
		if (ang[abc][i][l] < ang_condition) return -2;
		if (ang[abc][i][r] < ang_condition) return -2;

		if (ang[abc][i][l] > angMax) angMax = ang[abc][i][l];
		if (ang[abc][i][r] > angMax) angMax = ang[abc][i][r];
	}

	double vertangle[MaxNp] = {};
	for (int i = 0; i < numPoints; i++)
		vertangle[i] = atan(p[i][y] / p[i][z]) / Pi * 180.0;
	TdS_BubbleSort(vertangle, numPoints);
	for (int i = 1; i < numPoints; i++)
		if (vertangle[i] - vertangle[i - 1] < VertDiff) return -2;


	//Random noise
	if (noise > 0)
	{
		std::random_device rd;
		std::mt19937 e2(rd());
		std::normal_distribution<> dist(1, noise);
		for (int i = 0; i < numPoints; i++)
		{

			if (CORRELATED_NOISE)
			{
				double temp_noise;
				temp_noise = fabs(dist(e2));// TdS_RandomNumber(1-noise, 1+noise);
				ang[a][i][l] += temp_noise;
				ang[a][i][r] += temp_noise;

				temp_noise = fabs(dist(e2));// TdS_RandomNumber(1-noise, 1+noise);
				ang[b][i][l] += temp_noise;
				ang[b][i][r] += temp_noise;

				temp_noise = fabs(dist(e2));// TdS_RandomNumber(1-noise, 1+noise);
				ang[c][i][l] += temp_noise;
				ang[c][i][r] += temp_noise;
			}
			else
			{
				ang[a][i][l] += fabs(dist(e2));// TdS_RandomNumber(1-noise, 1+noise);
				ang[a][i][r] += fabs(dist(e2));// TdS_RandomNumber(1-noise, 1+noise);

				ang[b][i][l] += fabs(dist(e2));// TdS_RandomNumber(1-noise, 1+noise);
				ang[b][i][r] += fabs(dist(e2));// TdS_RandomNumber(1-noise, 1+noise);

				ang[c][i][l] += fabs(dist(e2));// TdS_RandomNumber(1-noise, 1+noise);
				ang[c][i][r] += fabs(dist(e2));// TdS_RandomNumber(1-noise, 1+noise);
			}
		}
	}



	double angA = TdS_Point3D_Angle(&(p[a][x]), &(p[b][x]), &(p[c][x]), NULL, 1); //deg
	double angB = TdS_Point3D_Angle(&(p[b][x]), &(p[c][x]), &(p[a][x]), NULL, 1); //deg
	double angC = TdS_Point3D_Angle(&(p[c][x]), &(p[a][x]), &(p[b][x]), NULL, 1); //deg



	double diff_dist, diff_verg;
	double debug_cost_min = SrereoEyeRecovery_Cost(angA, angB, ang, numPoints, distance_veridical, vergence_veridical, &diff_dist, &diff_verg);


	// Minimum angle is 1
	if (angA < 1 || angB < 1 || angC < 1) return -2;
	for (int i = 0; i < numPoints; i++) for (int abc = 0; abc < 3; abc++)
	{
		if (abc == i) continue;
		if (ang[abc][i][l] < 1) return -3;
		if (ang[abc][i][r] < 1) return -3;
	}
	

	// Exhaustive test: 180*180
	double degA_min = -1, degB_min = -1;
	double cost_min = 999999;
	ofstream fout;
	float zero = 0, minus1 = -1;
	
	if (SAVE_TABLE) //val_degA = 0;
	{
		fout.open("Output_table.bin", ios::out | ios::binary | ios::trunc);
		for (int degB = 0; degB <= MaxB; degB++) fout.write((char*)&minus1, sizeof(float));
	}

	//0 < val_degA < 180;
	for (int degA = 1; degA < MaxA; degA++)
	{
		float costArray[MaxB+1] = {};
		for (int degB = 0; degB <= MaxB; degB++) costArray[degB] = -1;

		for (int degB = 0; degB <= MaxB; degB++)
		{
			double val_degA = degA*180.0 / MaxA;
			double val_degB = degB*180.0 / MaxB;

			if (val_degA == 0 || val_degA >= 180 || val_degB == 0 || val_degB >= 180) continue;
			if (180 < val_degA + val_degB) continue;

			costArray[degB] = (float)(SrereoEyeRecovery_Cost(val_degA, val_degB, ang, numPoints, distance_veridical, vergence_veridical, NULL, NULL));

			if (costArray[degB] > 0 && costArray[degB] < cost_min)
			{
				degA_min = val_degA;
				degB_min = val_degB;
				cost_min = costArray[degB];
			}
		} // degB

		if (SAVE_TABLE)
		{
			for (int degB = 0; degB <= MaxB; degB++)
				fout.write((char*)&(costArray[degB]), sizeof(float));
		}

	} // degA

	if (SAVE_TABLE) //val_degA = 180;
	{
		for (int degB = 0; degB <= MaxB; degB++) fout.write((char*)&minus1, sizeof(float));
		fout.close();
	}
	

  // Compare recovered 3D scene with veridical 3D scene

	if (degA_min <= 0 || degB_min <= 0)
		return -4;

	double diff_ang, diff_cost;
	double temp_cost_min = SrereoEyeRecovery_Cost(degA_min, degB_min, ang, numPoints, distance_veridical, vergence_veridical, &diff_dist, &diff_verg);
	//double debug_cost_min = SrereoEyeRecovery_Cost(angA, angB, ang, numPoints, distance_veridical, vergence_veridical, &diff_dist, &diff_verg);
	diff_ang = sqrt( pow(angA-degA_min,2.0) + pow(angB-degB_min,2.0) );
	diff_cost = temp_cost_min - debug_cost_min;
	

	if (SAVE_TABLE)
	{
		printf("Min cost = %f at (%f,%f)\n", (float)temp_cost_min, (float)degA_min, (float)degB_min);
		printf("Min cost = %f at (%f,%f)\n", (float)debug_cost_min, (float)angA, (float)angB);
		printf("    Distance error = %f\n", (float)(diff_dist));
		printf("    Vergence error = %f\n", (float)(diff_verg));
		printf("    Angle    error = %f\n", (float)(diff_ang));
	}

	//printf("\n");
	//printf("Skew cost (%d,%d) = %f\n", (int)(degA_min), (int)(degB_min), temp_cost_min);
	//printf("Angle distortion (%d,%d) = %f\n", (int)(degA_min), (int)(degB_min), diff_veridical_recovered);
	//printf("\n");
	if (SAVE_SIMULATION)
	{
		//fout.open("Output_data.txt", ios::out | ios::app);
		fout.open(fout_name, ios::out | ios::app);
		fout << temp_cost_min << "\t";
		fout << diff_ang << "\t" << diff_dist << "\t" << diff_verg << "\t";
		fout << "(" << angA << "," << angB << ")" << "\t";
		fout << "(" << degA_min << "," << degB_min << ")" << endl;
		fout.close();
	}

	if(false) if (diff_ang > 30)
	{
		fout.open("LargeError.txt", ios::out | ios::app);
		fout << "diff_ang = " << diff_ang << endl;
		fout << "Cost:" << debug_cost_min << ", " << temp_cost_min << endl;
		fout << "(" << angA << "," << angB << ")" << ", (" << degA_min << "," << degB_min << ")" << endl;
		fout << "double px[MaxNp] = { ";
		for (int i = 0; i < numPoints-1; i++)
			fout << px[i] << ", ";
		fout << px[numPoints - 1] << " };" << endl;

		fout << "double py[MaxNp] = { ";
		for (int i = 0; i < numPoints-1; i++)
			fout << py[i] << ", ";
		fout << py[numPoints - 1] << " };" << endl;

		fout << "double pz[MaxNp] = { ";
		for (int i = 0; i < numPoints-1; i++)
			fout << pz[i] << ", ";
		fout << pz[numPoints - 1] << " };" << endl;
		fout << endl;
		fout.close();
	}
	return 0;
}

double SrereoEyeRecovery_Cost(	double degA, double degB, double(&visual_angles)[3][MaxNp][2], int numPoints,	// Step 1: Recover 3D scene & Compute cost
	double(&distance_veridical)[MaxNp][2], double(&vergence_veridical)[MaxNp], double *diff_dist, double *diff_verg)	// Step 2: Compare recovered 3D scene with veridical 3D scene
{
  // Step 1: Recover 3D scene & Compute cost

	if (numPoints > MaxNp) return -1;
	if (degA == 0 || degA >= 180 || degB == 0 || degB >= 180) return -1;
	if (180 < degA + degB) return -1;

	char x = 0, y = 1, z = 2;
	char a = 0, b = 1, c = 2;
	char l = 0, r = 1;
	double lAB = 1;
	double degC = 180 - degA - degB;
	double distA, distB, distC;
	double projection_vectors[MaxNp][2][4][3] = {}; //[MaxNp][LR][#Interpretations][XYZ]
	int num_interpretations[2];
	double interpretations_leye[10][3] = {}; //[#Interpretatins][ABC]
	double interpretations_reye[10][3] = {}; //[#Interpretatins][ABC]

	// Left eye 1/2
	num_interpretations[l] = P3P_MinkovSawada(degA, degB, visual_angles[a][b][l], visual_angles[b][c][l], visual_angles[c][a][l], interpretations_leye);
	if (num_interpretations[l] <= 0)
		return -1;
	// Right eye 1/2
	num_interpretations[r] = P3P_MinkovSawada(degA, degB, visual_angles[a][b][r], visual_angles[b][c][r], visual_angles[c][a][r], interpretations_reye);
	if (num_interpretations[r] <= 0)
		return -1;

	double pointA[3], pointB[3], pointC[3]; //[XYZ]
	double eyes_recovered[2][4][3] = {}; //[LR][#Interpretations][XYZ]
	int check[2][4] = {};

	// Left eye 2/2
	for (int il = 0; il < num_interpretations[l]; il++)
	{
		distA = interpretations_leye[il][a];
		distB = interpretations_leye[il][b];
		distC = interpretations_leye[il][c];
		check[l][il] = TdS_Tetrahedron3D_CartesianCoordinates(&(eyes_recovered[l][il][x]), pointA, pointB, pointC, degA, degB, degC,
																	visual_angles[a][b][l],
																	visual_angles[b][c][l],
																	visual_angles[c][a][l], lAB, distA, distB, distC);
		if (check[l][il] < 0)
			continue; // il
		TdS_Vector3D_Minus(&(projection_vectors[a][l][il][x]), &(eyes_recovered[l][il][x]), pointA, 1);
		TdS_Vector3D_Minus(&(projection_vectors[b][l][il][x]), &(eyes_recovered[l][il][x]), pointB, 1);
		TdS_Vector3D_Minus(&(projection_vectors[c][l][il][x]), &(eyes_recovered[l][il][x]), pointC, 1);
		for (int j = 3; j < numPoints; j++)
		{
			if (0 > ComputeAnotherProjectionLine(&(projection_vectors[j][l][il][x]), &(eyes_recovered[l][il][x]), pointA, pointB, pointC,
													visual_angles[a][j][l],
													visual_angles[b][j][l],
													visual_angles[c][j][l], 1))
			{
				check[l][il] = -1;
				break; // j (not il)
			}
		}// j
	}// il

	// Right eye 2/2
	for (int ir = 0; ir < num_interpretations[r]; ir++)
	{
		distA = interpretations_reye[ir][a];
		distB = interpretations_reye[ir][b];
		distC = interpretations_reye[ir][c];
		check[r][ir] = TdS_Tetrahedron3D_CartesianCoordinates(&(eyes_recovered[r][ir][x]), pointA, pointB, pointC, degA, degB, degC,
																	visual_angles[a][b][r],
																	visual_angles[b][c][r],
																	visual_angles[c][a][r], lAB, distA, distB, distC);
		if (check[r][ir] < 0)
			continue; // ir
		TdS_Vector3D_Minus(&(projection_vectors[a][r][ir][x]), &(eyes_recovered[r][ir][x]), pointA, 1);
		TdS_Vector3D_Minus(&(projection_vectors[b][r][ir][x]), &(eyes_recovered[r][ir][x]), pointB, 1);
		TdS_Vector3D_Minus(&(projection_vectors[c][r][ir][x]), &(eyes_recovered[r][ir][x]), pointC, 1);
		for (int j = 3; j < numPoints; j++)
		{
			if (0 > ComputeAnotherProjectionLine(&(projection_vectors[j][r][ir][x]), &(eyes_recovered[r][ir][x]), pointA, pointB, pointC,
													visual_angles[a][j][r],
													visual_angles[b][j][r],
													visual_angles[c][j][r], 1))
			{
				check[r][ir] = -1;
				break; // j (not ir)
			}
		}// j
	}// ir

	// Pairs of eyes
	double point_from_L[3], point_from_R[3];
	double dist_from_L, dist_from_R;
	int i_pair = 0;
	double cost4x4[16] = {};
	int pairs_ilir[16][2] = {};

	for (int il = 0; il < num_interpretations[l]; il++)
	{
		if (check[l][il] < 0) continue; // il

		for (int ir = 0; ir < num_interpretations[r]; ir++)
		{
			if (check[r][ir] < 0) continue; // ir
			pairs_ilir[i_pair][l] = il;
			pairs_ilir[i_pair][r] = ir;

			cost4x4[i_pair] = 0;
			for (int j = 3; j < numPoints; j++)
			{
				double error_intersection;
				error_intersection = TdS_Line3D_Intersection(point_from_L, point_from_R, &dist_from_L, &dist_from_R,
					&(projection_vectors[j][l][il][x]), &(eyes_recovered[l][il][x]),
					&(projection_vectors[j][r][ir][x]), &(eyes_recovered[r][ir][x]));
				if (error_intersection < 0 || dist_from_L < EPS || dist_from_R < EPS)
				{
					cost4x4[i_pair] = 999999.1;
					goto EndCheck;// next ir
				}
				else
					cost4x4[i_pair] += error_intersection / sqrt(fabs(dist_from_L) * fabs(dist_from_R));
			} // j

			if (cost4x4[i_pair] >= 999998) goto EndCheck;// next ir
			
			double ioaxis[3];
			TdS_Vector3D_Minus(ioaxis, &(eyes_recovered[l][il][x]), &(eyes_recovered[r][ir][x]), 1);	
			for (int j1 = 0; j1 < numPoints; j1++)for (int j2 = j1 + 1; j2 < numPoints; j2++)for (int j3 = j2 + 1; j3 < numPoints; j3++)
			{
				double weights[3];
				if (0 < TdS_SimultaneousEquations3(weights,
						projection_vectors[j1][l][il][x], projection_vectors[j2][l][il][x], projection_vectors[j3][l][il][x], -ioaxis[x],
						projection_vectors[j1][l][il][y], projection_vectors[j2][l][il][y], projection_vectors[j3][l][il][y], -ioaxis[y],
						projection_vectors[j1][l][il][z], projection_vectors[j2][l][il][z], projection_vectors[j3][l][il][z], -ioaxis[z]))
				{
					if ( (weights[x] > 0 && weights[y] > 0 && weights[z] > 0) || (weights[x] < 0 && weights[y] < 0 && weights[z] < 0) )
						cost4x4[i_pair] = 999999.3;
				} // if (0 < TdS_SimultaneousEquations3(weights,...

				if (cost4x4[i_pair] >= 999998) goto EndCheck;// j1,j2,j3

				if (0 > TdS_SimultaneousEquations3(weights,
						projection_vectors[j1][r][ir][x], projection_vectors[j2][r][ir][x], projection_vectors[j3][r][ir][x], -ioaxis[x],
						projection_vectors[j1][r][ir][y], projection_vectors[j2][r][ir][y], projection_vectors[j3][r][ir][y], -ioaxis[y],
						projection_vectors[j1][r][ir][z], projection_vectors[j2][r][ir][z], projection_vectors[j3][r][ir][z], -ioaxis[z]))
				{
					if ((weights[x] > 0 && weights[y] > 0 && weights[z] > 0) || (weights[x] < 0 && weights[y] < 0 && weights[z] < 0))
						cost4x4[i_pair] = 999999.5;
				} // if (0 > TdS_SimultaneousEquations3(weights,...

				if (cost4x4[i_pair] >= 999998) goto EndCheck;// j1,j2,j3
			}

		  EndCheck:
			i_pair++;
		} // ir
	}// il

	int i_pair_min;
	double cost_min = TdS_Min(i_pair, cost4x4, &i_pair_min);
	if (cost_min < 0 || cost_min >= 999998)
		return -1;

	if(diff_dist == NULL && diff_verg == NULL)
		return cost_min;

  // Step 2: Compare recovered 3D scene with veridical 3D scene
	int il_min = pairs_ilir[i_pair_min][l];
	int ir_min = pairs_ilir[i_pair_min][r];
	
	if (diff_dist != NULL)
	{
		double interocular_dist = TdS_Point3D_Distance(&(eyes_recovered[l][il_min][x]), &(eyes_recovered[r][ir_min][x]));
		double dist_err = 0;
		for (int j = 0; j < numPoints; j++)
		{
			TdS_Line3D_Intersection(point_from_L, point_from_R, &dist_from_L, &dist_from_R,
				&(projection_vectors[j][l][il_min][x]), &(eyes_recovered[l][il_min][x]),
				&(projection_vectors[j][r][ir_min][x]), &(eyes_recovered[r][ir_min][x]));

			dist_err += fabs((dist_from_L / interocular_dist - distance_veridical[j][l]) / distance_veridical[j][l]);
			dist_err += fabs((dist_from_R / interocular_dist - distance_veridical[j][r]) / distance_veridical[j][r]);
		}
		*diff_dist = dist_err / (double)(numPoints);
	}

	if (diff_verg != NULL)
	{
		double temp_ang;
		double verg_err = 0;
		for (int j = 0; j < numPoints; j++)
		{
			temp_ang = TdS_Vector3D_Angle(	&(projection_vectors[j][l][il_min][x]),
											&(projection_vectors[j][r][ir_min][x]), NULL, 1);
			verg_err += fabs(temp_ang - vergence_veridical[j]);
		}
		*diff_verg = verg_err / (double)(numPoints);
	}

	return cost_min;
}

int ComputeAnotherProjectionLine(double *vq, double *eye, double *pa, double *pb, double *pc, double anga, double angb, double angc, bool raddeg)
{
	double cosa, cosb, cosc;
	double va[3], vb[3], vc[3];
	double la, lb, lc;
	for (int xyz = 0; xyz < 3; xyz++)
	{
		va[xyz] = pa[xyz] - eye[xyz];
		vb[xyz] = pb[xyz] - eye[xyz];
		vc[xyz] = pc[xyz] - eye[xyz];
	}
	la = TdS_Vector3D_Length(va);
	lb = TdS_Vector3D_Length(vb);
	lc = TdS_Vector3D_Length(vc);
	if (la < EPS || lb < EPS || lc < EPS) return -1;// error

	TdS_Vector3D_Scale(va, 1.0 / la);
	TdS_Vector3D_Scale(vb, 1.0 / lb);
	TdS_Vector3D_Scale(vc, 1.0 / lc);

	if (raddeg)
	{
		cosa = cos(anga*Pi / 180.0);
		cosb = cos(angb*Pi / 180.0);
		cosc = cos(angc*Pi / 180.0);
	}
	else
	{
		cosa = cos(anga);
		cosb = cos(angb);
		cosc = cos(angc);
	}

	if (0 > TdS_SimultaneousEquations3(vq, va[0], va[1], va[2], -cosa, vb[0], vb[1], vb[2], -cosb, vc[0], vc[1], vc[2], -cosc))
		return -2;//error

	TdS_Vector3D_Normalize(vq);
	//double vqva = TdS_Vector3D_DotProduct(vq, va);
	//double vqvb = TdS_Vector3D_DotProduct(vq, vb);
	//double vqvc = TdS_Vector3D_DotProduct(vq, vc);
	//double consist0 = fabs(vqva-cosa)+fabs(vqvb-cosb)+fabs(vqvc-cosc);
	//double consist1 = fabs(vqva+cosa)+fabs(vqvb+cosb)+fabs(vqvc+cosc);
	//if (consist1 < consist0) TdS_Vector3D_Scale(vq, -1.0);

	return 1;
}

int SolvePolyEquation(double *x, double coef4, double coef3, double coef2, double coef1, double coef0)
{
	int res = -1;
	if (fabs(coef4) < 1e-10) // 3rd-polynomial or less
	{
		if (fabs(coef3) < 1e-10) // 2nd-polynomial or less
		{
			if (fabs(coef2) < 1e-10) // 1st-polynomial or less
			{
				if (fabs(coef1) < 1e-10)
					return -1; // error

							   // 1st-polynomial
				x[0] = -coef0 / coef1;
				return 1;
			}

			// 2nd-polynomial
			double b2 = coef1 / coef2;
			double c2 = coef0 / coef2;
			double determinant = b2 * b2 - 4 * c2;
			if (fabs(determinant) < 1e-10)
			{
				x[0] = -b2 / 2;
				return 1;
			}
			else if (determinant > 0)
			{
				x[0] = (-b2 - sqrt(determinant)) / 2;
				x[1] = (-b2 + sqrt(determinant)) / 2;
				return 2;
			}
			else
				return 0;
		}

		// 3rd-polynomial
		res = SolveP3_ver2(x, coef2 / coef3, coef1 / coef3, coef0 / coef3);
		TdS_BubbleSort(x, res);
	}
	else
		// 4th-polynomial
		res = SolveP4_ver2(x, coef3 / coef4, coef2 / coef4, coef1 / coef4, coef0 / coef4);

	// Remove identical roots of 3rd- and 4th-polynomials
	int i, j;
	if (res > 1)
	{
		for (i = 1, j = 1; i < res; i++)
			if (fabs(x[i - 1] - x[i]) > EPS) x[j++] = x[i];
		res = j;
	}

	return res;
}
