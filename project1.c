/* Project 1 - Microcanonical Ensemble Simulation
 * project1.c:  Main program
 *
 * gcc -Wall -lm -std=c99 project1.c
 *
 * Ben Niehoff
 * Jose Gonzalez
 *
 * Description:
 *
 * This program simulates a system of Argon atoms by explicitly
 * integrating the equations of motion of a collection of
 * N atoms.  By measuring energy fluctuations, the program determines
 * the specific heat and pressure of the system.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "project1.h"

int main() {	
	time_t Seed;
	//char dummy[3];
	
	int State;
	int last_rescale;
	int ActualIterations;
	
	clock_t Start, End;
	double Time;
		
	Seed = time(NULL);
	//Seed = 1226117015;
	srand(Seed);
	
	printf(":  Molecular Dynamics Simulation\n");
	printf(":  Running with random seed %ld\n", (long) Seed);
	printf(":...........................................\n\n");
	
	CellCount = 14;
	rho_s = 0.8;
	T_s = 0.9;
	
	printf("Settings:\n");
	printf("Cell Count:  %d,  rho_s:  %f,  T_s:  %f\n\n", CellCount, rho_s, T_s);
	
	/* Set up initial system */
	Initialize();
	
	//PositionCheck();
	
	//printf("Before adjusting P, L, and T:\n");
	//Diagnostic();
	
	CenterOfMass(Rcm);
	ZeroP();
	ZeroL();
	T_c = Temperature();
	ScaleTemperature();
	
	//printf("After adjusting P, L, and T:\n");
	//Diagnostic();
	
	//BlockCheck(0, 0, 0);
	
	Rmin2 = L_s * L_s;
	AllForcesAndPotential();
	printf("After calculating initial forces, potential:\n");
	Diagnostic();
	RminCheck();
	ForceCheck();
	
	remove("./data/temperature.dat");
	remove("./data/potential.dat");
	remove("./data/total_energy.dat");
	remove("./data/msd.dat");
			
	printf("Start simulation\n");
	
	last_rescale = 0;
	State = STARTUP;
	
	Rmin2 = L_s * L_s;
	Transfers = 0;
	
	subsample_index = 0;
	sample_index = 0;
	batch_index = 0;
	
	Start = clock();
	for (Iteration=1; Iteration <= ITERATIONS; Iteration++) {
		VelocityVerletFullStep();
		
		if ((Iteration % SAMPLE_RATE) == 0) {
			T_c = Temperature();
			
			T_array[batch_index] = T_c;
			U_array[batch_index] = U[subsample_index] / AtomCount;
			E_array[batch_index] = U[subsample_index] / AtomCount + 1.5 * T_c;
			MSD_array[batch_index] = MSD();
			
			T_sample[subsample_index] = T_c;
			E_sample[subsample_index] = E_array[batch_index];
						
			//Rmin2 = L_s * L_s;
			Transfers = 0;
			
			batch_index = (batch_index + 1) % BATCH_RATE;
			subsample_index = (subsample_index + 1) % SAMPLE_DEPTH;
		}
		
		if ((Iteration % (BATCH_RATE * SAMPLE_RATE)) == 0) {
			BatchFileWrite();
		}
		
		if ((Iteration % 1000) == 0) {
			printf("Iteration %ld\n", Iteration);
		}
		
		switch(State) {
		case STARTUP:
			if ((Iteration % 100) == 0) {
				T_c = Temperature();
				if (fabs((T_c - T_s) / T_s) > 0.10) {
					ScaleTemperature();
					//printf("Temperature scaled at iteration %d\n", Iteration);
					last_rescale = Iteration;
				}
			}
			if ((Iteration - last_rescale) > 5000) {
				sample_index = 0;
				printf("Equilibrium reached at %ld\n", Iteration);
				printf("Begin calculating thermodynamic quantities\n");
				State = EQUILIBRIUM;
			}
			break;
		case EQUILIBRIUM:			
			if (((Iteration % (SAMPLE_DEPTH * SAMPLE_RATE)) == 0) && (sample_index < SAMPLES)) {
				ThermoQuantities();
				sample_index++;
			}
			if (sample_index >= SAMPLES) {
				ActualIterations = Iteration;
				goto LoopExit;
			}
			if ((Iteration % 1000) == 0) {
				if (fabs((T_c - T_s) / T_s) > 20.0) {
					State = DIVERGE;
					ActualIterations = Iteration;
					goto LoopExit;
				}
			}
			break;
		case DIVERGE:
			ActualIterations = Iteration;
			goto LoopExit;
			break;
		}
				
	}
	
	ActualIterations = ITERATIONS;
	
	LoopExit:
	
	End = clock();
	printf("End simulation\n");
	
	
	BatchFileWrite();
	PositionsToFile();
	SpeedsToFile();
		
	Time = 1000.0 * (End - Start) / CLOCKS_PER_SEC;
	Time /= (double) ActualIterations;
	
	printf("Time per step is %f ms\n\n", Time);
	
	PrintThermoQuantities();
	
	//SpeedCheck();
	
	/* Free allocated Memory */
	FreeAll();
	
	return 0;
}

/* writes all collected data to files in batches */
void BatchFileWrite() {
	FILE *temp;
	FILE *pot;
	FILE *energy;
	FILE *msd;
	
	long Index;
	int max;
	
	temp = fopen("./data/temperature.dat", "a");
	pot = fopen("./data/potential.dat", "a");
	energy = fopen("./data/total_energy.dat", "a");
	msd = fopen("./data/msd.dat", "a");
	
	
	/* write data to files */
	if (batch_index > 0) {
		max = batch_index;
	} else {
		max = BATCH_RATE;
	}
	
	for (int i=0; i < max; i++) {
		Index = Iteration + (i - max) * SAMPLE_RATE;
		
		fprintf(temp, "%ld\t%e\n", Index, T_array[i]);
		fprintf(pot, "%ld\t%e\n", Index, U_array[i] / AtomCount);
		fprintf(energy, "%ld\t%e\n", Index, E_array[i]/ AtomCount);
		fprintf(msd, "%ld\t%e\n", Index, MSD_array[i]);
	}
	
	fclose(temp);
	fclose(pot);
	fclose(energy);
	fclose(msd);
	
	return;
}

/* calculates mean square displacement */
double MSD() {
	double MSD;
	
	MSD = 0;
	for (int i=0; i < AtomCount; i++) {
		MSD += (Rtrue[i][x] - Rinit[i][x]) * (Rtrue[i][x] - Rinit[i][x]);
		MSD += (Rtrue[i][y] - Rinit[i][y]) * (Rtrue[i][y] - Rinit[i][y]);
		MSD += (Rtrue[i][z] - Rinit[i][z]) * (Rtrue[i][z] - Rinit[i][z]);
	}
	
	MSD /= AtomCount;
	
	return sqrt(MSD);
}

/* calculates maximum velocity */
double MaxVelocity() {
	double vmax, v2;
	
	vmax = 0;
	for (int i=0; i < AtomCount; i++) {
		v2 = V[i][x] * V[i][x] + V[i][y] * V[i][y] + V[i][z] * V[i][z];
		vmax = (v2 > vmax) ? v2 : vmax;
	}
	
	return sqrt(vmax);
}

/* prints out pressure, specific heat */
void PrintThermoQuantities() {
	FILE *ef;
	FILE *cvf;
	FILE *pf;
	
	double Tsi;
	
	ef = fopen("./data/energy_vs_temperature.dat", "a");
	cvf = fopen("./data/cv_vs_temperature.dat", "a");
	pf = fopen("./data/pressure_vs_temperature.dat", "a");
		
	for (int i=0; i < sample_index; i++) {
		//Tsi = T[i] * epsilon / kB;
		fprintf(ef, "%e\t%e\n", T[i], E[i]);
		fprintf(cvf, "%e\t%e\n", T[i], Cv[i]);
		fprintf(pf, "%e\t%e\n", T[i], P[i]);
	}
	
	fclose(ef);
	fclose(cvf);
	fclose(pf);
	
	printf("Thermodynamic quantities recorded to files\n");
	
	return;
}

/* Do one iteration of integration algorithm on a single block */
void VelocityVerletFirstHalf(int bi, int bj, int bk) {
	int prev;
	//int dbx, dby, dbz;
	int ni, nj, nk;
	double dx, dy, dz;
	
	int next;
		
	prev = NULL_PREV;
	for (int i = Blocks[bi][bj][bk]; i != LIST_EMPTY; i = next) {
		//HalfKicks++;
		
		dx = V[i][x] * dt + 0.5 * F[i][x] * dt * dt;
		dy = V[i][y] * dt + 0.5 * F[i][y] * dt * dt;
		dz = V[i][z] * dt + 0.5 * F[i][z] * dt * dt;
				
		R[i][x] += dx;
		R[i][y] += dy;
		R[i][z] += dz;
		
		//PBC are not enforced on Rtrue; we want to see how far the atoms go
		Rtrue[i][x] += dx;
		Rtrue[i][y] += dy;
		Rtrue[i][z] += dz;
				
		//printf("R[%d] = (%f, %f, %f)\n", i, R[i][x], R[i][y], R[i][z]);
		
		if (R[i][x] > L_s) {
			R[i][x] -= L_s;
		} else if (R[i][x] < 0) {
			R[i][x] += L_s;
		}
		
		if (R[i][y] > L_s) {
			R[i][y] -= L_s;
		} else if (R[i][y] < 0) {
			R[i][y] += L_s;
		}
		
		if (R[i][z] > L_s) {
			R[i][z] -= L_s;
		} else if (R[i][z] < 0) {
			R[i][z] += L_s;
		}
		
		/*
		R[i][x] = fmod(R[i][x], L_s);
		R[i][y] = fmod(R[i][y], L_s);
		R[i][z] = fmod(R[i][z], L_s);
		*/
		/*
		R[i][x] += (R[i][x] < 0) * L_s;
		R[i][y] += (R[i][y] < 0) * L_s;
		R[i][z] += (R[i][z] < 0) * L_s;
		*/	
				
		ni = floor(R[i][x] / BlockSize);
		nj = floor(R[i][y] / BlockSize);
		nk = floor(R[i][z] / BlockSize);
				
		V[i][x] += 0.5 * F[i][x] * dt;
		V[i][y] += 0.5 * F[i][y] * dt;
		V[i][z] += 0.5 * F[i][z] * dt;
				
		next = LinkedList[i];
		
		//printf("prev = %d,  cur = %d,  next = %d\n", prev, i, next);
		
		if ( (ni != bi) || (nj != bj) || (nk != bk) ) {
			/*printf("\nInit:\n");
			BlockCheck(bi, bj, bk);
			HandoffCheck(ni, nj, nk);*/
			
			Transfers++;
			
			if (prev == NULL_PREV) {
				Tail(bi, bj, bk);
			} else {
				Remove(prev);
			}
			
			/*printf("\nRemove %d from Block:\n", i);
			BlockCheck(bi, bj, bk);
			HandoffCheck(ni, nj, nk);*/
			
			AppendHandoff(i, ni, nj, nk);
			
			/*printf("\nAppend %d to Handoff:\n", i);
			BlockCheck(bi, bj, bk);
			HandoffCheck(ni, nj, nk);*/
			
			//prev = prev;
			
		} else {
			prev = i;
		}
	}
		
	return;
}

/* second half-step */
void VelocityVerletSecondHalf(int bi, int bj, int bk) {
	for (int i = Blocks[bi][bj][bk]; i != LIST_EMPTY; i = LinkedList[i]) {
		//FullKicks++;
		
		V[i][x] += 0.5 * F[i][x] * dt;
		V[i][y] += 0.5 * F[i][y] * dt;
		V[i][z] += 0.5 * F[i][z] * dt;
	}
	
	return;
}

/* full Velocity-Verlet step on all blocks */
void VelocityVerletFullStep() {
	for (int i=0; i < BlocksPerSide; i++) {
		for (int j=0; j < BlocksPerSide; j++) {
			for (int k=0; k < BlocksPerSide; k++) {
				VelocityVerletFirstHalf(i, j, k);
			}
		}
	}
	
	for (int i=0; i < BlocksPerSide; i++) {
		for (int j=0; j < BlocksPerSide; j++) {
			for (int k=0; k < BlocksPerSide; k++) {
				Blocks[i][j][k] = Handoff[first][i][j][k];
			}
		}
	}
	
	AllForcesAndPotential();
	
	for (int i=0; i < BlocksPerSide; i++) {
		for (int j=0; j < BlocksPerSide; j++) {
			for (int k=0; k < BlocksPerSide; k++) {
				VelocityVerletSecondHalf(i, j, k);
			}
		}
	}
}

/* calculate the force as a function of r^2 */
double Force(double r2) {
	/* Instead of taking the square root of r^2, we use 
	 * a second-order Taylor approximation to sqrt() 
	 * centered at r^2 = 5.5.
	 * 
	 * As another shortcut, we note that the force must
	 * be multiplied by r_ij/|r_ij|, so we can make
	 * an extra division by r, thus giving us
	 * even powers of r in the denominator.
	 */
	
	double r6, r8, f;
	
	r6 = r2 * r2 * r2;
	r8 = r6 * r2;
	
	f = 24.0 * (2.0 / r6 - 1) / r8;
	
	//return f;
	
	// Subtract f(rc)/r
	f -= 1.0 / (-22.5504 + r2 * (-8.20014  + 0.248489 * r2));
	
	return f;
}

/* calculate the potential as a function of r^2 */
double Potential(double r2) {
	/* Instead of taking the square root of r^2, we use 
	 * a second-order Taylor approximation to sqrt() 
	 * centered at r^2 = 5.5.
	 */
	
	double r6, u;
	
	r6 = r2 * r2 * r2;
	
	u = 4.0 * (1.0/r6 - 1) / r6;
	
	// - u(rc)
	//u += 0.0163169;
	//return u;
	
	// - u(rc) + (rc - r) * f(rc)
	u += -0.0468836 + r2 * (0.0124721 - r2 * 0.000377942);
	
	return u;
}

/* calculate the forces (accelerations) and potential within a block 
 * Assumptions: All Forces have already been initialized to 0
 */
void InternalForcesAndPotential(int bi, int bj, int bk) {	
	double rij[3];
	double rij2;
	double fj;

	for (int i = Blocks[bi][bj][bk]; i != LIST_EMPTY; i = LinkedList[i]) {
		for (int j = LinkedList[i]; j != LIST_EMPTY; j = LinkedList[j]) {			
			rij[x] = R[j][x] - R[i][x];
			rij[y] = R[j][y] - R[i][y];
			rij[z] = R[j][z] - R[i][z];
			
			rij2 = rij[x] * rij[x] + rij[y] * rij[y] + rij[z] * rij[z];
			Rmin2 = (rij2 < Rmin2) ? rij2 : Rmin2;
						
			if (rij2 < (Rc * Rc)) {
				if (rij2 > Rsmall) {
					fj = Force(rij2);
					
					F[j][x] += fj * rij[x];
					F[j][y] += fj * rij[y];
					F[j][z] += fj * rij[z];
				
					F[i][x] -= fj * rij[x];
					F[i][y] -= fj * rij[y];
					F[i][z] -= fj * rij[z];
					
					if ((Iteration % SAMPLE_RATE) == 0) {
						Virial[subsample_index] += rij2 * fj;
						U[subsample_index] += Potential(rij2);
					}
					
					//U += Potential(rij2);
				} else {
					if ((Iteration % SAMPLE_RATE) == 0) {
						Virial[subsample_index] += rij2 * fj;
						U[subsample_index] += Umax;
					}
					//U += Umax;
				}
			}
			
		}
	}
	
	return;
}

/* calculate the forces and potential between two blocks */
void ExternalForcesAndPotential(int ai, int aj, int ak, int bi, int bj, int bk) {
	double rij[3];
	double rij2;
	double fj;
	
	// Corrections for periodic boundary conditions
	int sx, sy, sz;
	double dx, dy, dz;
	
	/* integer division returns zero unless blocks are on boundaries */
	sx = (ai - bi) / (BlocksPerSide - 1);
	sy = (aj - bj) / (BlocksPerSide - 1);
	sz = (ak - bk) / (BlocksPerSide - 1);
	
	/* these must be added to the radii of the B block to
	 * satisfy the periodic boundary conditions 
	 */
	dx = sx * L_s;
	dy = sy * L_s;
	dz = sz * L_s;
	
	//printf("sx = %d,  sy = %d,  sz = %d\n", sx, sy, sz);
	//printf("dx = %f,  dy = %f,  dz = %f\n", dx, dy, dz);
	
	for (int i = Blocks[ai][aj][ak]; i != LIST_EMPTY; i = LinkedList[i]) {
		for (int j = Blocks[bi][bj][bk]; j != LIST_EMPTY; j = LinkedList[j]) {
			rij[x] = R[j][x] - R[i][x] + dx;
			rij[y] = R[j][y] - R[i][y] + dy;
			rij[z] = R[j][z] - R[i][z] + dz;
			
			rij2 = rij[x] * rij[x] + rij[y] * rij[y] + rij[z] * rij[z];
			
			Rmin2 = (rij2 < Rmin2) ? rij2 : Rmin2;
			
			if (rij2 < (Rc * Rc)) {
				if (rij2 > Rsmall) {
					fj = Force(rij2);
					
					F[j][x] += fj * rij[x];
					F[j][y] += fj * rij[y];
					F[j][z] += fj * rij[z];
				
					F[i][x] -= fj * rij[x];
					F[i][y] -= fj * rij[y];
					F[i][z] -= fj * rij[z];
				
					
					if ((Iteration % SAMPLE_RATE) == 0) {
						Virial[subsample_index] += rij2 * fj;
						U[subsample_index] += Potential(rij2);
					}
					
					//U += Potential(rij2);
				} else {
					if ((Iteration % SAMPLE_RATE) == 0) {
						Virial[subsample_index] += rij2 * fj;
						U[subsample_index] += Umax;
					}
					//U += Umax;
				}
			}
			
		}
	}
	
	return;
}

/* calculates ALL forces (accelerations) and potential */
void AllForcesAndPotential() {
	int ip1, jp1, kp1, im1, jm1, km1;
	
	// clear Forces array
	for (int i=0; i < AtomCount; i++) {
		F[i][x] = 0;
		F[i][y] = 0;
		F[i][z] = 0;
	}
	
	//ForceCheck();
	
	if ((Iteration % SAMPLE_RATE) == 0) {
		// clear potential and virial
		U[subsample_index] = 0;
		Virial[subsample_index] = 0;
	}
	
	for (int i=0; i < BlocksPerSide; i++) {
		for (int j=0; j < BlocksPerSide; j++) {
			for (int k=0; k < BlocksPerSide; k++) {
				
				InternalForcesAndPotential(i, j, k);
								
				// nearest neighbor blocks
				// 13 of 26 cases
				ip1 = (i + 1) % BlocksPerSide;
				jp1 = (j + 1) % BlocksPerSide;
				kp1 = (k + 1) % BlocksPerSide;
				
				im1 = (i - 1 + BlocksPerSide) % BlocksPerSide;
				jm1 = (j - 1 + BlocksPerSide) % BlocksPerSide;
				km1 = (k - 1 + BlocksPerSide) % BlocksPerSide;
				
				ExternalForcesAndPotential(i, j, k,    i, j, kp1);
				ExternalForcesAndPotential(i, j, k,    i, jp1, km1);
				ExternalForcesAndPotential(i, j, k,    i, jp1, k);
				ExternalForcesAndPotential(i, j, k,    i, jp1, kp1);
				ExternalForcesAndPotential(i, j, k,    ip1, jm1, km1);
				ExternalForcesAndPotential(i, j, k,    ip1, jm1, k);
				ExternalForcesAndPotential(i, j, k,    ip1, jm1, kp1);
				ExternalForcesAndPotential(i, j, k,    ip1, j, km1);
				ExternalForcesAndPotential(i, j, k,    ip1, j, k);
				ExternalForcesAndPotential(i, j, k,    ip1, j, kp1);
				ExternalForcesAndPotential(i, j, k,    ip1, jp1, km1);
				ExternalForcesAndPotential(i, j, k,    ip1, jp1, k);
				ExternalForcesAndPotential(i, j, k,    ip1, jp1, kp1);
								
			}
		}
	}
	
	return;
}

/* calculates specific heat and pressure, records to arrays */
void ThermoQuantities() {
	double avgSqT, sqAvgT;
	double avgVir, V_s;
	double avgT;
	double avgE;
	
	avgSqT = 0;
	sqAvgT = 0;
	avgVir = 0;
	avgT = 0;
	avgE = 0;
	
	for (int i=0; i < SAMPLE_DEPTH; i++) {
		avgSqT += T_sample[i] * T_sample[i];
		//sqAvgT += T_sample[i];
		
		avgVir += Virial[i];
		
		avgT += T_sample[i];
		
		avgE += E_sample[i];
	}
		
	avgSqT /= SAMPLE_DEPTH;
	//sqAvgT /= SAMPLE_DEPTH;
	avgVir /= SAMPLE_DEPTH;
	avgT /= SAMPLE_DEPTH;
	avgE /= SAMPLE_DEPTH;
	
	sqAvgT = avgT * avgT;
	
	//T_c = Temperature();
	
	Cv[sample_index] = 1.5 / (1.0 - 1.5 * AtomCount * (avgSqT - sqAvgT) / sqAvgT);
	
	V_s = L_s * L_s * L_s;
	
	P[sample_index] = (AtomCount * avgT + 0.3333333 * avgVir) / V_s;
	
	T[sample_index] = avgT;
	
	E[sample_index] = avgE;
	
	return;
}

/* calculates sum of internal forces */
void ForceCheck() {
	Fsum[x] = 0;
	Fsum[y] = 0;
	Fsum[z] = 0;
	
	for (int i=0; i < AtomCount; i++) {
		Fsum[x] += F[i][x];
		Fsum[y] += F[i][y];
		Fsum[z] += F[i][z];
	}
	
	printf("Sum of internal forces = (%lf, %lf, %lf)\n\n", Fsum[x], Fsum[y], Fsum[z]);
	
	return;
}

void MemError() {
		fprintf(stderr, "\nError:  Not enough memory\n");
		exit(1);
	
	return;
}

/* set up initial arrays */
void Initialize() {
	printf("Initializing...\n");
	
	AtomCount = 4 * CellCount * CellCount * CellCount;
	printf("Number of particles:  %d\n", AtomCount);
	
	L_s = cbrt((double) AtomCount/rho_s);
	//L_s = CellCount * 2 * 1.12246;
	printf("Side length L_s:  %f\n", L_s);
	
	BlocksPerSide = floor(L_s/Rc);
	if (BlocksPerSide == 0) {
		BlocksPerSide = 1;
	}
	printf("Blocks per side:  %d\n", BlocksPerSide);
	
	BlockSize = L_s / (double) BlocksPerSide;
	printf("Block dimension:  %f\n", BlockSize);
	
	AllocateAll();
	
	/* place atoms in initial R */
	InitialConditions();
	
	for (int i=0; i < BlocksPerSide; i++) {
		for (int j=0; j < BlocksPerSide; j++) {
			for (int k=0; k < BlocksPerSide; k++) {
				Blocks[i][j][k] = BLOCK_EMPTY;
				
				//printf("Block (%d, %d, %d) initialized\n", i, j, k);
			}
		}
	}
	
	for (int i=0; i < AtomCount; i++) {
		LinkedList[i] = LIST_EMPTY;
	}
	
	PopulateLinkedList();
	
	Umax = Potential(Rsmall);
	
	printf("\n");
}

/* allocate all arrays */
int AllocateAll() {
	printf("Creating positions array...");
	R = (double **) malloc(AtomCount * sizeof(double *));
	if (R == NULL) {MemError();}
	for(int i=0; i < AtomCount; i++) {
		R[i] = (double *) malloc(3 * sizeof(double));
		if (R[i] == NULL) {MemError();}
	}
	printf("OK\n");
	
	printf("Creating initial positions array...");
	Rinit = (double **) malloc(AtomCount * sizeof(double *));
	if (Rinit == NULL) {MemError();}
	for(int i=0; i < AtomCount; i++) {
		Rinit[i] = (double *) malloc(3 * sizeof(double));
		if (Rinit[i] == NULL) {MemError();}
	}
	printf("OK\n");
	
	printf("Creating true positions array...");
	Rtrue = (double **) malloc(AtomCount * sizeof(double *));
	if (Rtrue == NULL) {MemError();}
	for(int i=0; i < AtomCount; i++) {
		Rtrue[i] = (double *) malloc(3 * sizeof(double));
		if (Rtrue[i] == NULL) {MemError();}
	}
	printf("OK\n");
	
	printf("Creating velocities array...");
	V = (double **) malloc(AtomCount * sizeof(double *));
	if (V == NULL) {MemError();}
	for(int i=0; i < AtomCount; i++) {
		V[i] = (double *) malloc(3 * sizeof(double));
		if (V[i] == NULL) {MemError();}
	}
	printf("OK\n");
		
	printf("Creating forces array...");
	F = (double **) malloc(AtomCount * sizeof(double *));
	if (F == NULL) {MemError();}
	for(int i=0; i < AtomCount; i++) {
		F[i] = (double *) malloc(3 * sizeof(double));
		if (F[i] == NULL) {MemError();}
	}
	printf("OK\n");
	
	printf("Creating linked list array...");
	LinkedList = (int *) malloc(AtomCount * sizeof(int));
	printf("OK\n");
		
	printf("Creating blocks array...");
	Blocks = (int ***) malloc(BlocksPerSide * sizeof(int **));
	if (Blocks == NULL) {
		MemError();
	} else {
		for(int i=0; i < BlocksPerSide; i++) {
			Blocks[i] = (int **) malloc(BlocksPerSide * sizeof(int *));
			if (Blocks[i] == NULL) {
				MemError();
			} else {
				for(int j=0; j < BlocksPerSide; j++) {
					Blocks[i][j] = (int *) malloc(BlocksPerSide * sizeof(int));
					if (Blocks[i][j] == NULL) {MemError();}
				}
			}
		}
	}
	printf("OK\n");
		
	printf("Creating handoff array...");
	for (int s = 0; s < 2; s++) {
		Handoff[s] = (int ***) malloc(BlocksPerSide * sizeof(int **));
		if (Handoff[s] == NULL) {
			MemError();
		} else {
			for(int i=0; i < BlocksPerSide; i++) {
				Handoff[s][i] = (int **) malloc(BlocksPerSide * sizeof(int *));
				if (Handoff[s][i] == NULL) {
					MemError();
				} else {
					for(int j=0; j < BlocksPerSide; j++) {
						Handoff[s][i][j] = (int *) malloc(BlocksPerSide * sizeof(int));
						if (Handoff[s][i][j] == NULL) {MemError();}
					}
				}
			}
		}
	}
	printf("OK\n");
		
	return 0;
}

/* frees all allocated memory */
void FreeAll() {
	printf("Freeing positions array...");
	for(int i=0; i < AtomCount; i++) {
		free(R[i]);
	}
	free(R);
	printf("OK\n");
	
	printf("Freeing initial positions array...");
	for(int i=0; i < AtomCount; i++) {
		free(Rinit[i]);
	}
	free(Rinit);
	printf("OK\n");
	
	printf("Freeing true positions array...");
	for(int i=0; i < AtomCount; i++) {
		free(Rtrue[i]);
	}
	free(Rtrue);
	printf("OK\n");
	
	printf("Freeing velocities array...");
	for(int i=0; i < AtomCount; i++) {
		free(V[i]);
	}
	free(V);
	printf("OK\n");
	
	printf("Freeing forces array...");
	for(int i=0; i < AtomCount; i++) {
		free(F[i]);
	}
	free(F);
	printf("OK\n");
	
	printf("Freeing linked list array...");
	free(LinkedList);
	printf("OK\n");
	
	printf("Freeing blocks array...");
	for(int i=0; i < BlocksPerSide; i++) {
		for(int j=0; j < BlocksPerSide; j++) {
			free(Blocks[i][j]);
		}
	
		free(Blocks[i]);
	}
	free(Blocks);
	printf("OK\n");
	
	printf("Freeing handoff array...");
	for (int s=0; s < 2; s++) {
		for(int i=0; i < BlocksPerSide; i++) {
			for(int j=0; j < BlocksPerSide; j++) {
				free(Handoff[s][i][j]);
			}
		
			free(Handoff[s][i]);
		}
	
		free(Handoff[s]);
	}
	printf("OK\n");
}

/* places atoms in initial R */
void InitialConditions() {
	int Index;
	double a_s;
	double r1, r2, r3, r4, r5, r6;
	
	a_s = L_s / ((double) CellCount);
	
	//a_s = 2 * 1.12246;
	//L_s = a_s * CellCount;
	
	printf("FCC Cell dimension a_s:  %f\n", a_s);
	printf("Setting up initial conditions...");
	
	for(int i=0; i < CellCount; i++) {
		for(int j=0; j < CellCount; j++) {
			for(int k=0; k < CellCount; k++) {
				Index = 4 * (i + CellCount * j + CellCount * CellCount * k);
				
				R[Index][x] = a_s * i;
				R[Index][y] = a_s * j;
				R[Index][z] = a_s * k;
				
				R[Index+1][x] = a_s * (i + 0.5);
				R[Index+1][y] = a_s * (j + 0.5);
				R[Index+1][z] = a_s * k;
				
				R[Index+2][x] = a_s * (i + 0.5);
				R[Index+2][y] = a_s * j;
				R[Index+2][z] = a_s * (k + 0.5);
				
				R[Index+3][x] = a_s * i;
				R[Index+3][y] = a_s * (j + 0.5);
				R[Index+3][z] = a_s * (k + 0.5);
				
				/* Use Box-Muller algorithm to generate velocities
				 * with a normal distribution of mean 0 and
				 * variance 1.
				 */
				
				r1 = rand()/((double)RAND_MAX + 1);
				r2 = rand()/((double)RAND_MAX + 1);
				r3 = rand()/((double)RAND_MAX + 1);
				r4 = rand()/((double)RAND_MAX + 1);
				r5 = rand()/((double)RAND_MAX + 1);
				r6 = rand()/((double)RAND_MAX + 1);
				
				V[Index][x] = sqrt(-2.0 * log(r1)) * cos(2.0 * PI * r2);
				V[Index][y] = sqrt(-2.0 * log(r3)) * cos(2.0 * PI * r4);
				V[Index][z] = sqrt(-2.0 * log(r5)) * cos(2.0 * PI * r6);
				
				r1 = rand()/((double)RAND_MAX + 1);
				r2 = rand()/((double)RAND_MAX + 1);
				r3 = rand()/((double)RAND_MAX + 1);
				r4 = rand()/((double)RAND_MAX + 1);
				r5 = rand()/((double)RAND_MAX + 1);
				r6 = rand()/((double)RAND_MAX + 1);
				
				V[Index+1][x] = sqrt(-2.0 * log(r1)) * cos(2.0 * PI * r2);
				V[Index+1][y] = sqrt(-2.0 * log(r3)) * cos(2.0 * PI * r4);
				V[Index+1][z] = sqrt(-2.0 * log(r5)) * cos(2.0 * PI * r6);
				
				r1 = rand()/((double)RAND_MAX + 1);
				r2 = rand()/((double)RAND_MAX + 1);
				r3 = rand()/((double)RAND_MAX + 1);
				r4 = rand()/((double)RAND_MAX + 1);
				r5 = rand()/((double)RAND_MAX + 1);
				r6 = rand()/((double)RAND_MAX + 1);
				
				V[Index+2][x] = sqrt(-2.0 * log(r1)) * cos(2.0 * PI * r2);
				V[Index+2][y] = sqrt(-2.0 * log(r3)) * cos(2.0 * PI * r4);
				V[Index+2][z] = sqrt(-2.0 * log(r5)) * cos(2.0 * PI * r6);
				
				r1 = rand()/((double)RAND_MAX + 1);
				r2 = rand()/((double)RAND_MAX + 1);
				r3 = rand()/((double)RAND_MAX + 1);
				r4 = rand()/((double)RAND_MAX + 1);
				r5 = rand()/((double)RAND_MAX + 1);
				r6 = rand()/((double)RAND_MAX + 1);
				
				V[Index+3][x] = sqrt(-2.0 * log(r1)) * cos(2.0 * PI * r2);
				V[Index+3][y] = sqrt(-2.0 * log(r3)) * cos(2.0 * PI * r4);
				V[Index+3][z] = sqrt(-2.0 * log(r5)) * cos(2.0 * PI * r6);
				
				/*
				V[Index+1][x] = rand()/((double)RAND_MAX + 1);
				V[Index+1][y] = rand()/((double)RAND_MAX + 1);
				V[Index+1][z] = rand()/((double)RAND_MAX + 1);
				
				V[Index+2][x] = rand()/((double)RAND_MAX + 1);
				V[Index+2][y] = rand()/((double)RAND_MAX + 1);
				V[Index+2][z] = rand()/((double)RAND_MAX + 1);
				
				V[Index+3][x] = rand()/((double)RAND_MAX + 1);
				V[Index+3][y] = rand()/((double)RAND_MAX + 1);
				V[Index+3][z] = rand()/((double)RAND_MAX + 1);
				*/
			}
		}
	}
	
	//copy to Rinit array, for calculating MSD
	for (int i=0; i < AtomCount; i++) {
		Rinit[i][x] = R[i][x];
		Rinit[i][y] = R[i][y];
		Rinit[i][z] = R[i][z];
		
		Rtrue[i][x] = R[i][x];
		Rtrue[i][y] = R[i][y];
		Rtrue[i][z] = R[i][z];
	}
	
	printf("OK\n");
}

// Puts origin at center
void Recenter() {	
	CenterOfMass(Rcm);
	
	for(int i=0; i < AtomCount; i++) {
		R[i][x] -= Rcm[x];
		R[i][y] -= Rcm[y];
		R[i][z] -= Rcm[z];
	}
		
	return;
}

/* zeroes the linear momentum */
void ZeroP() {
	double P[3];
	
	TotalP(P);
	
	for(int i=0; i < AtomCount; i++) {
		V[i][x] -= P[x];
		V[i][y] -= P[y];
		V[i][z] -= P[z];
	}
		
	return;
}

/* zeroes the angular momentum */
void ZeroL() {
	double L[3];
	double I[3][3];
	
	MomentOfInertia(I);
	
	TotalL(L);
	
	for(int i=0; i < AtomCount; i++) {
		//Note:  Inertia tensor is diagonal due to symmetry of initial conditions
		V[i][x] += ((R[i][y] - Rcm[y]) * L[z] - (R[i][z] - Rcm[z]) * L[y]) / I[x][x];
		V[i][y] += ((R[i][z] - Rcm[z]) * L[x] - (R[i][x] - Rcm[x]) * L[z]) / I[y][y];
		V[i][z] += ((R[i][x] - Rcm[x]) * L[y] - (R[i][y] - Rcm[y]) * L[x]) / I[z][z];
	}
	
	return;
}

 /* sets temperature to T_s by scaling V */
void ScaleTemperature() {
	//double T_c;
	double alpha;
	
	//T_c = Temperature();
	
	alpha = sqrt(T_s / T_c);
	
	for (int i=0; i < AtomCount; i++) {
		V[i][x] *= alpha;
		V[i][y] *= alpha;
		V[i][z] *= alpha;
	}
}

/* calculates CoM */
void CenterOfMass(double Rcm[3]) {
	Rcm[x] = 0;
	Rcm[y] = 0;
	Rcm[z] = 0;
	
	for(int i=0; i < AtomCount; i++) {
		Rcm[x] += R[i][x];
		Rcm[y] += R[i][y];
		Rcm[z] += R[i][z];
	}
	
	Rcm[x] = Rcm[x] / (double) AtomCount;
	Rcm[y] = Rcm[y] / (double) AtomCount;
	Rcm[z] = Rcm[z] / (double) AtomCount;
	
	return;
}

/* calculates total linear momentum */
void TotalP(double P[3]) {
	P[x] = 0;
	P[y] = 0;
	P[z] = 0;
	
	for (int i=0; i < AtomCount; i++) {
		P[x] += V[i][x];
		P[y] += V[i][y];
		P[z] += V[i][z];
	}
	
	P[x] = P[x] / (double) AtomCount;
	P[y] = P[y] / (double) AtomCount;
	P[z] = P[z] / (double) AtomCount;
	
	return;
}

/* calculates total angular momentum */
void TotalL(double L[3]) {
	L[x] = 0;
	L[y] = 0;
	L[z] = 0;
	
	for (int i=0; i < AtomCount; i++) {
		L[x] += (R[i][y] - Rcm[y]) * V[i][z] - (R[i][z] - Rcm[z]) * V[i][y];
		L[y] += (R[i][z] - Rcm[z]) * V[i][x] - (R[i][x] - Rcm[x]) * V[i][z];
		L[z] += (R[i][x] - Rcm[x]) * V[i][y] - (R[i][y] - Rcm[y]) * V[i][x];
	}
		
	return;
}

/* calculates moment of inertia */
void MomentOfInertia(double I[3][3]) {
	double A, B;
	
	for(int i=0; i < 3; i++) {
		for(int j=0; j < 3; j++) {
			I[i][j] = 0;
		}
	}
	
	for(int k=0; k < AtomCount; k++) {
		for(int i=0; i < 3; i++) {
			for(int j=0; j < 3; j++) {
				A = (i == j) * ((R[k][x] - Rcm[x]) * (R[k][x] - Rcm[x])
											+ (R[k][y] - Rcm[y]) * (R[k][y] - Rcm[y])
										  + (R[k][z] - Rcm[z]) * (R[k][z] - Rcm[z]));
				B = (R[k][i] - Rcm[i]) * (R[k][j] - Rcm[j]);
				I[i][j] += A - B;
			}
		}
	}
	
}

/* calculates temperature */
double Temperature() {
	double T_s;
	
	T_s = 0;
	for (int i=0; i < AtomCount; i++) {
		T_s += V[i][x] * V[i][x];
		T_s += V[i][y] * V[i][y];
		T_s += V[i][z] * V[i][z];
	}
	
	return T_s / (3.0 * (double) AtomCount);
}

/* sets up linked list array */
void PopulateLinkedList() {
	double rx, ry, rz;
	int bi, bj, bk;
	
	for (int i=0; i < AtomCount; i++) {
		rx = R[i][x] + Rcm[x];
		ry = R[i][y] + Rcm[y];
		rz = R[i][z] + Rcm[z];
		
		bi = floor(rx / BlockSize);
		bj = floor(ry / BlockSize);
		bk = floor(rz / BlockSize);
		
		AppendBlocks(i, bi, bj, bk);
	}
	
	//duplicate Blocks array to Handoff
	for (int i=0; i < BlocksPerSide; i++) {
		for (int j=0; j < BlocksPerSide; j++) {
			for (int k=0; k < BlocksPerSide; k++) {
				Handoff[first][i][j][k] = Blocks[i][j][k];
			}
		}
	}
	
	return;			
}

/* appends to linked list */
void AppendBlocks(int new, int i, int j, int k) {
	LinkedList[new] = Blocks[i][j][k];
	Blocks[i][j][k] = new;
	
	return;
}

/* appends to linked list and updates Handoff */
void AppendHandoff(int new, int i, int j, int k) {
	LinkedList[new] = Handoff[first][i][j][k];
	
	if (Handoff[first][i][j][k] == Blocks[i][j][k]) {
		Handoff[last][i][j][k] = new;
	}
	
	Handoff[first][i][j][k] = new;
	
	return;
}

/* appends Handoff sublist to Blocks linked list */
void AppendSublistBlocks(int head, int tail, int i, int j, int k) {
	LinkedList[tail] = Blocks[i][j][k];
	Blocks[i][j][k] = head;
	
	return;
}

/* removes from the middle of a list */
void Remove(int prev) {
	LinkedList[prev] = LinkedList[LinkedList[prev]];
	
	return;
}

/* removes first element from list */
void Tail(int i, int j, int k) {
	if (Handoff[first][i][j][k] == Blocks[i][j][k]) {
		Handoff[first][i][j][k] = LinkedList[Blocks[i][j][k]];
	} else {
		LinkedList[ Handoff[last][i][j][k] ] = LinkedList[ Blocks[i][j][k] ];
	}
	
	Blocks[i][j][k] = LinkedList[Blocks[i][j][k]];
	
	return;
}

/* list what atoms are in a block */
void BlockCheck(int bi, int bj, int bk) {
	int count;
	char dummy[15];
	
	printf("  Block (%d, %d, %d):  -1", bi, bj, bk);
	
	count = 0;
	for (int i = Blocks[bi][bj][bk]; i != LIST_EMPTY; i = LinkedList[i]) {
		printf(" -> %d", i);
		
		//stop execution of program on suspected infinite lists
		if ( (count > 100) ) {
			printf("\n");
			scanf("%s", dummy);
			return;
		}
		
		count++;
	}
	
	printf(" -> -1\n");
	
	return;
}

/* list what atoms are in a handoff block */
void HandoffCheck(int bi, int bj, int bk) {
	int count;
	char dummy[15];
	
	printf("  Handoff (%d, %d, %d):  -1", bi, bj, bk);
	
	count = 0;
	for (int i = Handoff[first][bi][bj][bk]; i != LIST_EMPTY; i = LinkedList[i]) {
		printf(" -> %d", i);
		
		//stop execution of program on suspected infinite lists
		if ( (count > 100) ) {
			printf("\n");
			scanf("%s", dummy);
			return;
		}
		
		count++;
	}
	
	printf(" -> -1\n");
	
	return;
}

/* checks to see if all lists are terminated */
void ListCheck() {
	int count;
	
	count = 0;
	for (int i=0; i < AtomCount; i++) {
		if (LinkedList[i] == LIST_EMPTY) {
			count++;
		}
	}
	
	printf("%d list terminations found.\n", count);
	
	return;
}

/* print current total energy */
void EnergyCheck() {
	double E_s;
	double KE;
	
	T_c = Temperature();
	
	KE = 1.5 * AtomCount * T_c;
	E_s = U[subsample_index] + KE;
	
	printf("Potential:  %lf     Kinetic:  %lf\n", U[subsample_index], KE);
	printf("Total energy:  %lf\n", E_s);
	
	return;
}

/* print out CoM, P, and L for diagnostic */
void Diagnostic() {
	double CM[3];
	double P[3];
	double L[3];
	double T_c, KE, E_tot;
	double T_err;
	
	CenterOfMass(CM);
	TotalP(P);
	TotalL(L);
	T_c = Temperature();
	T_err = 100.0 * (T_c - T_s) / T_s;
	KE = 1.5 * AtomCount * T_c;
	E_tot = U[subsample_index] + KE;
	
	printf("  CM: (%f, %f, %f)\n", CM[x], CM[y], CM[z]);
	printf("   P: (%f, %f, %f)\n", P[x], P[y], P[z]);
	printf("   L: (%f, %f, %f)\n", L[x], L[y], L[z]);
	printf("   T:  %f;   T_s:  %f  Error:  %f\n", T_c, T_s, T_err);
	printf("  PE:  %f;    KE:  %f\n", U[subsample_index], KE);
	printf("   E = KE + PE:  %f\n", E_tot);
	
	printf("\n");
	
	return;
}

/* prints out closest approach */
void RminCheck() {
	printf("Closest approach:  %lf\n", sqrt(Rmin2));
	
	return;
}

/* prints out all positions */
void PositionCheck() {
	printf("Positions:\n");
	for (int i=0; i < AtomCount; i++) {
		printf("  R[%d] = (%f, %f, %f)\n", i, R[i][x], R[i][y], R[i][z]);
	}
	
	return;
}

/* prints out all positions to file */
void PositionsToFile() {
	FILE *pfile;
	FILE *dim;
	
	dim = fopen("./data/box_dimension.dat", "w");
	fprintf(dim, "%e\n", L_s);
	fclose(dim);
	
	pfile = fopen("./data/final_positions.dat", "w");
	
	for (int i=0; i < AtomCount; i++) {
		fprintf(pfile, "%e\t%e\t%e\n", R[i][x], R[i][y], R[i][z]);
	}
	
	fclose(pfile);
	
	return;
}

/* prints out all speeds to file */
void SpeedsToFile() {
	double speed;
	FILE *sfile;
	
	sfile = fopen("./data/final_speeds.dat", "w");
	
	for (int i=0; i < AtomCount; i++) {
		speed = sqrt(V[i][x] * V[i][x] + V[i][y] * V[i][y] + V[i][z] * V[i][z]);
		//speed *= sigma / tau;
		fprintf(sfile, "%f\n", speed);
	}
	
	fclose(sfile);
	
	return;
}

/* prints out all velocities */
void VelocityCheck() {
	printf("Velocities:\n");
	for (int i=0; i < AtomCount; i++) {
		printf("  V[%d] = (%f, %f, %f)\n", i, V[i][x], V[i][y], V[i][z]);
	}
	
	return;
}

/* prints out all speeds */
void SpeedCheck() {
	double speed;
	
	printf("Speeds:\n");
	for (int i=0; i < AtomCount; i++) {
		speed = V[i][x] * V[i][x] + V[i][y] * V[i][y] + V[i][z] * V[i][z];
		speed = sqrt(speed);
		
		printf("%d:  %lf\n", i, speed);
	}
	
	return;
}


