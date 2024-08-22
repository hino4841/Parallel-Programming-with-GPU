
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <device_functions.h>
#include <stdio.h>
#include <cstdlib>
#include <math.h>
#include<time.h>


#include <cuda.h>
#include <cuda_runtime_api.h>




#define Number 1000
#define Delta_t 0.01

__global__
void Simulate(double* Vortex_p, double* Omega_v_p, double* VortexN_p, double* Omega_vN_p, double *Sigma_p)
{
	double radiika_p, t1_p, t2_p;
		double t3_p, Om22P_p, ssss_p, vxx_p, vyy_p, vzz_p;
		double dvxdxmov, dvxdymov, dvxdzmov;
		double  dvydxmov, dvydymov, dvydzmov;
		double dvzdxmov,  dvzdymov,  dvzdzmov;

	double Vxc, Vyc, Vzc, dssss_dr;
	
	int j =  3*(threadIdx.x + blockIdx.x*blockDim.x);

	if (true) {
		Vxc = 1.0;
		Vyc = 0;
		Vzc = 0;
		dvxdxmov = 0.0;
		dvxdymov = 0;
		dvxdzmov = 0;

		dvydxmov = 0;
		dvydymov = 0;
		dvydzmov = 0;

		dvzdxmov = 0;
		dvzdymov = 0;
		dvzdzmov = 0;

	}
	
	for (int i=0; i < Number; i++) {
		vxx_p = Vortex_p[j] - Vortex_p[i * 3];
		vyy_p = Vortex_p[j + 1] - Vortex_p[(i * 3) + 1];
		vzz_p = Vortex_p[j + 2] - Vortex_p[(i * 3) + 2];
		radiika_p = vxx_p*vxx_p + vyy_p*vyy_p + vzz_p*vzz_p;
		t1_p = vyy_p*Omega_v_p[(i * 3) + 2] - vzz_p*Omega_v_p[(i * 3) + 1];
		t2_p = vzz_p*Omega_v_p[i * 3] - vxx_p*Omega_v_p[(i * 3) + 2];
		t3_p = vxx_p*Omega_v_p[(i * 3) + 1] - vyy_p*Omega_v_p[i * 3];
		Om22P_p = 3.1416 / Sigma_p[i] / Sigma_p[i] / 2.0;
		ssss_p = exp(-radiika_p*Om22P_p);

		Vxc = Vxc + ssss_p*t1_p;
		Vyc = Vyc + ssss_p*t2_p;
		Vzc = Vzc + ssss_p*t3_p;

		dssss_dr = (-Om22P_p)*ssss_p;

		dvxdxmov = dssss_dr*vxx_p*t1_p + dvxdxmov;
		dvxdymov = dssss_dr*vyy_p*t1_p + Omega_v_p[(i * 3) + 2] * ssss_p + dvxdymov;
		dvxdzmov = dssss_dr*vzz_p*t1_p - Omega_v_p[(i * 3) + 1] * ssss_p + dvxdzmov;

		dvydxmov = dssss_dr*vxx_p*t2_p - Omega_v_p[(i * 3) + 2] * ssss_p + dvydxmov;
		dvydymov = dssss_dr*vyy_p*t2_p + dvydymov;
		dvydzmov = dssss_dr*vzz_p*t2_p + Omega_v_p[i * 3] * ssss_p + dvydzmov;

		dvzdxmov = dssss_dr*vxx_p*t3_p + Omega_v_p[(i * 3) + 1] * ssss_p + dvzdxmov;
		dvzdymov = dssss_dr*vyy_p*t3_p - Omega_v_p[i * 3] * ssss_p + dvzdymov;
		dvzdzmov = dssss_dr*vzz_p*t3_p + dvzdzmov;
	}

	if ( true) {
		
		VortexN_p[j] = Vortex_p[j] + Delta_t*Vxc;
		VortexN_p[j + 1] = Vortex_p[j + 1] + Delta_t*Vyc;
		VortexN_p[j + 2] = Vortex_p[j + 2] + Delta_t*Vzc;

		

		double domxdt, domydt, domzdt;
		domxdt = dvxdxmov*Omega_v_p[j] + dvxdymov*Omega_v_p[j + 1] + dvxdzmov*Omega_v_p[j + 2];
		domydt = dvydxmov*Omega_v_p[j] + dvydymov*Omega_v_p[j + 1] + dvydzmov*Omega_v_p[j + 2];
		domzdt = dvzdxmov*Omega_v_p[j] + dvzdymov*Omega_v_p[j + 1] + dvzdzmov*Omega_v_p[j + 2];


		Omega_vN_p[j] = Omega_v_p[j] + domxdt*Delta_t;
		Omega_vN_p[j + 1] = Omega_v_p[j + 1] + domydt*Delta_t;
		Omega_vN_p[j + 2] = Omega_v_p[j + 2] + domzdt*Delta_t;
		Vxc = 0, Vyc = 0, Vzc = 0;
		dvxdxmov = 0, dvxdymov = 0, dvxdzmov = 0;
		dvydxmov = 0, dvydymov = 0, dvydzmov = 0;
		dvzdxmov = 0, dvzdymov = 0, dvzdzmov = 0;

	}
		
	
	
}





int main()
{

	const int Ntime = 10000;
	//const double Delta_t = 0.01;
	const double Radius = 0.1;
	//const int Number = 10;
	//const double V_mean = 1.0;
	double Vortex[Number][3];
	double Omega_v[Number][3];
	double VortexN[Number][3];
	double Omega_vN[Number][3];
	
		double *Vortex_p=new double[Number*3];
		double *Omega_v_p = new double[Number * 3];
		double *VortexN_p = new double[Number * 3];
		double *Omega_vN_p = new double[Number * 3];

	double Sigma[Number];
	double *Sigma_p;
	double StatisticalMoments[4] = {0.000};
	double Amagni=0.0,Amagnit_old,Amagnit_new,Speed_max,Sigmas;
	double Energy=0;
	int Ncout=0;
	FILE *fp1,*fp2;
	fp1 = fopen("D:\\cudaa\\Velocities1.txt","w+");
	fp2 = fopen("D:\\cudaa\\MaxValue1.txt","w+");
	double Vx;
	double vxx, vyy, vzz;
	double *vxx_p, *vyy_p, *vzz_p;
	clock_t time0,time1;
	for (int ivorton = 0; ivorton < Number; ivorton++) {
		
		Vortex[ivorton][0] = (double)rand() / (double)RAND_MAX;
		Vortex[ivorton][1] = (double)rand() / (double)RAND_MAX;
		Vortex[ivorton][2] = (double)rand() / (double)RAND_MAX;

		Omega_v[ivorton][0] = (((double)rand() / (double)RAND_MAX) - 0.5);
		Omega_v[ivorton][1] = (((double)rand() / (double)RAND_MAX) - 0.5);
		Omega_v[ivorton][2] = (((double)rand() / (double)RAND_MAX) - 0.5);
		Sigma[ivorton] = Radius;
		

		//printf("%f",Vortex[ivorton][1]);

	}

	int counter = 0;
	for (int h = 0; h < Number; h++) {
		for (int w = 0; w < 3; w++) {
			Vortex_p[counter] = Vortex[h][w];
			Omega_v_p[counter] = Omega_v[h][w];
			VortexN_p[counter] = VortexN[h][w];
			Omega_vN_p[counter] = Omega_vN[h][w];
			counter++;
		}
	}


	




	time0 = clock();
	double radiika;
	double t1, t2, t3;
	double Om22P;
	double ssss,dssss_dr;
	double domxdt, domydt, domzdt,Replace=0;
	double Vxc, Vyc , Vzc ;
	double dvxdxmov , dvxdymov , dvxdzmov ;
	double  dvydxmov , dvydymov , dvydzmov ;
	double dvzdxmov , dvzdymov , dvzdzmov ;

	




	double *radiika_p;
	double *t1_p, *t2_p, *t3_p;
	double *Om22P_p;
	double *ssss_p, *dssss_dr_p;
	double *domxdt_p, *domydt_p, *domzdt_p, *Replace_p ;
	double *Vxc_p, *Vyc_p, *Vzc_p;
	double *dvxdxmov_p, *dvxdymov_p, *dvxdzmov_p;
	double  *dvydxmov_p, *dvydymov_p, *dvydzmov_p;
	double *dvzdxmov_p, *dvzdymov_p, *dvzdzmov_p;
	vxx_p = &vxx;
	cudaMalloc((void**)&Vortex_p, (Number * 3) * sizeof(double));
	cudaMalloc((void**)&Omega_v_p, (Number * 3) * sizeof(double));
	cudaMalloc((void**)&VortexN_p, (Number * 3) * sizeof(double));
	cudaMalloc((void**)&Omega_vN_p, (Number * 3) * sizeof(double));
	cudaMalloc((void**)&Sigma_p, (Number) * sizeof(double));
	

	

	for (int itime = 0; itime < Ntime; itime++) {
		printf("%*d %f %e %d \n ",4,itime,Amagni,Energy,Ncout);
		//cudaMalloc((void**)&domzdt_p, sizeof(double));

		cudaMemcpy(Vortex_p, Vortex, (Number * 3) * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(Omega_v_p, Omega_v, (Number * 3) * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(VortexN_p, VortexN, (Number * 3) * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(Omega_vN_p, Omega_vN, (Number * 3) * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(Sigma_p, Sigma, (Number) * sizeof(double), cudaMemcpyHostToDevice);
		

		dim3 dimBlock(Number, 1);
		dim3 dimGrid(Number, 1);


		Simulate <<<Number/32, 32 >>> (Vortex_p,Omega_v_p, VortexN_p, Omega_vN_p,Sigma_p);

		cudaMemcpy(Vortex, Vortex_p,  (Number * 3) * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(Omega_v,Omega_v_p, (Number * 3) * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(VortexN, VortexN_p,  (Number * 3) * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(Omega_vN, Omega_vN_p,  (Number * 3) * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(Sigma, Sigma_p,  (Number) * sizeof(double), cudaMemcpyDeviceToHost);
		
	




		/*
		for (int ivorton = 0; ivorton < Number;ivorton++) {
			double Vxc=V_mean, Vyc=0.0, Vzc=0.0;
			double dvxdxmov = 0.0, dvxdymov = 0.0, dvxdzmov = 0.0;
			double  dvydxmov = 0.0, dvydymov = 0.0, dvydzmov = 0.0;
			double dvzdxmov = 0.0, dvzdymov = 0.0, dvzdzmov = 0.0;

			for (int induced = 0; induced < Number; induced++) {
				 vxx = Vortex[ivorton][0] - Vortex[induced][0];
				 vyy = Vortex[ivorton][1] - Vortex[induced][1];
				 vzz = Vortex[ivorton][2] - Vortex[induced][2];
				radiika = vxx*vxx + vyy*vyy + vzz*vzz;
				t1 = vyy*Omega_v[induced][2] - vzz*Omega_v[induced][1];
				t2 = vzz*Omega_v[induced][0] - vzz*Omega_v[induced][2];
				t3 = vxx*Omega_v[induced][1] - vzz*Omega_v[induced][0];
				Om22P = 3.1416 / Sigma[induced] / Sigma[induced] / 2.0;
				ssss = exp(-radiika* Om22P);
				 
				Vxc = Vxc + ssss*t1;
				Vyc = Vyc + ssss*t2;
				Vzc = Vzc + ssss*t3;

				dssss_dr = (-Om22P)*ssss;

				dvxdxmov = dssss_dr*vxx*t1 + dvxdxmov;
				dvxdymov = dssss_dr*vyy*t1 + Omega_v[induced][2]*ssss + dvxdymov;
				dvxdzmov = dssss_dr*vzz*t1 - Omega_v[induced][1]*ssss + dvxdzmov;

				dvydxmov = dssss_dr*vxx*t2 - Omega_v[induced][2]*ssss + dvydxmov;
				dvydymov = dssss_dr*vyy*t2 + dvydymov;
				dvydzmov = dssss_dr*vzz*t2 + Omega_v[induced][0]*ssss + dvydzmov;

				dvzdxmov = dssss_dr*vxx*t3 + Omega_v[induced][1]*ssss + dvzdxmov;
				dvzdymov = dssss_dr*vyy*t3 - Omega_v[induced][0]*ssss + dvzdymov;
				dvzdzmov = dssss_dr*vzz*t3 + dvzdzmov;




			}

			VortexN[ivorton][0] = Vortex[ivorton][0] + Delta_t*Vxc;
			VortexN[ivorton][1] = Vortex[ivorton][1] + Delta_t*Vyc;
			VortexN[ivorton][2] = Vortex[ivorton][2] + Delta_t*Vzc;

		//	domxdt=dvxdxmov*Omega_v[ivorton][0]+dvxdymov*Omega_v[][]

			domxdt = dvxdxmov*Omega_v[ivorton][0] + dvxdymov*Omega_v[ivorton][1] + dvxdzmov*Omega_v[ivorton][2];
			domydt = dvydxmov*Omega_v[ivorton][0] + dvydymov*Omega_v[ivorton][1] + dvydzmov*Omega_v[ivorton][2];
			domzdt = dvzdxmov*Omega_v[ivorton][0] + dvzdymov*Omega_v[ivorton][1] + dvzdzmov*Omega_v[ivorton][2];
			Omega_vN[ivorton][0] = Omega_v[ivorton][0] + domxdt*Delta_t;
			Omega_vN[ivorton][1]= Omega_v[ivorton][1] + domydt*Delta_t;
			Omega_vN[ivorton][2] = Omega_v[ivorton][2] + domzdt*Delta_t;



		}
		*/
		Ncout = 0;
		for (int ivorton = 0; ivorton < Number; ivorton++) {
			Replace = 0.0;
			for (int kkk = 0; kkk < 3; kkk++) {
				if (VortexN[ivorton][kkk] < 0.0) {
					Replace = 1.0;
				}
				if (VortexN[ivorton][kkk] > 1.0) {
					Replace = 1.0;
				}

			}
			if (Replace == 1.0) {
				Ncout = Ncout + 1;

				VortexN[ivorton][0] = (double)rand() / (double)RAND_MAX;
				VortexN[ivorton][1] = (double)rand() / (double)RAND_MAX;
				VortexN[ivorton][2] = (double)rand() / (double)RAND_MAX;

				Omega_vN[ivorton][0] = (((double)rand() / (double)RAND_MAX) - 0.5);
				Omega_vN[ivorton][1] = (((double)rand() / (double)RAND_MAX) - 0.5);
				Omega_vN[ivorton][2] = (((double)rand() / (double)RAND_MAX) - 0.5);
				Sigma[ivorton] = Radius;

			}
		}
		Amagni = 0.0;
		for (int ivorton = 0; ivorton < Number; ivorton++) {
			Vortex[ivorton][0] = VortexN[ivorton][0];
			Vortex[ivorton][1] = VortexN[ivorton][1];
			Vortex[ivorton][2] = VortexN[ivorton][2];

			Amagnit_old = sqrt((Omega_v[ivorton][0] * Omega_v[ivorton][0]) +( Omega_v[ivorton][1] * Omega_v[ivorton][1]) +( Omega_v[ivorton][2] * Omega_v[ivorton][2]));

			Omega_v[ivorton][0] = Omega_vN[ivorton][0];
			Omega_v[ivorton][1] = Omega_vN[ivorton][1];
			Omega_v[ivorton][2] = Omega_vN[ivorton][2];

			Amagnit_new= sqrt((Omega_v[ivorton][0] * Omega_v[ivorton][0]) + (Omega_v[ivorton][1] * Omega_v[ivorton][1]) + (Omega_v[ivorton][2] * Omega_v[ivorton][2]));
			Sigma[ivorton] = Sigma[ivorton] * sqrt(Amagnit_old / Amagnit_new);
			if (Amagnit_new >= Amagni) {
				Amagni = Amagnit_new;
				Energy = (Amagnit_new*Amagnit_new)*(pow(Sigma[ivorton], 5));
				Speed_max = Amagnit_new*Sigma[ivorton];
				Sigmas = Sigma[ivorton];
			}
		}

		//file write
		fprintf(fp1, "%f %f %f %f %f \n", itime*Delta_t, Amagni, Energy, Speed_max, Sigmas);
		Vx = 0.0;
		for (int induced = 0; induced < Number; induced++) {
			vxx = 0.5 - Vortex[induced][0];
			vyy = 0.5 - Vortex[induced][1];
			vzz = 0.5 - Vortex[induced][2];
			radiika = vxx*vxx + vyy*vyy + vzz*vzz;
			t1 = vyy*	Omega_v[induced][2] - vzz*Omega_v[induced][1];
			Om22P = 3.1416 / Sigma[induced] / Sigma[induced] / 2.0;
			ssss = exp(-radiika*Om22P);
			Vx = Vx + ssss*t1;

		}
		fprintf(fp2, "%f %f \n", itime*Delta_t, Vx);
		for (int ier = 0; ier < 4; ier++) {
			StatisticalMoments[ier] = StatisticalMoments[ier]  +pow(Vx, ier);
		}
		
	}
		cudaFree(Vortex_p);
		cudaFree(Omega_v_p);
		cudaFree(VortexN_p);
		cudaFree(Omega_vN_p);
		cudaFree(Sigma_p);
		

	fclose(fp1);
	fclose(fp2);
	time1 = clock();
	printf("Time taken for execution with GPU acceleration (~31blocks & 32 threads)= %f sec(s)", (double)(time1 - time0) / CLOCKS_PER_SEC);

    
}

