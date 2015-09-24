#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex_math.h>

#define TOINDEX(a,b,c) ((c)*sizeXY*sizeXY + (b)*sizeXY + (a)%sizeXY)

typedef struct microscopeParams
{

	double NA;
	double NoilReal;
	double NoilDesign;

	double NcovReal;
	double NcovDesign;

	double NspecReal;
	double NspecDesign;

	double ToilReal;
	double ToilDesign;

	double TcovReal;
	double TcovDesign;

	double TspecReal;
	double TspecDesign;
	
	double zdReal;
	double zdDesign;

	double M;

}microscopeParams;

gsl_complex OPD(microscopeParams params, double zdEffective, double rho, double deltaZ)
{
	double *diff, diffImag, diffReal;

	double secOrderOilReal =    1-gsl_pow_2(params.NA*rho/params.NoilReal);
	
	double secOrderSpecReal =   1-gsl_pow_2(params.NA*rho/params.NspecReal);
	double secOrderCovReal =    1-gsl_pow_2(params.NA*rho/params.NcovReal);
	
	double secOrderOilDesign =  1-gsl_pow_2(params.NA*rho/params.NoilDesign);
	double secOrderSpecDesign = 1-gsl_pow_2(params.NA*rho/params.NspecDesign);
	double secOrderCovDesign =  1-gsl_pow_2(params.NA*rho/params.NcovDesign);

	diffImag=0, diffReal=0;
	
	if(secOrderOilReal < 0)
	{
		diff = &diffImag;
	}else{
		diff = &diffReal;
	}
	
	double T, N;
	
	T = (params.zdDesign-params.zdReal)*pow(zdEffective,2)*params.NoilReal;
	N = (params.zdDesign*params.zdReal*gsl_pow_2(params.NA));
	
	*diff =  params.NoilReal*(deltaZ+T/N);
	*diff -= params.NspecReal*params.TspecReal * gsl_pow_2(params.NoilReal/params.NspecReal);
	*diff -= params.NcovReal*params.TcovReal * gsl_pow_2(params.NoilReal/params.NcovReal);
	*diff += params.NcovDesign*params.TcovDesign * gsl_pow_2(params.NoilReal/params.NcovDesign);
	*diff += params.NoilDesign*params.ToilDesign * gsl_pow_2(params.NoilReal/params.NoilDesign);

	*diff *= sqrt(fabs(secOrderOilReal));
	
	diffReal += gsl_pow_2(zdEffective*rho)*(params.zdDesign-params.zdReal)/(2*params.NoilReal*params.zdDesign*params.zdReal);

	if(secOrderSpecReal < 0)
	{
		diff = &diffImag;
	}else{
		diff = &diffReal;
	}
	
	*diff += params.NspecReal*params.TspecReal * sqrt(fabs(secOrderSpecReal));

	if(secOrderCovReal < 0)
	{
		diff = &diffImag; 
	}else{
		diff = &diffReal;
	}	
	
	*diff += params.NcovReal*params.TcovReal * sqrt(fabs(secOrderCovReal));

	if(secOrderCovDesign < 0)	
	{
		diff = &diffImag;
	}else{
		diff = &diffReal;
	}	
	
	*diff -= params.NcovDesign*params.TcovDesign * sqrt(fabs(secOrderCovDesign));

	if(secOrderOilDesign < 0)
	{
		diff = &diffImag;
	}else{
		diff = &diffReal;
	}	
	
	*diff -= params.NoilDesign*params.ToilDesign * sqrt(fabs(secOrderOilDesign));

	return gsl_complex_rect(diffReal,diffImag);

}

double Intensity(microscopeParams params, double lambda, double R, double deltaZ) 
// we calculate the intensity as 
// a function of the radial distance
// and the defocus deltaZ, using
// Kirchhoffs diffraction integeral 
{

	double rho, dRho=0.01;
	gsl_complex integral = gsl_complex_rect(0,0); 

	double k=2*M_PI/lambda;

	double besselJ0;
	gsl_complex OpticalPathDiff;
	gsl_complex i = gsl_complex_rect(0,1);
	gsl_complex stuffInExponent; 
	gsl_complex dInt;
	
	double zdEffective = params.zdDesign*params.NA/sqrt(pow(params.M,2)-pow(params.NA,2)); 

	for(rho=0.0;rho<=1.0;rho+=dRho)
	{
		besselJ0 = gsl_sf_bessel_J0(k*zdEffective*params.M*rho*R/params.zdReal);
		OpticalPathDiff =  OPD(params, zdEffective, rho,deltaZ);
		stuffInExponent = gsl_complex_mul_real( gsl_complex_mul( i,OpticalPathDiff ),k);
		dInt = gsl_complex_mul_real( gsl_complex_exp(stuffInExponent),besselJ0*rho*dRho );
		integral = gsl_complex_add(integral,dInt);
	}
	return 1/gsl_pow_2(params.zdReal)*gsl_complex_abs2(integral);
}

void normalizePSF(double normFactor, double *image, int sizeZ, int sizeXY)
{

	int x,y,z;
	
	for(x=0;x<sizeXY;x++)
	{
		for(y=0;y<sizeXY;y++)
		{
			for(z=0;z<sizeZ;z++)
			{
				image[TOINDEX(x,y,z)] /= normFactor;
			}
		}
	}
}

int  calculate3dPSF(microscopeParams params, 		// input 1
					long sizeZ, 					// input 2
					long sizeXY, 					// input 3
					double lambda, 					// input 4
					double oversampelingFactor, 	// input 5
					double **image, 				// output 1
					double **scales)				// output 2
{  
	double Rstep = lambda/(4*params.NA*oversampelingFactor); 
							// the distance that is represented 
							// by two neigbouring pixels 
							// in the radial direction, 
							// calculated with nyquist

	double Zstep = lambda/(2*( params.NoilReal-sqrt( gsl_pow_2(params.NoilReal) - gsl_pow_2(params.NA))) * oversampelingFactor );
	
	Zstep = fabs(Zstep);
	
	double defocus=0, R;
	long x,y,z=0;
	int calcFlag=1;
	
	double Zstart = -0.5*Zstep*sizeZ;
	
	double normFactor=0; 
				// the factor we will be dividing our PSF with in order that the integral over the 
			   	// the PSF be one

	defocus = Zstart;
	
	(*image) = (double*)malloc(sizeZ*sizeXY*sizeXY*sizeof(double));
	
	if (*image == 0)
	{
		printf("Error allocating memory!\n");
		return 1;
	}
	
	(*scales) = (double*)malloc(2*sizeof(double));
	
	if (*scales == 0)
	{
		printf("Error allocating memory!\n");
		return 1;
	}

	for(z=0;z<sizeZ;z++)
	{
		for(y=0;y<sizeXY;y++)
		{
			for(x=0;x<sizeXY;x++)
			{
				if (y > sizeXY/2 )
				{
					(*image)[TOINDEX(x,y,z)] = (*image)[TOINDEX(x,sizeXY-y,z)];
					calcFlag=0;
				}

				if (x > sizeXY/2  && calcFlag)
				{
					(*image)[TOINDEX(x,y,z)] = (*image)[TOINDEX(sizeXY-x,y,z)];
					calcFlag=0;
				}

				if (y > x && calcFlag)
				{
					(*image)[TOINDEX(x,y,z)] = (*image)[TOINDEX(y,x,z)];
					calcFlag=0;
				}

				if(calcFlag)
				{
					R = sqrt(gsl_pow_2(x-sizeXY/2)+gsl_pow_2(y-sizeXY/2)) * Rstep;
					(*image)[TOINDEX(x,y,z)] = Intensity(params, lambda, R, defocus);
				}
				calcFlag=1;
				normFactor += (*image)[TOINDEX(x,y,z)];
			}
		}
		defocus+=Zstep;
	}

	normalizePSF(normFactor, *image, sizeZ, sizeXY);
	
	(*scales)[0] = Rstep;
	(*scales)[1] = Zstep;
	
	return 0;
	
}

