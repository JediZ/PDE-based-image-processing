/*
Functions related to level-set:s

Author: Jarno Ralli
E-mail: jarno@ralli.fi
Web: www.jarnoralli.fi (or www.jarnoralli.com)

If you use this code, please reference (some) of my papers available at http://www.jarnoralli.fi

The program is delivered as it is and the author accepts no responsibility what so ever of its use and/or results.
Errors and suggestions are kindly to be communicated to the author.

(C) Copyright 2011-2014, Jarno Ralli

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "levelsetSolvers.h"

/*Maximum and minimum values that the level-set function may obtain*/
#define PMIN -5.0f
#define PMAX 5.0f
#define GRADNORM_ZERO_CHECK

#if !defined( isnan )
  #define isnan(x)((x)!=(x))
#endif

#define max(A,B) (((A)>(B))?(A):(B))
#define maxP2(A) (((A)>(0.0f))?(A*A):(0.0f))
#define min(A,B) (((A)<(B))?(A):(B))
#define minP2(A) (((A)<(0.0f))?(A*A):(0.0f))

/*Semi-implicit and AOS (Additive Operator Splitting) schemes
$ u{k+1} = \left( I - \tau \sum\limits_{l=1}^m A_l(u^k) \right)^{-1} u^k $: semi-implicit scheme
$ u^{k+1} = \frac{1}{m} \sum\limits_{l=1}^{m} \left( I - m \tau A_l(u^k) \right)^{-1} u^k: AOS scheme
*/

/*
//---------------------------------------------------------------------------------------
//---Additive Operator Splitting (AOS) scheme for Active Contour Models
//---Uses tridiagonal matrix algorithm (TDMA) aka Thomas algorithm.
//---2D (spatial) regularization with 4-neighbours
//---------------------------------------------------------------------------------------
*/
void	CV_AOS_4_2d(	struct matrixM *PHI_out,
			struct matrixM *PHI_in,
			struct matrixM *D_in,		/*Data*/
			struct matrixM *GradNorm_in,	/*Derived from N = (\nabla \Phi)/|\nabla \Phi|*/
			struct matrixM *Diff_in,	/*Diffusivity, inside DIV*/
			float tau,
			float nu )
{
	float Cp[MAX_BUF_SIZE], Dp[MAX_BUF_SIZE];
	int nrows = PHI_in->dimElems[0];
	int ncols = PHI_in->dimElems[1];
	int nframes = 1;
	int maxLength;
	int iterations = 1;
	int iter;

	if(PHI_in->ndims>2)
		nframes = PHI_in->dimElems[2];
	
	/*
	------------------------
	-Vertical line relaxing
	------------------------
	*/
	CV_TDMA_Column4( PHI_out, PHI_in, D_in, GradNorm_in, Diff_in, tau, nu, Cp, Dp, ncols, nrows, nframes );

	/*
	--------------------------
	-Horizontal line relaxing
	--------------------------
	*/
	CV_TDMA_Row4( PHI_out, PHI_in, D_in, GradNorm_in, Diff_in, tau, nu, Cp, Dp, ncols, nrows, nframes );
	
	/*Reinitialize back to a signed distance function*/
	/*reinit( PHI_out, 0.25f );*/

}

/*
//---------------------------------------------------------------------------------------
//---Additive Operator Splitting (AOS) scheme for Active Contour Models
//---Uses tridiagonal matrix algorithm (TDMA) aka Thomas algorithm.
//---2D (spatial) regularization with 4-neighbours
//--- Slightly optimized version using Open MP when solving for more than 1 level-set equation
//---------------------------------------------------------------------------------------
*/
void	CV_AOSOMP_4_2d(	struct matrixM *PHI_out,
			struct matrixM *PHI_in,
			struct matrixM *D_in,		/*Data*/
			struct matrixM *GradNorm_in,	/*Derived from N = (\nabla \Phi)/|\nabla \Phi|*/
			struct matrixM *Diff_in,	/*Diffusivity, inside DIV*/
			float tau,
			float nu )
{
	int nrows = PHI_in->dimElems[0];
	int ncols = PHI_in->dimElems[1];
	int nframes = 1;

	if(PHI_in->ndims>2)
		nframes = PHI_in->dimElems[2];
	
	/*
	------------------------
	-Vertical line relaxing
	------------------------
	*/
	CV_TDMA_Column4_omp( PHI_out, PHI_in, D_in, GradNorm_in, Diff_in, tau, nu, ncols, nrows, nframes );

	/*
	--------------------------
	-Horizontal line relaxing
	--------------------------
	*/
	CV_TDMA_Row4_omp( PHI_out, PHI_in, D_in, GradNorm_in, Diff_in, tau, nu, ncols, nrows, nframes );
	
	/*Reinitialize back to a signed distance function*/
	/*reinit( PHI_out, 0.25f );*/
}


/*
//---------------------------------------------------------------------------------------
//---Additive Operator Splitting (AOS) scheme for Active Contour Models
//---With convective derivative term
//---Uses tridiagonal matrix algorithm (TDMA) aka Thomas algorithm.
//---2D (spatial) regularization with 4-neighbours
//---------------------------------------------------------------------------------------
*/
void	AC_AOS_4_2d(	struct matrixM *PHI_out,
			struct matrixM *PHI_in,
			struct matrixM *D_in,		/*Data*/
			struct matrixM *GradNorm_in,	/*Derived from N = (nabla PHI)/|nabla PHI|*/
			struct matrixM *Diff_in,	/*Diffusivity, inside DIV*/
			float tau,			/*tau = delta t = time step*/
			float nu, 			/*diffusivity coefficient*/
			float Treinit)			/*t=0:0.25:Treinit...reinitialisation steps*/
{
	/*float *Cp, *Dp;		See below*/
	float Cp[MAX_BUF_SIZE], Dp[MAX_BUF_SIZE];
	int nrows = PHI_in->dimElems[0];
	int ncols = PHI_in->dimElems[1];
	int nframes = 1;

	if(PHI_in->ndims>2)
		nframes = PHI_in->dimElems[2];
	
	/*
	------------------------
	-Vertical line relaxing
	------------------------
	*/
	AC_TDMA_column4( PHI_out, PHI_in, D_in, GradNorm_in, Diff_in, tau, nu, Cp, Dp, ncols, nrows, nframes );

	/*
	--------------------------
	-Horizontal line relaxing
	--------------------------
	*/
	AC_TDMA_row4( PHI_out, PHI_in, D_in, GradNorm_in, Diff_in, tau, nu, Cp, Dp, ncols, nrows, nframes );
	
	/*Reinitialize back to a signed distance function*/
	reinit( PHI_out, Treinit );
	
}

/*
//-------------------------------------------------------------------
//--- TDMA vertical line solving - 4 neighbours - harmonic averaging
//--- Slightly optimized version using Open MP
//-------------------------------------------------------------------
*/
void CV_TDMA_Column4_omp(	struct matrixM *PHI_out, 
				struct matrixM *PHI_in, 
				struct matrixM *D_in, 
				struct matrixM *GradNorm_in, 
				struct matrixM *Diff_in, 
				float tau, 
				float nu, 
				int ncols, 
				int nrows, 
				int nframes)
{
omp_set_dynamic(1);
#pragma omp parallel
{
	int k, frameoffset=ncols*nrows, myFrames=nframes;

	/*tau = 4.0f*tau;	/*This comes from the AOS scheme: u^{t+1} = \sum_{l=1}^{m} ( mI - m^{2} \tau A_{l}(u^t) )^{-1} u^t*/
	
	#pragma omp for
	for(k=0;k<myFrames;k++)
	{
	/*mexPrintf("Num threads %d, thread ID %d.\n", omp_get_num_threads(), omp_get_thread_num());*/
	int i,j,pos;
	float a,b,c,d,div,Diffp,Diffn,temp1,temp2;
	float cp[MAX_BUF_SIZE], dp[MAX_BUF_SIZE];
	  
	for(j=0;j<=ncols-1;j++)
	{
		/*---FORWARD SWEEP---*/
		/*---First row---*/
		i = 0;
		pos = j*nrows+k*frameoffset;
		/*Harmonic averaging - next element*/
		Diffn = ( (temp1=Diff_in->data[pos]+Diff_in->data[pos+1])>0.0f )?(2*tau*GradNorm_in->data[pos]/temp1):(0.0f);
		
		b = 	2.0f + nu*Diffn;
		cp[i] = -nu*Diffn/b;
		dp[i] = (PHI_in->data[pos] + tau*GradNorm_in->data[pos]*D_in->data[pos])/b;

		/*---Middle rows---*/
		for(i=1;i<=nrows-2;i++)
		{
			pos++;
			/*Harmonic averaging - next and previous elements*/
			Diffn = ( (temp1=Diff_in->data[pos]+Diff_in->data[pos+1])>0.0f )?(2*tau*GradNorm_in->data[pos]/temp1):(0.0f);
			Diffp = ( (temp1=Diff_in->data[pos]+Diff_in->data[pos-1])>0.0f )?(2*tau*GradNorm_in->data[pos]/temp1):(0.0f);

			a = 	-nu*Diffp;
			div = 1 / (2.0f + nu*(Diffn+Diffp) - cp[i-1]*a);
			cp[i] = -nu*Diffn*div;
			dp[i] = (PHI_in->data[pos] + tau*GradNorm_in->data[pos]*D_in->data[pos] - dp[i-1]*a)*div;
		}
		/*---Last row---*/
		pos++;
		/*Harmonic averaging - previous element*/
		Diffp = ( (temp1=Diff_in->data[pos]+Diff_in->data[pos-1])>0.0f )?(2*tau*GradNorm_in->data[pos]/temp1):(0.0f);
		
		a = 	-nu*Diffp;
		dp[i] = (PHI_in->data[pos] + tau*GradNorm_in->data[pos]*D_in->data[pos]-dp[i-1]*a) / (2.0f + nu*Diffp-cp[i-1]*a);
		
		#if defined(GRADNORM_ZERO_CHECK)
			/*---BACKWARD SWEEP---*/
			temp1 = dp[i];
			
			for(i=nrows-2;i>=0;i--)
			{
				pos--;

				temp2 = dp[i] - cp[i]*temp1;

				/*TODO: This is not very graceful!!!*/
				if(Diff_in->data[pos+1]!=0.0f)
					PHI_out->data[pos+1] += temp1;
				else
					PHI_out->data[pos+1] = PHI_in->data[pos+1];

				temp1 = temp2;	
				/*Limit the results*/
				if(PHI_out->data[pos+1]>PMAX) PHI_out->data[pos+1] = PMAX;
				if(PHI_out->data[pos+1]<PMIN) PHI_out->data[pos+1] = PMIN;
			}
			PHI_out->data[pos] += temp1;
			/*Limit the results*/
			if(PHI_out->data[pos]>PMAX) PHI_out->data[pos] = PMAX;
			if(PHI_out->data[pos]<PMIN) PHI_out->data[pos] = PMIN;
		#else
			/*---BACKWARD SWEEP---*/
			temp1 = dp[i];
			
			for(i=nrows-2;i>=0;i--)
			{
				pos--;
				
				temp2 = dp[i] - cp[i]*temp1;
				PHI_out->data[pos+1] += temp1;
				temp1 = temp2;
				/*Limit the results*/
				if(PHI_out->data[pos+1]>PMAX) PHI_out->data[pos+1] = PMAX;
				if(PHI_out->data[pos+1]<PMIN) PHI_out->data[pos+1] = PMIN;
			}
			PHI_out->data[pos] += temp1;
			/*Limit the results*/
			if(PHI_out->data[pos]>PMAX) PHI_out->data[pos] = PMAX;
			if(PHI_out->data[pos]<PMIN) PHI_out->data[pos] = PMIN;
		#endif
	}
	}
}
}

/*
//---------------------------------------------------------------------
//--- TDMA horizontal line solving - 4 neighbours - harmonic averaging
//--- Slightly optimized version using Open MP
//---------------------------------------------------------------------
*/
void CV_TDMA_Row4_omp(	struct matrixM *PHI_out, 
			struct matrixM *PHI_in, 
			struct matrixM *D_in, 
			struct matrixM *GradNorm_in, 
			struct matrixM *Diff_in, 
			float tau, 
			float nu, 
			int ncols, 
			int nrows, 
			int nframes)
{
omp_set_dynamic(1);
#pragma omp parallel
{
	int k,frameoffset=ncols*nrows, myFrames=nframes;;

	/*tau = 4.0f*tau;	/*This comes from the AOS scheme: u^{t+1} = \sum_{l=1}^{m} ( mI - m^{2} \tau A_{l}(u^t) )^{-1} u^t*/
	
	#pragma omp for
	for(k=0;k<myFrames;k++)
	{
	int i,j,pos;
	float a,b,c,d,div,Diffp,Diffn,temp1,temp2;
	float cp[MAX_BUF_SIZE], dp[MAX_BUF_SIZE];
	
	for(i=0;i<=nrows-1;i++)
	{
		/*---FORWARD SWEEP---*/
		/*---First column---*/
		j = 0;
		pos = j*nrows+i+k*frameoffset;
		/*Harmonic averaging - next element*/
		Diffn = ( (temp1=Diff_in->data[pos]+Diff_in->data[pos+nrows])>0.0f )?(2*tau*GradNorm_in->data[pos]/temp1):(0.0f);
	
		/*a = 	-nu*Diffp*/
		b = 	2.0f + nu*(Diffn);
		cp[j] = -nu*Diffn/b;
		dp[j] = (PHI_in->data[pos] + tau*GradNorm_in->data[pos]*D_in->data[pos])/b;

		/*---Middle columns---*/
		for(j=1;j<=ncols-2;j++)
		{
			pos += nrows;
			/*Harmonic averaging - next and previous elements*/
			Diffn = ( (temp1=Diff_in->data[pos]+Diff_in->data[pos+nrows])>0.0f )?(2*tau*GradNorm_in->data[pos]/temp1):(0.0f);
			Diffp = ( (temp1=Diff_in->data[pos]+Diff_in->data[pos-nrows])>0.0f )?(2*tau*GradNorm_in->data[pos]/temp1):(0.0f);
	
			a = 	-nu*Diffp;
			div = 1 / (2.0f + nu*(Diffn+Diffp) - cp[j-1]*a);
			cp[j] = -nu*Diffn*div;
			dp[j] = (PHI_in->data[pos] + tau*GradNorm_in->data[pos]*D_in->data[pos] - dp[j-1]*a)*div;
		}

		/*---Last column---*/
		pos += nrows;
		/*Harmonic averaging - previous element*/
		Diffp = ( (temp1=Diff_in->data[pos]+Diff_in->data[pos-nrows])>0.0f )?(2*tau*GradNorm_in->data[pos]/temp1):(0.0f);
	
		a = 	-nu*Diffp;
		b = 	2.0f + nu*(Diffp);
		/*c = 	-nu*Diffn*/
		d = 	PHI_in->data[pos] + tau*GradNorm_in->data[pos]*D_in->data[pos];
		
		dp[j] = (d-dp[j-1]*a) / (b-cp[j-1]*a);

		/*---BACKWARD SWEEP---*/
		#if defined(GRADNORM_ZERO_CHECK)
			/*---BACKWARD SWEEP---*/
			temp1 = dp[j];
			
			for(j=ncols-2;j>=0;j--)
			{
				pos -= nrows;

				temp2 = dp[j] - cp[j]*temp1;

				/*TODO: This is not very graceful!!!*/
				if(Diff_in->data[pos+nrows]!=0.0f)
					PHI_out->data[pos+nrows] += temp1;
				else
					PHI_out->data[pos+nrows] = PHI_in->data[pos+nrows];

				temp1 = temp2;	
				/*Limit the results*/
				if(PHI_out->data[pos+nrows]>PMAX) PHI_out->data[pos+nrows] = PMAX;
				if(PHI_out->data[pos+nrows]<PMIN) PHI_out->data[pos+nrows] = PMIN;
			}
			PHI_out->data[pos] += temp1;
			/*Limit the results*/
			if(PHI_out->data[pos]>PMAX) PHI_out->data[pos] = PMAX;
			if(PHI_out->data[pos]<PMIN) PHI_out->data[pos] = PMIN;
		#else
			temp1 = PHI_out->data[pos];
			PHI_out->data[pos] =  dp[j];
			/*Limit the results*/
			if(PHI_out->data[pos]>PMAX) PHI_out->data[pos] = PMAX;
			if(PHI_out->data[pos]<PMIN) PHI_out->data[pos] = PMIN;

			for(j=ncols-2;j>=0;j--)
			{
				pos -= nrows;

				temp2 = PHI_out->data[pos];
				PHI_out->data[pos] = dp[j] - cp[j]*PHI_out->data[pos+nrows];
				PHI_out->data[pos+nrows] += temp1;
				temp1 = temp2;
				/*Limit the results*/
				if(PHI_out->data[pos+nrows]>PMAX) PHI_out->data[pos+nrows] = PMAX;
				if(PHI_out->data[pos+nrows]<PMIN) PHI_out->data[pos+nrows] = PMIN;
			}
			PHI_out->data[pos] += temp1;
			/*Limit the results*/
			if(PHI_out->data[pos]>PMAX) PHI_out->data[pos] = PMAX;
			if(PHI_out->data[pos]<PMIN) PHI_out->data[pos] = PMIN;
		#endif
	}
	}
}
}

/*
//-------------------------------------------------------------------
//--- TDMA vertical line solving - 4 neighbours - harmonic averaging
//-------------------------------------------------------------------
*/
void CV_TDMA_Column4(		struct matrixM *PHI_out, 
				struct matrixM *PHI_in, 
				struct matrixM *D_in, 
				struct matrixM *GradNorm_in, 
				struct matrixM *Diff_in, 
				float tau, 
				float nu, 
				float *cp,
				float *dp,
				int ncols, 
				int nrows, 
				int nframes)
{
	int i,j,k, pos, frameoffset=ncols*nrows;
	float a,b,c,d,div,Diffp,Diffn,temp1,temp2;

	/*tau = 4.0f*tau;	/*This comes from the AOS scheme: u^{t+1} = \sum_{l=1}^{m} ( mI - m^{2} \tau A_{l}(u^t) )^{-1} u^t*/

	for(k=0;k<nframes;k++)
	{
	for(j=0;j<=ncols-1;j++)
	{
		/*---FORWARD SWEEP---*/
		/*---First row---*/
		i = 0;
		pos = j*nrows+k*frameoffset;
		/*Harmonic averaging - next element*/
		Diffn = ( (temp1=Diff_in->data[pos]+Diff_in->data[pos+1])>0.0f )?(2*tau*GradNorm_in->data[pos]/temp1):(0.0f);
		
		/*a = 	0.0f;*/
		b = 	2.0f + nu*Diffn;
		c = 	-nu*Diffn;
		d = 	PHI_in->data[pos] + tau*GradNorm_in->data[pos]*D_in->data[pos];
			
		cp[i] = c/b;
		dp[i] = d/b;

		/*---Middle rows---*/
		for(i=1;i<=nrows-2;i++)
		{
			pos++;
			/*Harmonic averaging - next and previous elements*/
			Diffn = ( (temp1=Diff_in->data[pos]+Diff_in->data[pos+1])>0.0f )?(2*tau*GradNorm_in->data[pos]/temp1):(0.0f);
			Diffp = ( (temp1=Diff_in->data[pos]+Diff_in->data[pos-1])>0.0f )?(2*tau*GradNorm_in->data[pos]/temp1):(0.0f);

			a = 	-nu*Diffp;
			b = 	2.0f + nu*(Diffn+Diffp);
			c = 	-nu*Diffn;
			d = 	PHI_in->data[pos] + tau*GradNorm_in->data[pos]*D_in->data[pos];
		
			div = 1 / (b-cp[i-1]*a);
			cp[i] = c*div;
			dp[i] = (d-dp[i-1]*a)*div;
		}
		/*---Last row---*/
		pos++;
		/*Harmonic averaging - previous element*/
		Diffp = ( (temp1=Diff_in->data[pos]+Diff_in->data[pos-1])>0.0f )?(2*tau*GradNorm_in->data[pos]/temp1):(0.0f);
		
		a = 	-nu*Diffp;
		b = 	2.0f + nu*Diffp;
		/*c = 	-nu*Diffn*/
		d = 	PHI_in->data[pos] + tau*GradNorm_in->data[pos]*D_in->data[pos];
		
		dp[i] = (d-dp[i-1]*a) / (b-cp[i-1]*a);
		
		#if defined(GRADNORM_ZERO_CHECK)
			/*---BACKWARD SWEEP---*/
			temp1 = dp[i];
			
			for(i=nrows-2;i>=0;i--)
			{
				pos--;

				temp2 = dp[i] - cp[i]*temp1;

				/*TODO: This is not very graceful!!!*/
				if(Diff_in->data[pos+1]!=0.0f)
					PHI_out->data[pos+1] += temp1;
				else
					PHI_out->data[pos+1] = PHI_in->data[pos+1];

				temp1 = temp2;	
				/*Limit the results*/
				if(PHI_out->data[pos+1]>PMAX) PHI_out->data[pos+1] = PMAX;
				if(PHI_out->data[pos+1]<PMIN) PHI_out->data[pos+1] = PMIN;
			}
			PHI_out->data[pos] += temp1;
			/*Limit the results*/
			if(PHI_out->data[pos]>PMAX) PHI_out->data[pos] = PMAX;
			if(PHI_out->data[pos]<PMIN) PHI_out->data[pos] = PMIN;
		#else
			/*---BACKWARD SWEEP---*/
			temp1 = dp[i];
			
			for(i=nrows-2;i>=0;i--)
			{
				pos--;
				
				temp2 = dp[i] - cp[i]*temp1;
				PHI_out->data[pos+1] += temp1;
				temp1 = temp2;
				/*Limit the results*/
				if(PHI_out->data[pos+1]>PMAX) PHI_out->data[pos+1] = PMAX;
				if(PHI_out->data[pos+1]<PMIN) PHI_out->data[pos+1] = PMIN;
			}
			PHI_out->data[pos] += temp1;
			/*Limit the results*/
			if(PHI_out->data[pos]>PMAX) PHI_out->data[pos] = PMAX;
			if(PHI_out->data[pos]<PMIN) PHI_out->data[pos] = PMIN;
		#endif
	}
	}
}

/*
//---------------------------------------------------------------------
//--- TDMA horizontal line solving - 4 neighbours - harmonic averaging
//---------------------------------------------------------------------
*/
void CV_TDMA_Row4(	struct matrixM *PHI_out, 
			struct matrixM *PHI_in, 
			struct matrixM *D_in, 
			struct matrixM *GradNorm_in, 
			struct matrixM *Diff_in, 
			float tau, 
			float nu, 
			float *cp,
			float *dp,
			int ncols, 
			int nrows, 
			int nframes)
{
	int i,j,k,pos,frameoffset=ncols*nrows;
	float a,b,c,d,div,Diffn,Diffp,temp1,temp2;

	/*tau = 4.0f*tau;	/*This comes from the AOS scheme: u^{t+1} = \sum_{l=1}^{m} ( mI - m^{2} \tau A_{l}(u^t) )^{-1} u^t*/
	
	for(k=0;k<nframes;k++)
	{
	for(i=0;i<=nrows-1;i++)
	{
		/*---FORWARD SWEEP---*/
		/*---First column---*/
		j = 0;
		pos = j*nrows+i+k*frameoffset;
		/*Harmonic averaging - next element*/
		Diffn = ( (temp1=Diff_in->data[pos]+Diff_in->data[pos+nrows])>0.0f )?(2*tau*GradNorm_in->data[pos]/temp1):(0.0f);
	
		/*a = 	-nu*Diffp*/
		b = 	2.0f + nu*(Diffn);
		c = 	-nu*Diffn;
		d = 	PHI_in->data[pos] + tau*GradNorm_in->data[pos]*D_in->data[pos];
				
		cp[j] = c/b;
		dp[j] = d/b;

		/*---Middle columns---*/
		for(j=1;j<=ncols-2;j++)
		{
			pos+=nrows;
			/*Harmonic averaging - next and previous elements*/
			Diffn = ( (temp1=Diff_in->data[pos]+Diff_in->data[pos+nrows])>0.0f )?(2*tau*GradNorm_in->data[pos]/temp1):(0.0f);
			Diffp = ( (temp1=Diff_in->data[pos]+Diff_in->data[pos-nrows])>0.0f )?(2*tau*GradNorm_in->data[pos]/temp1):(0.0f);
	
			a = 	-nu*Diffp;
			b = 	2.0f + nu*(Diffn+Diffp);
			c = 	-nu*Diffn;
			d = 	PHI_in->data[pos] + tau*GradNorm_in->data[pos]*D_in->data[pos];
				
			div = 1 / (b-cp[j-1]*a);
			cp[j] = c*div;
			dp[j] = (d-dp[j-1]*a)*div;
		}

		/*---Last column---*/
		pos+=nrows;
		/*Harmonic averaging - previous element*/
		Diffp = ( (temp1=Diff_in->data[pos]+Diff_in->data[pos-nrows])>0.0f )?(2*tau*GradNorm_in->data[pos]/temp1):(0.0f);
	
		a = 	-nu*Diffp;
		b = 	2.0f + nu*(Diffp);
		/*c = 	-nu*Diffn*/
		d = 	PHI_in->data[pos] + tau*GradNorm_in->data[pos]*D_in->data[pos];
		
		dp[j] = (d-dp[j-1]*a) / (b-cp[j-1]*a);

		/*---BACKWARD SWEEP---*/
		#if defined(GRADNORM_ZERO_CHECK)
			/*---BACKWARD SWEEP---*/
			temp1 = dp[j];
			
			for(j=ncols-2;j>=0;j--)
			{
				pos -= nrows;

				temp2 = dp[j] - cp[j]*temp1;

				/*TODO: This is not very graceful!!!*/
				if(Diff_in->data[pos+nrows]!=0.0f)
					PHI_out->data[pos+nrows] += temp1;
				else
					PHI_out->data[pos+nrows] = PHI_in->data[pos+nrows];

				temp1 = temp2;	
				/*Limit the results*/
				if(PHI_out->data[pos+nrows]>PMAX) PHI_out->data[pos+nrows] = PMAX;
				if(PHI_out->data[pos+nrows]<PMIN) PHI_out->data[pos+nrows] = PMIN;
			}
			PHI_out->data[pos] += temp1;
			/*Limit the results*/
			if(PHI_out->data[pos]>PMAX) PHI_out->data[pos] = PMAX;
			if(PHI_out->data[pos]<PMIN) PHI_out->data[pos] = PMIN;
		#else
			temp1 = PHI_out->data[pos];
			PHI_out->data[pos] =  dp[j];
			/*Limit the results*/
			if(PHI_out->data[pos]>PMAX) PHI_out->data[pos] = PMAX;
			if(PHI_out->data[pos]<PMIN) PHI_out->data[pos] = PMIN;

			for(j=ncols-2;j>=0;j--)
			{
				pos -= nrows;

				temp2 = PHI_out->data[pos];
				PHI_out->data[pos] = dp[j] - cp[j]*PHI_out->data[pos+nrows];
				PHI_out->data[pos+nrows] += temp1;
				temp1 = temp2;
				/*Limit the results*/
				if(PHI_out->data[pos+nrows]>PMAX) PHI_out->data[pos+nrows] = PMAX;
				if(PHI_out->data[pos+nrows]<PMIN) PHI_out->data[pos+nrows] = PMIN;
			}
			PHI_out->data[pos] += temp1;
			/*Limit the results*/
			if(PHI_out->data[pos]>PMAX) PHI_out->data[pos] = PMAX;
			if(PHI_out->data[pos]<PMIN) PHI_out->data[pos] = PMIN;
		#endif
	}
	}
}

/*
//------------------------------------------------------------------------------------------
//--- TDMA vertical line solving - 4 neighbours - harmonic averaging - active contour model
//------------------------------------------------------------------------------------------
*/
void AC_TDMA_column4(	struct matrixM *PHI_out, 
			struct matrixM *PHI_in, 
			struct matrixM *D_in, 
			struct matrixM *GradNorm_in, 
			struct matrixM *Diff_in, 
			float tau, 
			float nu, 
			float *cp,
			float *dp,
			int ncols, 
			int nrows, 
			int nframes)
{
	int i,j,k, pos, frameoffset=ncols*nrows;
	float a,b,c,d,div,Diffp,Diffn,temp1,temp2;

	/*tau = 4.0f*tau;	/*This comes from the AOS scheme: u^{t+1} = \sum_{l=1}^{m} ( mI - m^{2} \tau A_{l}(u^t) )^{-1} u^t*/

	for(k=0;k<nframes;k++)
	{
	for(j=0;j<=ncols-1;j++)
	{
		/*---FORWARD SWEEP---*/
		/*---First row---*/
		i = 0;
		pos = j*nrows+k*frameoffset;
		Diffn = ( (temp1=Diff_in->data[pos]+Diff_in->data[pos+1])>0.0f )?(2*tau*GradNorm_in->data[pos]/temp1):(0.0f);
		
		/*a = 	0.0f;*/
		b = 	2.0f + nu*Diffn;
		c = 	-nu*Diffn;
		d = 	PHI_in->data[pos] + tau*D_in->data[pos];
			
		cp[i] = c/b;
		dp[i] = d/b;

		/*---Middle rows---*/
		for(i=1;i<=nrows-2;i++)
		{
			pos = j*nrows+i+k*frameoffset;
			Diffn = ( (temp1=Diff_in->data[pos]+Diff_in->data[pos+1])>0.0f )?(2*tau*GradNorm_in->data[pos]/temp1):(0.0f);
			Diffp = ( (temp1=Diff_in->data[pos]+Diff_in->data[pos-1])>0.0f )?(2*tau*GradNorm_in->data[pos]/temp1):(0.0f);

			a = 	-nu*Diffp;
			b = 	2.0f + nu*(Diffn+Diffp);
			c = 	-nu*Diffn;
			d = 	PHI_in->data[pos] + tau*D_in->data[pos];
		
			div = 1 / (b-cp[i-1]*a);
			cp[i] = c*div;
			dp[i] = (d-dp[i-1]*a)*div;
		}
		/*---Last row---*/
		pos = j*nrows+i+k*frameoffset;
		Diffp = ( (temp1=Diff_in->data[pos]+Diff_in->data[pos-1])>0.0f )?(2*tau*GradNorm_in->data[pos]/temp1):(0.0f);
		
		a = 	-nu*Diffp;
		b = 	2.0f + nu*Diffp;
		/*c = 	-nu*Diffn*/
		d = 	PHI_in->data[pos] + tau*D_in->data[pos];
		
		dp[i] = (d-dp[i-1]*a) / (b-cp[i-1]*a);
		
		/*---BACKWARD SWEEP---*/
		temp1 = PHI_out->data[pos];
		PHI_out->data[pos] =  dp[i];
		/*Limit the results*/
		/*if(PHI_out->data[pos]>PMAX) PHI_out->data[pos] = PMAX;
		if(PHI_out->data[pos]<PMIN) PHI_out->data[pos] = PMIN;*/
		
		for(i=nrows-2;i>=0;i--)
		{
			pos = j*nrows+i+k*frameoffset;
			

			temp2 = PHI_out->data[pos];
			/*TODO: This is not very graceful!!!*/
			if(Diff_in->data[pos]!=0.0f)
				PHI_out->data[pos] = dp[i] - cp[i]*PHI_out->data[pos+1];
			else
			{
				PHI_out->data[pos] = PHI_in->data[pos];
				temp2=0.0f;
			}
				PHI_out->data[pos+1] += temp1;
				temp1 = temp2;
				/*Limit the results*/
				/*if(PHI_out->data[pos+1]>PMAX) PHI_out->data[pos+1] = PMAX;
				if(PHI_out->data[pos+1]<PMIN) PHI_out->data[pos+1] = PMIN;*/
		}
		PHI_out->data[pos] += temp1;
		/*Limit the results*/
		/*if(PHI_out->data[pos]>PMAX) PHI_out->data[pos] = PMAX;
		if(PHI_out->data[pos]<PMIN) PHI_out->data[pos] = PMIN;*/
	}
	}
}

/*
//--------------------------------------------------------------------------------------------
//--- TDMA horizontal line solving - 4 neighbours - harmonic averaging - active contour model
//--------------------------------------------------------------------------------------------
*/
void AC_TDMA_row4(	struct matrixM *PHI_out, 
			struct matrixM *PHI_in, 
			struct matrixM *D_in, 
			struct matrixM *GradNorm_in, 
			struct matrixM *Diff_in, 
			float tau, 
			float nu, 
			float *cp,
			float *dp,
			int ncols, 
			int nrows, 
			int nframes)
{
	int i,j,k,pos,frameoffset=ncols*nrows;
	float a,b,c,d,div,Diffn,Diffp,temp1,temp2;

	/*tau = 4.0f*tau;	/*This comes from the AOS scheme: u^{t+1} = \sum_{l=1}^{m} ( mI - m^{2} \tau A_{l}(u^t) )^{-1} u^t*/
	
	for(k=0;k<nframes;k++)
	{
	for(i=0;i<=nrows-1;i++)
	{
		/*---FORWARD SWEEP---*/
		/*---First column---*/
		j = 0;
		pos = j*nrows+i+k*frameoffset;
		Diffn = ( (temp1=Diff_in->data[pos]+Diff_in->data[pos+nrows])>0.0f )?(2*tau*GradNorm_in->data[pos]/temp1):(0.0f);
	
		/*a = 	-nu*Diffp*/
		b = 	2.0f + nu*(Diffn);
		c = 	-nu*Diffn;
		d = 	PHI_in->data[pos] + tau*D_in->data[pos];
				
		cp[j] = c/b;
		dp[j] = d/b;

		/*---Middle columns---*/
		for(j=1;j<=ncols-2;j++)
		{
			pos = j*nrows+i+k*frameoffset;
			Diffn = ( (temp1=Diff_in->data[pos]+Diff_in->data[pos+nrows])>0.0f )?(2*tau*GradNorm_in->data[pos]/temp1):(0.0f);
			Diffp = ( (temp1=Diff_in->data[pos]+Diff_in->data[pos-nrows])>0.0f )?(2*tau*GradNorm_in->data[pos]/temp1):(0.0f);
	
			a = 	-nu*Diffp;
			b = 	2.0f + nu*(Diffn+Diffp);
			c = 	-nu*Diffn;
			d = 	PHI_in->data[pos] + tau*D_in->data[pos];
				
			div = 1 / (b-cp[j-1]*a);
			cp[j] = c*div;
			dp[j] = (d-dp[j-1]*a)*div;
		}

		/*---Last column---*/
		pos = j*nrows+i+k*frameoffset;
		Diffp = ( (temp1=Diff_in->data[pos]+Diff_in->data[pos-nrows])>0.0f )?(2*tau*GradNorm_in->data[pos]/temp1):(0.0f);
	
		a = 	-nu*Diffp;
		b = 	2.0f + nu*(Diffp);
		/*c = 	-nu*Diffn*/
		d = 	PHI_in->data[pos] + tau*D_in->data[pos];
		
		dp[j] = (d-dp[j-1]*a) / (b-cp[j-1]*a);

		/*---BACKWARD SWEEP---*/
		temp1 = PHI_out->data[pos];
		PHI_out->data[pos] = dp[j];
		/*Limit the results*/
		/*if(PHI_out->data[pos]>PMAX) PHI_out->data[pos] = PMAX;
		if(PHI_out->data[pos]<PMIN) PHI_out->data[pos] = PMIN;*/

		for(j=ncols-2;j>=0;j--)
		{
			pos = j*nrows+i+k*frameoffset;
			temp2 = PHI_out->data[pos];

			/*TODO: This is not very graceful!!!*/
			if(Diff_in->data[pos]!=0.0f)
				PHI_out->data[pos] = dp[j] - cp[j]*PHI_out->data[pos+nrows];
			else
			{
				temp2 = 0.0f;
				PHI_out->data[pos+nrows] = PHI_in->data[pos+nrows];
			}
				PHI_out->data[pos+nrows] += temp1;
				temp1 = temp2;
				/*Limit the results*/
				/*if(PHI_out->data[pos+nrows]>PMAX) PHI_out->data[pos+nrows] = PMAX;
				if(PHI_out->data[pos+nrows]<PMIN) PHI_out->data[pos+nrows] = PMIN;*/
		}
		PHI_out->data[pos] += temp1;
		/*Limit the results*/
		/*if(PHI_out->data[pos]>PMAX) PHI_out->data[pos] = PMAX;
		if(PHI_out->data[pos]<PMIN) PHI_out->data[pos] = PMIN;*/
	}
	}
}

/*
//-------------------------
//--- Vertical convolution 
//-------------------------
*/
void VerticalConv(float *Result, float *I, float operator[2], unsigned int rows, unsigned int cols, unsigned int frames)
{
omp_set_dynamic(1);
#pragma omp parallel
{
	int i,j,k;
	int pos, ppos, npos;
	float temp_result[2];

	#pragma omp for
	for(k=0;k<frames;k++)
	{
	for(j=0;j<cols;j++)
	{
		/*position, previous position and next position*/
		pos = k*rows*cols+j*rows; ppos = pos-1; npos = pos+1;
		temp_result[0] = I[pos]*operator[0];
		temp_result[1] = I[pos+1]*operator[1];
		Result[pos] = temp_result[0]+temp_result[1];

		for(i=0;i<rows-2;i++)
		{
			/*position, previous position and next position*/
			pos++; ppos++; npos++;
			temp_result[0] = I[ppos]*operator[0];
			temp_result[1] = I[npos]*operator[1];
			Result[pos] = temp_result[0]+temp_result[1];
			
		}
		/*position, previous position and next position*/
		pos++; ppos++; npos++;
		temp_result[0] = I[ppos]*operator[0];
		temp_result[1] = I[pos]*operator[1];
		Result[pos] = temp_result[0]+temp_result[1];
	}
	}
}
}
/*
//----------------------------
//--- Horizontal convolution 
//----------------------------
*/
void HorizontalConv(float *Result, float *I, float operator[2], unsigned int rows, unsigned int cols, unsigned int frames)
{
omp_set_dynamic(1);
#pragma omp parallel
{
	int i,j,k;
	int pos, ppos, npos;
	float temp_result[2];

	#pragma omp for
	for(k=0;k<frames;k++)
	{
		for(i=0;i<rows;i++)
		{
			/*position, previous position and next position*/
			pos = k*rows*cols+i; ppos = pos-rows; npos = pos+rows;
			temp_result[0] = I[pos]*operator[0];
			temp_result[1] = I[npos]*operator[1];
			Result[pos] = temp_result[0] + temp_result[1];

			for(j=0;j<cols-2;j++)
			{	
				/*position, previous position and next position*/
				pos+=rows; ppos+=rows; npos+=rows;
				temp_result[0] = I[ppos]*operator[0];
				temp_result[1] = I[npos]*operator[1];
				Result[pos] = temp_result[0] + temp_result[1];
			}
			/*position and previous position*/
			pos+=rows; ppos+=rows; npos+=rows;
			temp_result[0] = I[ppos]*operator[0];
			temp_result[1] = I[pos]*operator[1];
			Result[pos] = temp_result[0] + temp_result[1];
		}
	}
}
}

/*
//---------------------------------------------------------------------------------------
//---Reinitialization of the level-set function to a signed distance function
//---TODO: this code is both inefficient and ugly...fix it!!!
//---------------------------------------------------------------------------------------
*/
void	reinit(	struct matrixM *PHI_inout,
		float T)
{
		float *PHIx, *PHIy;	/*Normal derivatives*/
		float *PHIxuw, *PHIyuw;	/*Upwind derivatives*/
		float *S;		/*Blurred sign function*/
		unsigned int nrows = PHI_inout->dimElems[0];
		unsigned int ncols = PHI_inout->dimElems[1];
		int nframes = 1;
		unsigned int elems;
		int i;
		float operator[]={-0.5f,0.5f};
		float t;
		
#if defined(SSE)
		__m128 *PHI_data = (__m128 *)PHI_inout->data, tempm, constm;
#endif		

		/*Number of 'frames'*/
		if(PHI_inout->ndims>2)
			nframes = PHI_inout->dimElems[2];
		elems = nrows*ncols*nframes;
		
		/*Reinitialisation equation: \Phi_t + S( \Phi_0 )( | \nabla \Phi | -1 ) = 0*/
		
#if !defined(SSE)
	/*Allocate space for derivatives and the blurred sign function*/
	if( (PHIx = (float*)calloc( elems, sizeof(float) ))==NULL )
	{
		printf("Error reserving space for 'PHIx'\n");
		return;
	}
	if( (PHIy = (float*)calloc( elems, sizeof(float) ))==NULL )
	{
		printf("Error reserving space for 'PHIy'\n");
		return;
	}
	if( (PHIxuw = (float*)calloc( elems, sizeof(float) ))==NULL )
	{
		printf("Error reserving space for 'PHIxuw'\n");
		return;
	}
	if( (PHIyuw = (float*)calloc( elems, sizeof(float) ))==NULL )
	{
		printf("Error reserving space for 'PHIyuw'\n");
		return;
	}
	if( (S = (float*)calloc( elems, sizeof(float) ))==NULL )
	{
		printf("Error reserving space for 'S'\n");
		return;
	}

	/*Time step(s)*/
	for(t=0.0f;t<T;t+=0.25f)
	{
		HorizontalConv(	PHIx, PHI_inout->data, operator, nrows, ncols, nframes);
		VerticalConv( PHIy, PHI_inout->data, operator, nrows, ncols, nframes);

		/*Implementation depends on if SSE and/or GASM is defined (in level_sets_c.h)...see the code below*/
		blurredSignFunction( S, PHI_inout->data, PHIx, PHIy, nrows, ncols, nframes);
		godunovUpwind( PHIxuw, PHIyuw, S, PHI_inout->data, nrows, ncols, nframes);

		#pragma omp parallel for
		for(i=0;i<elems;i++)
			PHI_inout->data[i] = PHI_inout->data[i] + 0.25f*( S[i]-S[i]*sqrt(PHIxuw[i] + PHIyuw[i]) );*/
	}
#else
	constm = _mm_set_ps1( 0.25f );

	/*Allocate space for derivatives and the blurred sign function*/
	if( (PHIx = (float*)_mm_malloc( elems*sizeof(float), 16 ))==NULL )
	{
		printf("Error reserving space for 'PHIx'\n");
		return;
	}
	if( (PHIy = (float*)_mm_malloc( elems*sizeof(float), 16 ))==NULL )
	{
		printf("Error reserving space for 'PHIy'\n");
		return;
	}
	if( (PHIxuw = (float*)_mm_malloc( elems*sizeof(float), 16 ))==NULL )
	{
		printf("Error reserving space for 'PHIxuw'\n");
		return;
	}
	if( (PHIyuw = (float*)_mm_malloc( elems*sizeof(float), 16 ))==NULL )
	{
		printf("Error reserving space for 'PHIyuw'\n");
		return;
	}
	if( (S = (float*)_mm_malloc( elems*sizeof(float), 16 ))==NULL )
	{
		printf("Error reserving space for 'S'\n");
		return;
	}
	
	/*Time step(s)*/
	for(t=0.0f;t<T;t+=0.25f)
	{
		HorizontalConv(	PHIx, PHI_inout->data, operator, nrows, ncols, nframes);
		VerticalConv( PHIy, PHI_inout->data, operator, nrows, ncols, nframes);

		/*Implementation depends on if SSE and/or GASM is defined (in level_sets_c.h)...see the code below*/
		blurredSignFunction( S, PHI_inout->data, PHIx, PHIy, nrows, ncols, nframes);
		godunovUpwind( PHIxuw, PHIyuw, S, PHI_inout->data, nrows, ncols, nframes);

		for(i=0;i+3<elems;i+=4)
		{
			  tempm = _mm_add_ps( *(__m128*)&PHIxuw[i], *(__m128*)&PHIyuw[i] );	/*PHIxuw[i] + PHIyuw[i]*/
			  tempm = _mm_sqrt_ps( tempm );						/*sqrt(PHIxuw[i] + PHIyuw[i])*/
			  tempm = _mm_mul_ps( tempm, *(__m128*)&S[i] );				/*S[i]*sqrt(PHIxuw[i] + PHIyuw[i])*/
			  tempm = _mm_sub_ps( *(__m128*)&S[i], tempm );				/*S[i]-S[i]*sqrt(PHIxuw[i] + PHIyuw[i])*/
			  tempm = _mm_mul_ps( constm, tempm );					/*0.25( S[i]-S[i]*sqrt(PHIxuw[i] + PHIyuw[i]) )*/
			  *(__m128*)&PHI_inout->data[i] = _mm_add_ps( tempm, *(__m128*)&PHI_inout->data[i] );
		}
		for(;i<elems;i++)
		{
			PHI_inout->data[i] = PHI_inout->data[i] + 0.25f*( S[i]-S[i]*sqrt(PHIxuw[i] + PHIyuw[i]) );
		}

	}
#endif

#if !defined(SSE)
	free(PHIx);
	free(PHIy);
	free(PHIxuw);
	free(PHIyuw);
	free(S);
#else
	_mm_free(PHIx);
	_mm_free(PHIy);
	_mm_free(PHIxuw);
	_mm_free(PHIyuw);
	_mm_free(S);
#endif
}

#if defined(SSE) && defined(GASASM)
  /*#warning "blurredSignFunction: GAS+SSE"*/
  /*
  //------------------------------
  //--- Faster due to ASM and SSE
  //--- Programmed using GAS
  //--- Only works with IA32 or compatible...does not work in Intel 64
  //------------------------------
  */
  void blurredSignFunction(float *S, float *PHI, float *PHIx, float *PHIy, unsigned int rows, unsigned int cols, unsigned int frames)
  {
    
	  /*	Peng, Merriman, Osher, Zhao and Kang, "A PDE-Based Fast Local Level Set Method"
		S( \Phi_0 ) = \dfrac{\Phi_0}{ \sqrt(\Phi_0^2 + (\Delta x)^2 }			*/
    
	  unsigned int elems=rows*cols-1;
	  char mathStatus[512] __attribute__((aligned(16)));
	  
	  __asm (
			  "fxsave %[mathStatus]\n\t"	/*Store status of the FPU, MMX and SSE*/
			  "finit\n\t"
			  "pxor %%xmm7, %%xmm7\n\t"	/*Sets XMM7 to zero...used later on for comparisons*/
			  "movl %[elems], %%ecx\n\t"	/*Store number of elems to ECX*/
			  
			  /*Store pointers to registers*/
			  "movl %[PHI], 	%%eax\n\t"
			  "movl %[PHIx],	%%esi\n\t"
			  "movl %[PHIy],	%%edi\n\t"
			  "movl %[S],	%%edx\n\t"
	  
			  "loop_me:\n\t"
		  
				  /*S[i] = PHI[i]/sqrt( PHI[i]*OM[i] + sqrt(PHIx[i]*PHIx[i]+PHIy[i]*PHIy[i]) );*/
				  "movss (%%eax,%%ecx,4), %%xmm0\n\t"	/*OM into xmm0*/

				  "movss %%xmm0,		%%xmm1\n\t"		/*xmm0 into xmm1*/
				  "movss (%%esi,%%ecx,4),	%%xmm2\n\t"		/*PHIx into xmm2*/
				  "movss %%xmm2,		%%xmm3\n\t"		/*xmm2 into xmm3*/
				  "movss (%%edi,%%ecx,4),	%%xmm4\n\t"		/*PHIy into xmm4*/
				  "movss %%xmm4,		%%xmm5\n\t"		/*xmm4 into xmm5*/

				  "mulss %%xmm1, %%xmm0\n\t"	/*PHI^2*/
				  "mulss %%xmm3, %%xmm2\n\t"	/*PHIx^2*/
				  "mulss %%xmm5, %%xmm4\n\t"	/*PHIy^2*/

				  "addss %%xmm4, %%xmm2\n\t"	/*PHIx^2+PHIy^2*/
				  "sqrtss %%xmm2,  %%xmm2\n\t"	/*sqrt(PHIx^2+PHIy^2)*/

				  "addss %%xmm2, %%xmm0\n\t"	/*PHI^2+sqrt(PHIx^2+PHIy^2)*/
				  "comiss %%xmm0, %%Xmm7\n\t"	/*Compare if xmm0 is equivalent to xmm7 (containing 0)*/
				  "je zero\n\t"
				  "jmp notzero\n\t"

				  "zero:\n\t"
					  /*"pxor %%xmm0, %%xmm0\n\t"*/		/*Set xmm0 to 0*/
					  "movss (%%eax),	%%xmm0\n\t"	/*PHI into xmm0*/
					  "jmp store\n\t"
		  
				  "notzero:\n\t"

					  "rsqrtss %%xmm0,%%xmm0\n\t"	/*1/sqrt( PHI^2+sqrt(PHIx^2+PHIy^2) )*/
					  "mulss %%xmm1, %%xmm0\n\t"	/*PHI/sqrt( PHI^2+sqrt(PHIx^2+PHIy^2) )*/

				  "store:\n\t"
					  "movss %%xmm0, (%%edx,%%ecx,4)\n\t"
				  
				  "subl $1, %%ecx\n\t"
				  
				  "jae loop_me\n\t"

		  "fxrstor %[mathStatus]\n\t"	/*Return status of FPU, MMX and SSE*/

		  
		  :	[mathStatus]	"=m"(mathStatus)
		  :	[PHI]		"m"(PHI),
			[PHIx]		"m"(PHIx),
			[PHIy]		"m"(PHIy),
			[S]		"m"(S),
			[elems]		"m"(elems)
		  :"%eax","%ecx","%esi","%edi","%edx"
		  );
  }
#elif defined(SSE) && !defined(GASASM)
  /*#warning "blurredSignFunction: SSE"*/
  void blurredSignFunction(float *S, float *PHI, float *PHIx, float *PHIy, unsigned int rows, unsigned int cols, unsigned int frames)
  {
	  /*	Peng, Merriman, Osher, Zhao and Kang, "A PDE-Based Fast Local Level Set Method"
		S( \Phi_0 ) = \dfrac{\Phi_0}{ \sqrt(\Phi_0^2 + (\Delta x)^2 }			*/
	  
	  unsigned int i, elems = rows*cols*frames;
	  unsigned int sseloops = elems/4;
	  unsigned int remainloops = elems-4*sseloops;
	  unsigned int offset = sseloops*4;
	  
	  __m128 m1, m2, m3, m4_flt_eps;
  
	  __m128 *Sm = (__m128*) S;
	  __m128 *PHIm = (__m128*) PHI;
	  __m128 *PHIxm = (__m128*) PHIx;
	  __m128 *PHIym = (__m128*) PHIy;

	  m4_flt_eps = _mm_set_ps1( FLT_EPSILON );

	/*Process 4 data at a time*/
	for(i=0;i<sseloops;i++)
	{
		  /*S[i] = PHI[i]/sqrt( PHI[i]*PHI[i] + sqrt(PHIx[i]*PHIx[i]+PHIy[i]*PHIy[i]) + FLT_EPSILON );*/
		  /*S[i] = PHI[i]*InvSqrt( PHI[i]*PHI[i] + sqrt(PHIx[i]*PHIx[i]+PHIy[i]*PHIy[i]) + FLT_EPSILON );*/
		  m1 = _mm_mul_ps(*PHIm, *PHIm);
		  m2 = _mm_mul_ps(*PHIxm, *PHIxm);
		  m3 = _mm_mul_ps(*PHIym, *PHIym);
	  
		  m2 = _mm_add_ps(m2, m3);		/*PHIx^2 + PHIy^2*/
		  m2 = _mm_add_ps(m2, m4_flt_eps);	/*PHIx^2 + PHIy^2 + FLT_EPSILON*/
		  
		  m2 = _mm_sqrt_ps(m2);			/*sqrt( PHIx^2 + PHIy^2 + FLT_EPSILON ) */
		  m1 = _mm_add_ps(m2, m1);		/*PHI^2 + sqrt( PHIx^2 + PHIy^2 + FLT_EPSILON ) */
		  
		  m1 = _mm_rsqrt_ps(m1);		/*1/sqrt( PHI^2 + sqrt( PHIx^2 + PHIy^2 + FLT_EPSILON ) ) */
		  *Sm = _mm_mul_ps(*PHIm, m1);		/*PHI * 1/sqrt( PHI^2 + sqrt( PHIx^2 + PHIy^2 + FLT_EPSILON ) ) */

		  Sm++;
		  PHIm++;
		  PHIxm++;
		  PHIym++;
		}
	  /*Process the remaining data one at a time*/
	  for(i=0;i<remainloops;i++)
	  {
		  /*S[i] = PHI[i]/sqrt( PHI[i]*PHI[i] + sqrt(PHIx[i]*PHIx[i]+PHIy[i]*PHIy[i]) + FLT_EPSILON );*/
		  /*S[i] = PHI[i]*InvSqrt( PHI[i]*PHI[i] + sqrt(PHIx[i]*PHIx[i]+PHIy[i]*PHIy[i]) + FLT_EPSILON );*/
		  m1 = _mm_load_ss(&PHI[offset+i]);
		  m2 = _mm_load_ss(&PHIx[offset+i]);
		  m2 = _mm_load_ss(&PHIy[offset+i]);
		  
		  m1 = _mm_mul_ss(m1, m1);
		  m2 = _mm_mul_ss(m2, m2);
		  m3 = _mm_mul_ps(m3, m3);
		  
		  m2 = _mm_add_ss(m2, m3);		/*PHIx^2 + PHIy^2*/
		  m2 = _mm_add_ss(m2, m4_flt_eps);	/*PHIx^2 + PHIy^2 + FLT_EPSILON*/
		  
		  m2 = _mm_sqrt_ss(m2);			/*sqrt( PHIx^2 + PHIy^2 + FLT_EPSILON ) */
		  m1 = _mm_add_ss(m2, m1);		/*PHI^2 + sqrt( PHIx^2 + PHIy^2 + FLT_EPSILON ) */
		  
		  m2 = _mm_load_ss(&PHI[offset+i]);
		  m1 = _mm_rsqrt_ss(m1);		/*1/sqrt( PHI^2 + sqrt( PHIx^2 + PHIy^2 + FLT_EPSILON ) ) */
		  m1 = _mm_mul_ss(m1, m2);		/*PHI * 1/sqrt( PHI^2 + sqrt( PHIx^2 + PHIy^2 + FLT_EPSILON ) ) */
		  _mm_store_ss(&S[offset+i], m1);
	  }
	  
  }
#else
  /*#warning "blurredSignFunction: C"*/
  /*
  //-------------------------------------
  //--- Calculate blurred sign function
  //--- !!!SLOW DUE TO 2 SQRT:S!!!
  //-------------------------------------
  */
  void blurredSignFunction(float *S, float *PHI, float *PHIx, float *PHIy, unsigned int rows, unsigned int cols, unsigned int frames)
  {
	  /*	Peng, Merriman, Osher, Zhao and Kang, "A PDE-Based Fast Local Level Set Method"
		S( \Phi_0 ) = \dfrac{\Phi_0}{ \sqrt(\Phi_0^2 + (\Delta x)^2 }			*/
	  
	  unsigned int i, elems=rows*cols*frames-1;
	  float temp_result[3];

	  #pragma omp parallel for
	  for(i=0;i<=elems;i++)
	  {
		  /*S[i] = PHI[i]/sqrt( PHI[i]*PHI[i] + sqrt(PHIx[i]*PHIx[i]+PHIy[i]*PHIy[i]) + FLT_EPSILON );*/
		  /*S[i] = PHI[i]*InvSqrt( PHI[i]*PHI[i] + sqrt(PHIx[i]*PHIx[i]+PHIy[i]*PHIy[i]) + FLT_EPSILON );*/
		  temp_result[0] = PHI[i]*PHI[i];
		  temp_result[1] = PHIx[i]*PHIx[i];
		  temp_result[2] = PHIy[i]*PHIy[i];

		  temp_result[0] += FLT_EPSILON;
		  temp_result[1] += temp_result[2];
		  temp_result[0] += sqrt(temp_result[1]);

		  S[i] = PHI[i]/sqrt( temp_result[0] );
	  }
  }
#endif

/*
//-------------------------------------------------------------------------------------
//--- Godunov's method for discretising the hyperbolic $S(\Phi_0) | \Delta x |$ term
//-------------------------------------------------------------------------------------
*/
void godunovUpwind(	float *PHI2dx_out,
			float *PHI2dy_out,
			float *a_in,
			float *PHI_in,
			unsigned int rows,
			unsigned int cols,
			unsigned int frames)
{
omp_set_dynamic(1);
#pragma omp parallel
{
	/*Godunov's method using formulation by Rouy and Tourin
	if S>0
		PHIx^2 = max(max(PHIx_bd,0)^2,min(PHIx_fd,0).^2);
	if S<0
		PHIx^2 = max(min(PHIx_bd,0)^2,max(PHIx_fd,0).^2);*/

	int i,j,k,pos;
	float PHIx_bd; /*backward difference, horizontal*/
	float PHIx_fd; /*forward difference, horizontal*/
	float PHIy_bd; /*backward difference, vertical*/
	float PHIy_fd; /*forward difference, vertical*/
  
	#pragma omp for
	for(k=0;k<frames;k++)
	{
	for(j=0;j<cols;j++)
	{
		for(i=0;i<rows;i++)
		{
			pos = k*rows*cols+j*rows+i;

			/*Horizontal derivatives*/
			if(j==0)
			{
				PHIx_fd = PHI_in[pos+rows] - PHI_in[pos] ;
				PHIx_bd = 0.0f;
			}
			else if(j==cols-1)
			{
				PHIx_fd = 0.0f;
				PHIx_bd = PHI_in[pos] - PHI_in[pos-rows];
			}
			else
			{
				PHIx_fd = PHI_in[pos+rows] - PHI_in[pos] ;
				PHIx_bd = PHI_in[pos] - PHI_in[pos-rows];
			}
			/*Vertical derivatives*/
			if(i==0)
			{
				PHIy_fd = PHI_in[pos+1] - PHI_in[pos];
				PHIy_bd = 0.0f;
			}
			else if(i==rows-1)
			{
				PHIy_fd = 0.0f;
				PHIy_bd = PHI_in[pos] - PHI_in[pos-1];
			}
			else
			{
				PHIy_fd = PHI_in[pos+1] - PHI_in[pos];
				PHIy_bd = PHI_in[pos] - PHI_in[pos-1];
			}
			/*Rouy and Tourin formulation*/
			if(a_in[pos]>0.0f){
				PHI2dx_out[pos] = max( maxP2(PHIx_bd), minP2(PHIx_fd) );
				PHI2dy_out[pos] = max( maxP2(PHIy_bd), minP2(PHIy_fd) );
			}
			else{
				PHI2dx_out[pos] = max( minP2(PHIx_bd), maxP2(PHIx_fd) );
				PHI2dy_out[pos] = max( minP2(PHIy_bd), maxP2(PHIy_fd) );
			}
		}
	}
	}
}
}

/*
//-----------------------------------------
//--- Fast 'unprecise' inverse square root
//-----------------------------------------
*/
float InvSqrt (float x){
	/*John Carmack's method*/
	/*See http://www.beyond3d.com/content/articles/8/*/
	float xhalf = 0.5f*x;
	int i = *(int*)&x;
	
	i = 0x5f3759df - (i>>1);
	x = *(float*)&i;
	x = x*(1.5f - xhalf*x*x);

	return x;
}

/* THIS SHOULD BE OBSOLETE!!!!
//-------------------------------------
//--- Calculate blurred sign function
//-------------------------------------
*/
/*void blurredSignFunctionSSE(float *S, float *OM, float *PHIx, float *PHIy, unsigned int rows, unsigned int cols)
{
	unsigned int i, elems=rows*cols-1;
	char mathStatus[512] __attribute__((aligned(16)));
	float eps = FLT_EPSILON;

	/*Store status of the FPU, MMX and SSE*/
/*	__asm__ ( 	
			"fxsave %0\n\t"
			"finit\n\t"
			"pxor %%xmm7, %%xmm7\n\t"	/*Sets XMM7 to zero...used later on for comparisons*/
/*			:"=m"(mathStatus)
		);
	for(i=0;i<=elems;i++)
	{
		__asm (		
				/*S[i] = OM[i]/sqrt( OM[i]*OM[i] + sqrt(PHIx[i]*PHIx[i]+PHIy[i]*PHIy[i]) );*/
				/*Calculate pointers to OM[i]*/
/*				"movl %[OM], 	%%eax\n\t"
				"movl %[PHIx],	%%ecx\n\t"
				"movl %[PHIy],	%%edx\n\t"
				
				"leal (%%eax,%[i],4),	%%eax\n\t"
				"leal (%%ecx,%[i],4),	%%ecx\n\t"
				"leal (%%edx,%[i],4),	%%edx\n\t"
				
				"movss (%%eax),	%%xmm0\n\t"	/*OM into xmm0*/
/*				"movss %%xmm0,	%%xmm1\n\t"	/*xmm0 into xmm1*/
/*				"movss (%%ecx),	%%xmm2\n\t"	/*PHIx into xmm2*/
/*				"movss %%xmm2,	%%xmm3\n\t"	/*xmm2 into xmm3*/
/*				"movss (%%edx),	%%xmm4\n\t"	/*PHIy into xmm4*/
/*				"movss %%xmm4,	%%xmm5\n\t"	/*xmm4 into xmm5*/

/*				"mulss %%xmm1, %%xmm0\n\t"	/*OM^2*/
/*				"mulss %%xmm3, %%xmm2\n\t"	/*PHIx^2*/
/*				"mulss %%xmm5, %%xmm4\n\t"	/*PHIy^2*/

/*				"addss %%xmm4, %%xmm2\n\t"	/*PHIx^2+PHIy^2*/
/*				"sqrtss %%xmm2,  %%xmm2\n\t"	/*sqrt(PHIx^2+PHIy^2)*/

/*				"addss %%xmm2, %%xmm0\n\t"	/*OM^2+sqrt(PHIx^2+PHIy^2)*/
/*				"comiss %%xmm0, %%Xmm7\n\t"	/*Compare if xmm0 is equivalent to xmm7 (containing 0)*/
/*				"je zero\n\t"
				"jmp notzero\n\t"

				"zero:\n\t"
					/*"pxor %%xmm0, %%xmm0\n\t"*/	/*Set xmm0 to 0*/
/*					"movss (%%eax),	%%xmm0\n\t"	/*OM into xmm0*/
/*					"jmp store\n\t"
		
				"notzero:\n\t"
					"rsqrtss %%xmm0,%%xmm0\n\t"	/*1/sqrt( OM^2+sqrt(PHIx^2+PHIy^2) )*/
/*					"mulss %%xmm1, %%xmm0\n\t"	/*OM/sqrt( OM^2+sqrt(PHIx^2+PHIy^2) )*/

/*				"store:\n\t"
					"movl %[S], 	%%eax\n\t"
					"movss %%xmm0, (%%eax,%[i],4)\n\t"
				
			:
			:	[OM]		"m"(OM),
				[PHIx]		"m"(PHIx),
				[PHIy]		"m"(PHIy),
				[S]		"m"(S),
				[i]		"r"(i),	
				[EPS]		"m"(eps)
			:"%eax","%ecx","%edx"
			);
	}

	/*Return status of FPU, MMX and SSE*/
/*	__asm__ ( 	
			"fxrstor %0\n\t"
			:"=m"(mathStatus)
		);

}
*/
