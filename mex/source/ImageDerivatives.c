/*
Functions related to image derivative calculation

Author: Jarno Ralli
E-mail: jarno@ralli.fi

If you use this code, please reference (some) of my papers available at http://atc.ugr.es/~jarnor

Copyright 2011, Jarno Ralli

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

#include "imageDerivatives_c.h"

#if !defined( isnan )
  #define isnan(x)((x)!=(x))
#endif

void NOP(float *Result, float *It0, float *It1, float *operator, unsigned int rows, unsigned int cols)
{
	memcpy( Result, It1, rows*cols*sizeof(float) );
}
/*
//-----------------------------------------------------------
//--- Temporal convolution without angle wrap-around handling
//-----------------------------------------------------------
*/
void TemporalConvWO2(float *Result, float *It0, float *It1, float *operator, unsigned int rows, unsigned int cols)
{
	unsigned int i,j;
	unsigned int pos;
	unsigned int frame_offset = rows*cols;

	for(j=0;j<cols;j++)
	{
		pos = j*rows;
		
		for(i=0;i<rows;i++)
		{
			Result[pos] = It0[pos]*operator[0] + It1[pos]*operator[1];
			pos++;
		}
	}
}
/*
//-----------------------------------------------------------
//--- Vertical convolution without angle wrap-around handling
//-----------------------------------------------------------
*/
void VerticalConvWO5(float *Result, float *It0, float *It1, float *operator, unsigned int rows, unsigned int cols)
{
	unsigned int i,j;
	unsigned int pos,cntr=0;
	unsigned int frame_offset = rows*cols;

	float temp_result[5];

	for(j=0;j<cols;j++)
	{
		pos = j*rows;

		temp_result[0] = It1[pos]*operator[0];
		temp_result[1] = It1[pos]*operator[1];
		temp_result[2] = It1[pos]*operator[2];
		temp_result[3] = It1[pos+1]*operator[3];
		temp_result[4] = It1[pos+2]*operator[4];
		Result[pos] = temp_result[0]+temp_result[1]+temp_result[2]+temp_result[3]+temp_result[4];

		temp_result[0] = It1[pos]*operator[0];
		temp_result[1] = It1[pos]*operator[1];
		temp_result[2] = It1[pos+1]*operator[2];
		temp_result[3] = It1[pos+2]*operator[3];
		temp_result[4] = It1[pos+3]*operator[4];
		Result[pos+1] = temp_result[0]+temp_result[1]+temp_result[2]+temp_result[3]+temp_result[4];
		
		for(i=0;i<rows-4;i++)
		{
			temp_result[0] = It1[pos]*operator[0];
			temp_result[1] = It1[pos+1]*operator[1];
			temp_result[2] = It1[pos+2]*operator[2];
			temp_result[3] = It1[pos+3]*operator[3];
			temp_result[4] = It1[pos+4]*operator[4];
			
			Result[pos+2] = temp_result[0]+temp_result[1]+temp_result[2]+temp_result[3]+temp_result[4];
			pos++;
		}

		temp_result[0] = It1[pos]*operator[0];
		temp_result[1] = It1[pos+1]*operator[1];
		temp_result[2] = It1[pos+2]*operator[2];
		temp_result[3] = It1[pos+3]*operator[3];
		temp_result[4] = It1[pos+3]*operator[4];
		Result[pos+2] = temp_result[0]+temp_result[1]+temp_result[2]+temp_result[3]+temp_result[4];
		pos++;

		temp_result[0] = It1[pos]*operator[0];
		temp_result[1] = It1[pos+1]*operator[1];
		temp_result[2] = It1[pos+2]*operator[2];
		temp_result[3] = It1[pos+2]*operator[3];
		temp_result[4] = It1[pos+2]*operator[4];
		Result[pos+2] = temp_result[0]+temp_result[1]+temp_result[2]+temp_result[3]+temp_result[4];
	}
	
}
/*
//-----------------------------------------------------------
//--- Horizontal convolution without angle wrap-around handling
//-----------------------------------------------------------
*/
void HorizontalConvWO5(float *Result, float *It0, float *It1, float *operator, unsigned int rows, unsigned int cols)
{
	unsigned int i,j;
	unsigned int pos, ppos;
	unsigned int frame_offset = rows*cols;

	float temp_result[5];

	for(i=0;i<rows;i++)
	{
		pos = i;

		ppos = pos;
		temp_result[0] = It1[ppos];
		temp_result[1] = It1[ppos];
		temp_result[2] = It1[ppos];
		temp_result[3] = It1[ppos+=rows];
		temp_result[4] = It1[ppos+=rows];
		temp_result[0] = temp_result[0] * operator[0];
		temp_result[1] = temp_result[1] * operator[1];
		temp_result[2] = temp_result[2] * operator[2];
		temp_result[3] = temp_result[3] * operator[3];
		temp_result[4] = temp_result[4] * operator[4];
		Result[pos] = temp_result[0]+temp_result[1]+temp_result[2]+temp_result[3]+temp_result[4];

		ppos = pos;
		temp_result[0] = It1[ppos];
		temp_result[1] = It1[ppos];
		temp_result[2] = It1[ppos+=rows];
		temp_result[3] = It1[ppos+=rows];
		temp_result[4] = It1[ppos+=rows];
		temp_result[0] = temp_result[0] * operator[0];
		temp_result[1] = temp_result[1] * operator[1];
		temp_result[2] = temp_result[2] * operator[2];
		temp_result[3] = temp_result[3] * operator[3];
		temp_result[4] = temp_result[4] * operator[4];
		Result[pos+rows] = temp_result[0]+temp_result[1]+temp_result[2]+temp_result[3]+temp_result[4];

		for(j=0;j<cols-4;j++)
		{	
			ppos = pos;
			temp_result[0] = It1[ppos];
			temp_result[1] = It1[ppos+=rows];
			temp_result[2] = It1[ppos+=rows];
			temp_result[3] = It1[ppos+=rows];
			temp_result[4] = It1[ppos+=rows];

			temp_result[0] = temp_result[0] * operator[0];
			temp_result[1] = temp_result[1] * operator[1];
			temp_result[2] = temp_result[2] * operator[2];
			temp_result[3] = temp_result[3] * operator[3];
			temp_result[4] = temp_result[4] * operator[4];

			Result[pos+2*rows] = temp_result[0]+temp_result[1]+temp_result[2]+temp_result[3]+temp_result[4];
			pos+=rows;
		}

		ppos = pos;
		temp_result[0] = It1[ppos];
		temp_result[1] = It1[ppos+=rows];
		temp_result[2] = It1[ppos+=rows];
		temp_result[3] = It1[ppos+=rows];
		temp_result[4] = It1[ppos];
		temp_result[0] = temp_result[0] * operator[0];
		temp_result[1] = temp_result[1] * operator[1];
		temp_result[2] = temp_result[2] * operator[2];
		temp_result[3] = temp_result[3] * operator[3];
		temp_result[4] = temp_result[4] * operator[4];
		Result[pos+2*rows] = temp_result[0]+temp_result[1]+temp_result[2]+temp_result[3]+temp_result[4];
		pos+=rows;

		ppos = pos;
		temp_result[0] = It1[ppos];
		temp_result[1] = It1[ppos+=rows];
		temp_result[2] = It1[ppos+=rows];
		temp_result[3] = It1[ppos];
		temp_result[4] = It1[ppos];
		temp_result[0] = temp_result[0] * operator[0];
		temp_result[1] = temp_result[1] * operator[1];
		temp_result[2] = temp_result[2] * operator[2];
		temp_result[3] = temp_result[3] * operator[3];
		temp_result[4] = temp_result[4] * operator[4];
		Result[pos+2*rows] = temp_result[0]+temp_result[1]+temp_result[2]+temp_result[3]+temp_result[4];
	}
	
}

/*
//-----------------------------------------------------------
//--- Vertical convolution without angle wrap-around handling
//-----------------------------------------------------------
*/
void VerticalConvWO3(float *Result, float *It0, float *It1, float *operator, unsigned int rows, unsigned int cols)
{
	unsigned int i,j;
	unsigned int pos,cntr=0;
	unsigned int frame_offset = rows*cols;

	float temp_result[3];

	for(j=0;j<cols;j++)
	{
		pos = j*rows;

		temp_result[0] = It1[pos]*operator[0];
		temp_result[1] = It1[pos]*operator[1];
		temp_result[2] = It1[pos+1]*operator[2];
		Result[pos] = temp_result[0]+temp_result[1]+temp_result[2];

		for(i=0;i<rows-2;i++)
		{
			temp_result[0] = It1[pos]*operator[0];
			temp_result[1] = It1[pos+1]*operator[1];
			temp_result[2] = It1[pos+2]*operator[2];
			
			Result[pos+1] = temp_result[0]+temp_result[1]+temp_result[2];
			pos++;
		}

		temp_result[0] = It1[pos]*operator[0];
		temp_result[1] = It1[pos+1]*operator[1];
		temp_result[2] = It1[pos+1]*operator[2];
			
		Result[pos+1] = temp_result[0]+temp_result[1]+temp_result[2];
	}
}
/*
//-----------------------------------------------------------
//--- Horizontal convolution without angle wrap-around handling
//-----------------------------------------------------------
*/
void HorizontalConvWO3(float *Result, float *It0, float *It1, float *operator, unsigned int rows, unsigned int cols)
{
	unsigned int i,j;
	unsigned int pos, ppos;
	unsigned int frame_offset = rows*cols;

	float temp_result[3];

	for(i=0;i<rows;i++)
	{
		pos = i;

		ppos = pos;
		temp_result[0] = It1[ppos];
		temp_result[1] = It1[ppos];
		temp_result[2] = It1[ppos+=rows];
		temp_result[0] = temp_result[0] * operator[0];
		temp_result[1] = temp_result[1] * operator[1];
		temp_result[2] = temp_result[2] * operator[2];
		Result[pos] = temp_result[0]+temp_result[1]+temp_result[2];

		for(j=0;j<cols-2;j++)
		{	
			ppos = pos;
			temp_result[0] = It1[ppos];
			temp_result[1] = It1[ppos+=rows];
			temp_result[2] = It1[ppos+=rows];

			temp_result[0] = temp_result[0] * operator[0];
			temp_result[1] = temp_result[1] * operator[1];
			temp_result[2] = temp_result[2] * operator[2];

			Result[pos+rows] = temp_result[0]+temp_result[1]+temp_result[2];
			pos+=rows;
		}

		ppos = pos;
		temp_result[0] = It1[ppos];
		temp_result[1] = It1[ppos+=rows];
		temp_result[2] = It1[ppos];
		temp_result[0] = temp_result[0] * operator[0];
		temp_result[1] = temp_result[1] * operator[1];
		temp_result[2] = temp_result[2] * operator[2];
		Result[pos+rows] = temp_result[0]+temp_result[1]+temp_result[2];
	}
}

/*
//------------------------------------------------------------------
//--- First order image derivatives using Simoncelli-like operators
//------------------------------------------------------------------
*/
void fstSimoncelli_c(	struct matrixM *Idt,
			struct matrixM *Idx,
			struct matrixM *Idy,
			struct matrixM *It0,
			struct matrixM *It1,
			struct matrixM *temp,
			unsigned int size,	/*Size of the convolution operator, 3 or 5*/
			float *Smoother,
			float *SpatialDerivator,
			float *TemporalDerivator
		)
{
		unsigned int rows = It0->dimElems[0];
		unsigned int cols = It0->dimElems[1];
		unsigned int frames = 1;
		unsigned int k, offset;

		/*Function pointers to "operations"*/
		void (*IdtPreProcFunc)	(float *Result, float *It0, float *It1, float *operator, unsigned int rows, unsigned int cols);
		void (*IdtProcFunc)	(float *Result, float *It0, float *It1, float *operator, unsigned int rows, unsigned int cols);
		
		void (*IdxPreProcFunc)	(float *Result, float *It0, float *It1, float *operator, unsigned int rows, unsigned int cols);
		void (*IdxProcFunc)	(float *Result, float *It0, float *It1, float *operator, unsigned int rows, unsigned int cols);

		void (*IdyPreProcFunc)	(float *Result, float *It0, float *It1, float *operator, unsigned int rows, unsigned int cols);
		void (*IdyProcFunc)	(float *Result, float *It0, float *It1, float *operator, unsigned int rows, unsigned int cols);
		
		/*Bind the functions*/
		if(size==3)
		{
			/*Preprocessing*/
			IdtPreProcFunc = &NOP;
			IdxPreProcFunc = &VerticalConvWO3;
			IdyPreProcFunc = &HorizontalConvWO3;
				
			/*Actual processing*/
			IdtProcFunc = &TemporalConvWO2;
			IdxProcFunc = &HorizontalConvWO3;
			IdyProcFunc = &VerticalConvWO3;
		}
		else if(size==5)
		{
			/*Preprocessing*/
			IdtPreProcFunc = &NOP;
			IdxPreProcFunc = &VerticalConvWO5;
			IdyPreProcFunc = &HorizontalConvWO5;
			/*Actual processing*/
			IdtProcFunc = &TemporalConvWO2;
			IdxProcFunc = &HorizontalConvWO5;
			IdyProcFunc = &VerticalConvWO5;
		}
		else
		{
			mexErrMsgTxt("Only convolution sizes 3 and 5 are implemented!");
		}
		
		if( It0->ndims>2 )
			frames = It0->dimElems[2];

		/*Handle all the frames*/
		for(k=0;k<frames;k++)
		{
			offset = k*rows*cols;
	
			/*Temporal derivative Idt (no preprocessing)*/
			IdtProcFunc	( &Idt->data[offset], &It0->data[offset], &It1->data[offset], TemporalDerivator, rows, cols );
	
			/*Horizontal spatial derivative Idx*/
			IdxPreProcFunc	( &temp->data[offset], &It0->data[offset], &It1->data[offset], Smoother, rows, cols );
			IdxProcFunc	( &Idx->data[offset], &It0->data[offset], &temp->data[offset], SpatialDerivator, rows, cols );
	
			/*Vertical spatial derivative Idy*/
			IdyPreProcFunc	( &temp->data[offset], &It0->data[offset], &It1->data[offset], Smoother, rows, cols );
			IdyProcFunc	( &Idy->data[offset], &It0->data[offset], &temp->data[offset], SpatialDerivator, rows, cols );

		}
}
/*
//-------------------------------------------------------------------
//--- Second order image derivatives using Simoncelli-like operators
//-------------------------------------------------------------------
*/
void sndSimoncelli_c(	struct matrixM *Idxt,
			struct matrixM *Idyt,
			struct matrixM *Idxx,
			struct matrixM *Idyy,
			struct matrixM *Idxy,
			struct matrixM *It0,
			struct matrixM *It1,
			struct matrixM *temp1,
			struct matrixM *temp2,
			unsigned int size,	/*Size of the convolution operator => 5*/
			float *Smoother,
			float *fstSpatialDerivator,
			float *sndSpatialDerivator,
			float *TemporalDerivator
		)
{
		unsigned int rows = It0->dimElems[0];
		unsigned int cols = It0->dimElems[1];
		unsigned int frames = 1;
		unsigned int k, offset;
		
		/*Function pointers to "operations"*/
		void (*IdxtPreProc1Func)(float *Result, float *It0, float *It1, float *operator, unsigned int rows, unsigned int cols);
		void (*IdxtPreProc2Func)(float *Result, float *It0, float *It1, float *operator, unsigned int rows, unsigned int cols);
		void (*IdxtProcFunc)	(float *Result, float *It0, float *It1, float *operator, unsigned int rows, unsigned int cols);

		void (*IdytPreProc1Func)(float *Result, float *It0, float *It1, float *operator, unsigned int rows, unsigned int cols);
		void (*IdytPreProc2Func)(float *Result, float *It0, float *It1, float *operator, unsigned int rows, unsigned int cols);
		void (*IdytProcFunc)	(float *Result, float *It0, float *It1, float *operator, unsigned int rows, unsigned int cols);
		
		void (*IdxxPreProcFunc)	(float *Result, float *It0, float *It1, float *operator, unsigned int rows, unsigned int cols);
		void (*IdxxProcFunc)	(float *Result, float *It0, float *It1, float *operator, unsigned int rows, unsigned int cols);

		void (*IdyyPreProcFunc)	(float *Result, float *It0, float *It1, float *operator, unsigned int rows, unsigned int cols);
		void (*IdyyProcFunc)	(float *Result, float *It0, float *It1, float *operator, unsigned int rows, unsigned int cols);

		void (*IdxyPreProcFunc)	(float *Result, float *It0, float *It1, float *operator, unsigned int rows, unsigned int cols);
		void (*IdxyProcFunc)	(float *Result, float *It0, float *It1, float *operator, unsigned int rows, unsigned int cols);
		
		/*Bind the functions*/
		/*Preprocessing*/
		IdxtPreProc1Func = &VerticalConvWO5;
		IdxtPreProc2Func = &HorizontalConvWO5;
		IdytPreProc1Func = &HorizontalConvWO5;
		IdytPreProc2Func = &VerticalConvWO5;

		IdxxPreProcFunc = &VerticalConvWO5;
		IdyyPreProcFunc = &HorizontalConvWO5;
		IdxyPreProcFunc = &HorizontalConvWO5;

		/*Actual processing*/
		IdxtProcFunc = &TemporalConvWO2;
		IdytProcFunc = &TemporalConvWO2;

		IdxxProcFunc = &HorizontalConvWO5;
		IdyyProcFunc = &VerticalConvWO5;
		IdxyProcFunc = &VerticalConvWO5;
		
		if( It0->ndims>2 )
			frames = It0->dimElems[2];
		
		for(k=0;k<frames;k++)
		{
			offset = k*rows*cols;

			/*Temporal derivatives*/
			IdxtPreProc1Func( &temp2->data[offset], &It0->data[offset], &It0->data[offset], Smoother, rows, cols );
			IdxtPreProc2Func( &temp1->data[offset], &temp2->data[offset], &temp2->data[offset], fstSpatialDerivator, rows, cols );
			IdxtPreProc1Func( &Idxt->data[offset], &It1->data[offset], &It1->data[offset], Smoother, rows, cols );
			IdxtPreProc2Func( &temp2->data[offset], &Idxt->data[offset], &Idxt->data[offset], fstSpatialDerivator, rows, cols );
			IdxtProcFunc	( &Idxt->data[offset], &temp1->data[offset], &temp2->data[offset], TemporalDerivator, rows, cols );
	
			IdytPreProc1Func( &temp2->data[offset], &It0->data[offset], &It0->data[offset], Smoother, rows, cols );
			IdytPreProc2Func( &temp1->data[offset], &temp2->data[offset], &temp2->data[offset], fstSpatialDerivator, rows, cols );
			IdytPreProc1Func( &Idyt->data[offset], &It1->data[offset], &It1->data[offset], Smoother, rows, cols );
			IdytPreProc2Func( &temp2->data[offset], &Idyt->data[offset], &Idyt->data[offset], fstSpatialDerivator, rows, cols );
			IdytProcFunc	( &Idyt->data[offset], &temp1->data[offset], &temp2->data[offset], TemporalDerivator, rows, cols );
	
			/*Horizontal spatial derivative Idxx*/
			IdxxPreProcFunc	( &temp1->data[offset], &It0->data[offset], &It1->data[offset], Smoother, rows, cols );
			IdxxProcFunc	( &Idxx->data[offset], &It0->data[offset], &temp1->data[offset], sndSpatialDerivator, rows, cols );
	
			/*Vertical spatial derivative Idyy*/
			IdyyPreProcFunc	( &temp1->data[offset], &It0->data[offset], &It1->data[offset], Smoother, rows, cols );
			IdyyProcFunc	( &Idyy->data[offset], &It0->data[offset], &temp1->data[offset], sndSpatialDerivator, rows, cols );
	
			/*Cross spatial derivative Idxy*/
			IdxyPreProcFunc	( &temp1->data[offset], &It0->data[offset], &It1->data[offset], fstSpatialDerivator, rows, cols );
			IdxyProcFunc	( &Idxy->data[offset], &It0->data[offset], &temp1->data[offset], fstSpatialDerivator, rows, cols);

		}
}
