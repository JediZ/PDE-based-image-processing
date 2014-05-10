/*
Functions related to image interpolation

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

#include "imageInterpolation.h"

#ifdef _MSC_VER
	#define INFINITY (DBL_MAX+DBL_MAX)
	#define NAN (INFINITY-INFINITY)
#endif

/*
 //-------------------
 //--- INTERPOLATE ---
 //-------------------

 //Interpolates along the dimensions in the following order: stack, column and row.
 //This way the multiplicator values need to be calculated only once per X- Y-position.
 */
void bilinInterp2( struct matrixM *Iout,
		   struct matrixM *Iin,
		   struct matrixM *X,
		   struct matrixM *Y,
		   float NaN)
{
	unsigned int i, j, k;
	unsigned int nrows = Iin->dimElems[0];
	unsigned int ncols = Iin->dimElems[1];
	unsigned int nframes = 1;
	
	unsigned int pos;
	unsigned int iopos, ioposx1, ioposy1, ioposx1y1, inpos;
	unsigned int frame_offset, current_frame_offset;

	unsigned int x, y;		/*Pixel base position (x,y)*/
	float x_f, y_f;			/*Pixel 'fraction' position in unit square*/
	float w_00, w_10, w_01, w_11;
	
	#ifndef NAN
		float NAN = 0.0f/0.0f;		/*Okay...let's do it the ugly way...depending on the compiler, might cause problems!*/
	#endif

	frame_offset = nrows*ncols;
	
	if(Iout->ndims>2)
	    nframes = Iout->dimElems[2];

	/*--- Columns ---*/
	for(j=0;j<ncols;j++)
	{
		/*--- Rows ---*/
		for(i=0;i<nrows;i++)
		{
			/*--- Interpolation multiplicator values ---*/
			/* Current position index */
			pos = j*nrows + i;
			/* Base position of pixel in question */
			x = (unsigned int)floor( X->data[pos] - 1.0f );
			y = (unsigned int)floor( Y->data[pos] - 1.0f );

			/* Check legimity of position */
			if( (x>=0) && (x<ncols) && (y>=0) && (y<nrows) )
			{
			      /*Position inside unit square*/
			      x_f = X->data[pos] - 1.0f - (float)x;
			      y_f = Y->data[pos] - 1.0f - (float)y;

			      /*Interpolation multiplicators*/
			      w_00 = (1.0f-x_f) * (1.0f-y_f);
			      w_10 = x_f * (1.0f-y_f);
			      w_01 = (1.0f-x_f) * y_f;
			      w_11 = x_f * y_f;

			      /*--- Frames ---*/
			      for(k=0;k<nframes;k++) 			
			      {
				      current_frame_offset = k*frame_offset;

				     /*Old image position indices*/
				      iopos =		(current_frame_offset + nrows*x + y);
				      ioposx1 = iopos;
				      ioposy1 = iopos;
				      ioposx1y1 = iopos;
				      if(x<ncols-1)			ioposx1 += nrows;
				      if(y<nrows-1)			ioposy1 += 1;
				      if( (x<ncols-1) && (y<nrows-1) )	ioposx1y1 = ioposx1y1 + nrows + 1;
				      
				      /*These might point outside of the image, depending on the for-loops
				      ioposx1 =		(current_frame_offset + nrows*(x+1) + y );
				      ioposy1 =		(current_frame_offset + nrows*x + y + 1);
				      ioposx1y1 =	(current_frame_offset + nrows*(x+1) + y + 1);*/
				      
				      /*New image position index*/
				      inpos =		(current_frame_offset + pos);

				      Iout->data[inpos] = w_00 * Iin->data[iopos] 
							  + w_10 * Iin->data[ioposx1]
							  + w_01 * Iin->data[ioposy1]
							  + w_11 * Iin->data[ioposx1y1];
			      }
			}
			else
			{
			      /*--- Frames ---*/
			      for(k=0;k<nframes;k++) 			
			      {
				      current_frame_offset = k*frame_offset;
				      inpos =	(unsigned int)(current_frame_offset + pos);

				      Iout->data[inpos] = NaN;
			      }
			  
			}
 		}
 	}
 }
