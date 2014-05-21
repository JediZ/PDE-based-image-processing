/*
Functions related to RANSAC

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

#include "ransac.h"

void RANSAC( struct matrixM *A_in,				/*Data to be fitted A_in * model_out = B_in*/
	     struct matrixM *B_in,				/*Target*/
	     struct matrixM *model_out,				/*The best model found by RANSAC*/
	     struct matrixM *model_in,				/*A model participating in the competition*/
	     struct matrixM *error_out,				/*Residual (or fitting error) of each data for the best model*/
	     /*function pointer as an argument*/
	     void (*createModel) (	struct matrixM *,	/* A */
					struct matrixM *,	/* B */
					struct matrixM *,	/* model_out */
					struct matrixM *,	/* model_in */
					struct matrixM *,	/* error_out */
					unsigned int *,		/* rand_set */
					unsigned int ),		/* n */
	     float err_thr, 					/*Error threshold for determining if the data fits the model*/
	     float min_set_size,	 			/*Minimum number of close data (percentage; 0.0-1.0) values required to assert that a model fits well to data*/
	     unsigned int iter, 				/*Number of iterations*/
	     unsigned int n )					/*Number of data used for generating the model*/
{
    /*
    Generic RANSAC
    createModel() function creates the actual model and calculates the corresponding error vector based on the model and the data: pointers
    to the function is given as parameter.
    
    Storage and pointers: model_out and model_storage are used to store the model(s) while current_model and best_model point to the
    storage objects so that best_model always points to the storage object containing the best model so far. Same goes for error_out
    error_storage and current_error and best_error
    */
    
    unsigned int *rand_set = NULL, i, j, nr_inliers, abs_min_set_size, data_length, best_inlr_nr_inliers;
    unsigned int found_licit_model;
    float current_error_sum, best_error_sum, err_thr2;
    
    struct matrixM *model_storage, *error_storage;						/*Storage for model and error*/
    struct matrixM *bestinlr_model_storage, *bestinlr_error_storage;				/*Storage for "best inlier" model*/
    struct matrixM *best_model, *best_error, *current_model, *current_error, *matrixMptr;	/*Pointers to storage objects*/
    
    /*Definitions etc.*/
    err_thr2 = err_thr*err_thr;
    data_length = A_in->dimElems[0];
    /*TODO: round()-function is not necessarily included in math.h...seems it was included in C99. Nicer solution should be used!
    abs_min_set_size = (unsigned int)round( min_set_size*(float)data_length );*/
    abs_min_set_size = (unsigned int)( min_set_size*(float)data_length + 0.5f );
    
    /*Initialize storages*/
    model_storage = (struct matrixM*)malloc( 1*sizeof(struct matrixM) );
    error_storage = (struct matrixM*)malloc( 1*sizeof(struct matrixM) );
    bestinlr_model_storage = (struct matrixM*)malloc( 1*sizeof(struct matrixM) );
    bestinlr_error_storage = (struct matrixM*)malloc( 1*sizeof(struct matrixM) );

    model_storage->ndims = model_out->ndims;			/*Taking a short cut...can cause problems!*/
    model_storage->dimElems = model_out->dimElems;		/*Taking a short cut...can cause problems!*/
    model_storage->data = (float*)malloc( model_out->dimElems[0]*model_out->dimElems[1]*sizeof(float) );
    
    error_storage->ndims = error_out->ndims;			/*Taking a short cut...can cause problems!*/
    error_storage->dimElems = error_out->dimElems;		/*Taking a short cut...can cause problems!*/
    error_storage->data = (float*)malloc( error_out->dimElems[0]*error_out->dimElems[1]*sizeof(float) );
    
    bestinlr_model_storage->ndims = model_out->ndims;		/*Taking a short cut...can cause problems!*/
    bestinlr_model_storage->dimElems = model_out->dimElems;	/*Taking a short cut...can cause problems!*/
    bestinlr_model_storage->data = (float*)malloc( model_out->dimElems[0]*model_out->dimElems[1]*sizeof(float) );
    
    bestinlr_error_storage->ndims = error_out->ndims;		/*Taking a short cut...can cause problems!*/
    bestinlr_error_storage->dimElems = error_out->dimElems;	/*Taking a short cut...can cause problems!*/
    bestinlr_error_storage->data = (float*)malloc( error_out->dimElems[0]*error_out->dimElems[1]*sizeof(float) );
    
    /*Initialize model and error pointers*/
    best_model = model_out;
    best_error = error_out;
    best_error_sum = FLT_MAX;
    best_inlr_nr_inliers = 0;
    found_licit_model = 0;

    current_model = model_storage;
    current_error = error_storage;
    
    /*Reserve space for the random set vector*/
    rand_set = (unsigned int*)malloc( n*sizeof(unsigned int) );
 
    /*---GIVEN MODEL---*/
    /*If a model is given (model_in) calculate error vector and sum of the error vector
    and use these as a starting point for the search*/
    if( model_in->data != NULL )
    {

      nr_inliers = 0;
      current_error_sum = 0.0f;
      
      createModel( A_in, B_in, best_model, model_in, best_error, rand_set, n );
      
      /*Go through all the members and find number of inliers and sum of the errors for the inliers*/
      for(j=0;j<data_length;j++)
      {
	  /*printf("%f\n",best_error->data[j]);*/
	  if(best_error->data[j]<=err_thr2)
	  {	
	    current_error_sum += best_error->data[j];
	    nr_inliers++;
	  }
      }
      
      /*Update best model error*/
      if( nr_inliers>=abs_min_set_size )
      {
	  /*Update sum of the best error*/
	  /*printf("updating error sum to:%f\n",current_error_sum);*/
	  best_error_sum = current_error_sum;
	  found_licit_model = 1;
      }else
      {
	  /*printf("updating error sum to:FLT_MAX\n");*/
	  best_error_sum = FLT_MAX;
      }

    }
    /*printf("best_error_sum:%f\n",best_error_sum);
    printf("---END INITIAL MODEL---\n");*/
    
    /*---RANSAC MAIN LOOP---*/
    i = 0;
    while(i<iter)
    {
	
	/*printf("\nmain loop\n");*/
	nr_inliers = 0;
	current_error_sum = 0.0f;
      
	/*Generate a random vector pointing to data members*/
	randVect( rand_set, 0, data_length-1, n );
	/*Generate a model based on the chosen data members and calculate the error vector*/
	createModel( A_in, B_in, current_model, NULL, current_error, rand_set, n );
	
	/*Go through all the members and find number of inliers and sum of the errors for the inliers*/
	for(j=0;j<data_length;j++)
	{
	  /*printf("%f ", current_error->data[j]);*/

	  if(current_error->data[j]<=err_thr2)
	  {	
		current_error_sum += current_error->data[j];
		nr_inliers++;
	  }
	}
	/*Update best model IF sufficiently many inliers AND error smaller than the best error so far*/
	if( (nr_inliers>=abs_min_set_size)&&(current_error_sum<best_error_sum) )
	{
	      found_licit_model = 1;
	  
	      /*Swap pointers so that the best_model points to the current_model*/
	      matrixMptr = best_model;
	      best_model = current_model;
	      current_model = matrixMptr;
	      /*Swap pointers so that the best_error points to the current_error*/
	      matrixMptr = best_error;
	      best_error = current_error;
	      current_error = matrixMptr;
	  
	      /*Update sum of the best error*/
	      best_error_sum = current_error_sum;
	}else if( (nr_inliers>=best_inlr_nr_inliers)&&(found_licit_model==0) )
	{
	      best_inlr_nr_inliers = nr_inliers;
	      memcpy( bestinlr_model_storage->data, current_model->data, model_out->dimElems[0]*model_out->dimElems[1]*sizeof(float) );
	      memcpy( bestinlr_error_storage->data, current_error->data, error_out->dimElems[0]*error_out->dimElems[1]*sizeof(float) );
	}

	i++;
	
    }

    if( found_licit_model == 1 )
    {
 
	memcpy( model_out->data, best_model->data, model_out->dimElems[0]*model_out->dimElems[1]*sizeof(float) );
	memcpy( error_out->data, best_error->data, error_out->dimElems[0]*error_out->dimElems[1]*sizeof(float) );
    }
    else
    {
	/*mexPrintf("hups...using bestlinr!\n");*/
	memcpy( model_out->data, bestinlr_model_storage->data, model_out->dimElems[0]*model_out->dimElems[1]*sizeof(float) );
	memcpy( error_out->data, bestinlr_error_storage->data, error_out->dimElems[0]*error_out->dimElems[1]*sizeof(float) );
    }

    /*Householding*/
    if( rand_set != NULL )		free( rand_set );
    if( model_storage->data != NULL )	free( model_storage->data );
    if( model_storage != NULL)		free( model_storage );
    if( error_storage->data != NULL)	free( error_storage->data );
    if( error_storage != NULL)		free( error_storage );

}

void randVect(unsigned int *vect, unsigned int min, unsigned int max, unsigned int length)
{
  static char init = 0;
  unsigned int i = 0;
  
  if( vect==NULL  )
    vect = (unsigned int*)malloc( length*sizeof(unsigned int) );
  
  if (init == 0)
  {
    srand(time(NULL));
    init = 1;
  }

  /*
   * Formula:  
   *    rand() % N   <- To get a number between 0 - N-1
   *    Then add the result to min, giving you 
   *    a random number between min - max.
   */  
  while(i<length)
  {	vect[i] = (rand() % (max - min + 1) + min);
	i++;
  }

}

