/*
gcc -fPIC -c urs_ext.c -I/usr/include/python2.5 -o urs_ext.o -Wall -O
gcc -shared urs_ext.o  -o urs_ext.so
*/

#include "structure.h"
#include "math.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <float.h>
#include <time.h>
#include <stdint.h>
#include "anuga_typedefs.h"

#define MAX_FILE_NAME_LENGTH 128
#define NODATA 99.0
#define EPSILON  0.00001

#define DEBUG 0

#define POFFSET 5 //Number of site_params

/* Header state shared between _read_mux2_headers and _read_mux2. This
   used to be four file-scope statics, which made concurrent read_mux2 calls
   corrupt each other; now _read_mux2 owns one of these on its stack. */
struct mux2_headers {
    int32_t *fros;           // First recorded output step
    int32_t *lros;           // Last recorded output step
    struct tgsrwg *mytgs0;
    anuga_int numDataMax;
};

static void free_mux2_headers(struct mux2_headers *hd)
{
    free(hd->fros);   hd->fros = NULL;
    free(hd->lros);   hd->lros = NULL;
    free(hd->mytgs0); hd->mytgs0 = NULL;
}


/*The MUX file format 


*/



/////////////////////////////////////////////////////////////////////////
//Auxiliary functions
void fillDataArray(int32_t ista, int32_t total_number_of_stations, int32_t nt, int32_t ig, int32_t *nst, 
                   int32_t *nft, float *data, int32_t *istart_p, 
                   int32_t *istop_p, float *muxData)
{
    int32_t it, last_it, jsta;
    anuga_int offset=0;


    last_it = -1;
    /* Make arrays of starting and finishing time steps for the tide gauges */
    /* and fill them from the file */

    /* Update start and stop timesteps for this gauge */
    if (nst[ista]!= -1)
    {
        if(*istart_p == -1)
        {
            *istart_p = nst[ista];
        }
        else
        {
            *istart_p = ((nst[ista] < *istart_p) ? nst[ista] : *istart_p);
        }
    }

    if (nft[ista] != -1)
    {
        if (*istop_p == -1)
        {
            *istop_p = nft[ista];
        }
        else
        {
            *istop_p = ((nft[ista] < *istop_p) ? nft[ista] : *istop_p);
        }
    }     

    if (ig == -1 || nst[ista] == -1) /* currently ig==-1 => nst[ista]==-1 */
    {
        /* gauge never started recording, or was outside of all grids, 
        fill array with 0 */
        for(it = 0; it < nt; it++)
        {
            data[it] = 0.0;
        }
    }   
    else
    {
        for(it = 0; it < nt; it++)
        {
            last_it = it;
            /* skip t record of data block */
            offset++;
            /* skip records from earlier tide gauges */
            for(jsta = 0; jsta < ista; jsta++)
                if(it + 1 >= nst[jsta] && it + 1 <= nft[jsta])
                    offset++;

            /* deal with the tide gauge at hand */
            if(it + 1 >= nst[ista] && it + 1 <= nft[ista])
            {
                /* gauge is recording at this time */
                memcpy(data + it, muxData + offset, sizeof(float));

                //printf("%d: muxdata=%f\n", it, muxData[offset]); 
                //printf("data[%d]=%f, offset=%d\n", it, data[it], offset); 
                offset++;
            }
            else if (it + 1 < nst[ista])
            {
                /* gauge has not yet started recording */
                data[it] = 0.0;
            }   
            else
                /* gauge has finished recording */
            {
                data[it] = NODATA;
                break;
            }

            /* skip records from later tide gauges */
            for(jsta = ista + 1; jsta < total_number_of_stations; jsta++)
                if(it + 1 >= nst[jsta] && it+1 <= nft[jsta])
                    offset++;
        }

        if(last_it < nt - 1)
            /* the loop was exited early because the gauge had 
            finished recording */
            for(it = last_it+1; it < nt; it++)
                data[it] = NODATA;
    }
} 


char isdata(float x)
{
    if(x < NODATA + EPSILON && NODATA < x + EPSILON)
    {
        return 0;
    }
    else
    {
        return 1;  
    }
}


anuga_int getNumData(const int32_t *fros, const int32_t *lros, const int32_t total_number_of_stations)
/* calculates the number of data in the data block of a mux file */
/* based on the first and last recorded output steps for each gauge */ 
{
    int32_t ista, last_output_step;
    anuga_int numData = 0;

    last_output_step = 0;   
    for(ista = 0; ista < total_number_of_stations; ista++){

        if(*(fros + ista) != -1)
        {
            numData += *(lros + ista) - *(fros + ista) + 1;
            last_output_step = (last_output_step < *(lros+ista) ? 
                *(lros+ista):last_output_step);
        }   
    }
    numData += last_output_step*total_number_of_stations; /* these are the t records */
    return numData;
}

/////////////////////////////////////////////////////////////////////////
//Internal Functions
int32_t _read_mux2_headers(int32_t numSrc, 
                       char **muxFileNameArray, 
                       int32_t* total_number_of_stations,
                       int32_t* number_of_time_steps,
                       double* delta_t,
                       struct mux2_headers *hd,
                       int32_t verbose)
{
    FILE *fp = NULL;
    int32_t numsta, i, j;
    struct tgsrwg *mytgs=0;
    char *muxFileName;                                                                  
    char susMuxFileName;
    anuga_int numData;
    size_t elements_read; // fread return value
    int32_t rc = 0;      // error code handed to the fail path
    int32_t block_size;

    /* Allocate space for the names and the weights and pointers to the data*/

    /* Check that the input files have mux2 extension*/
    susMuxFileName = 0;
    for(i = 0; i < numSrc; i++)
    { 
        muxFileName = muxFileNameArray[i];
        if(!susMuxFileName && strcmp(muxFileName + strlen(muxFileName) - 4, 
            "mux2") != 0)
        {
            susMuxFileName = 1;
            break;
        }
    }

    if(susMuxFileName)
    {
        printf("\n**************************************************************************\n");
        printf("   WARNING: This program operates only on multiplexed files in mux2 format\n"); 
        printf("   At least one input file name does not end with mux2\n");
        printf("   Check your results carefully!\n");
        printf("**************************************************************************\n\n");
    }   

    if (verbose)
    {
        printf("Reading mux header information\n");
    }

    if (numSrc <= 0)
    {
        fprintf(stderr, "No mux2 source files given\n");
        return -1;
    }
    if (hd->fros != NULL || hd->lros != NULL || hd->mytgs0 != NULL)
    {
        fprintf(stderr, "_read_mux2_headers: the header struct must start empty\n");
        return -4;
    }

    // Loop over all sources, read headers and check compatibility
    for (i = 0; i < numSrc; i++)
    {
        muxFileName = muxFileNameArray[i];

        // Open the mux file
        if((fp = fopen(muxFileName, "rb")) == NULL)
        {
            char *err_msg = strerror(errno);

            fprintf(stderr, "cannot open file '%s': %s\n", muxFileName, err_msg);
            rc = -1; goto fail;
        }

        // First file: read the station count and allocate the header arrays.
        // Keyed on the arrays not existing yet rather than on i == 0, which is
        // the same thing on every real run but lets gcc's analyser see that
        // the allocation cannot happen twice.
        if (hd->mytgs0 == NULL)
        {
            elements_read = fread(total_number_of_stations, sizeof(int32_t), 1, fp);
            if ((int32_t) elements_read == 0 && ferror(fp)){
                fprintf(stderr, "Error reading total number of stations\n");
                rc = -2; goto fail;
            }

            hd->fros = (int32_t*) malloc(*total_number_of_stations*numSrc*sizeof(int32_t));
            hd->lros = (int32_t*) malloc(*total_number_of_stations*numSrc*sizeof(int32_t));

            hd->mytgs0 = (struct tgsrwg*) malloc(*total_number_of_stations*sizeof(struct tgsrwg));
            mytgs = (struct tgsrwg*) malloc(*total_number_of_stations*sizeof(struct tgsrwg));
            if (!hd->fros || !hd->lros || !hd->mytgs0 || !mytgs)
            {
                fprintf(stderr, "Could not allocate header arrays for %d stations\n",
                        (int)*total_number_of_stations);
                rc = -3; goto fail;
            }

            block_size = *total_number_of_stations*sizeof(struct tgsrwg);
            elements_read = fread(hd->mytgs0, block_size , 1, fp);
            if ((int32_t) elements_read == 0 && ferror(fp)){
                fprintf(stderr, "Error reading hd->mytgs0\n");
                rc = -2; goto fail;
            }
        }
        else
        {
            // Check that the mux files are compatible
            elements_read = fread(&numsta, sizeof(int32_t), 1, fp);
            if ((int32_t) elements_read == 0 && ferror(fp)){
                fprintf(stderr, "Error reading numsta\n");
                rc = -2; goto fail;
            }

            if(numsta != *total_number_of_stations)
            {
                fprintf(stderr,"%s has different number of stations to %s\n", 
                    muxFileName, 
                    muxFileNameArray[0]);
                rc = -1; goto fail;   
            }

            if (mytgs == NULL)   // cannot happen after the first file; keeps the analyser honest
            {
                rc = -4; goto fail;
            }
            block_size = numsta*sizeof(struct tgsrwg);
            elements_read = fread(mytgs, block_size, 1, fp); 
            if ((int32_t) elements_read == 0 && ferror(fp)){
                fprintf(stderr, "Error reading mgtgs\n");
                rc = -2; goto fail;
            }	    


            for (j = 0; j < numsta; j++)
            {
                if (mytgs[j].dt != hd->mytgs0[j].dt)
                {
                    fprintf(stderr, "%s has different sampling rate to %s\n", 
                        muxFileName, 
                        muxFileNameArray[0]);
                    rc = -1; goto fail;            
                }   
                if (mytgs[j].nt != hd->mytgs0[j].nt)
                {
                    fprintf(stderr, "%s has different series length to %s\n", 
                        muxFileName, 
                        muxFileNameArray[0]);
                    rc = -1; goto fail;            
                }

                if (mytgs[j].nt != hd->mytgs0[0].nt)
                {
                    printf("Station 0 has different series length to Station %d\n", j); 
                }
            }
        }

        /* Read the start and stop times for this source */
        elements_read = fread(hd->fros + i*(*total_number_of_stations), 
            *total_number_of_stations*sizeof(int32_t), 1, fp);
        if ((int32_t) elements_read == 0 && ferror(fp)){
            fprintf(stderr, "Error reading start times\n");
            rc = -3; goto fail;
        }	    


        elements_read = fread(hd->lros + i*(*total_number_of_stations), 
            *total_number_of_stations*sizeof(int32_t), 1, fp);
        if ((int32_t) elements_read == 0 && ferror(fp)){
            fprintf(stderr, "Error reading stop times\n");
            rc = -3; goto fail;
        }	    	      

        /* Compute the size of the data block for this source */
        numData = getNumData(hd->fros + i*(*total_number_of_stations), 
            hd->lros + i*(*total_number_of_stations), 
            (*total_number_of_stations));

        /* Sanity check */
        if (numData < 0)
        {
            fprintf(stderr,"Size of data block appears to be negative!\n");
            rc = -1; goto fail;
        }

        if (hd->numDataMax < numData)
        {
            hd->numDataMax = numData;
        }

        fclose(fp);
        fp = NULL;
    }


    // Store time resolution and number of timesteps    
    // These are the same for all stations as tested above, so 
    // we take the first one.
    *delta_t = (double)hd->mytgs0[0].dt;
    *number_of_time_steps = hd->mytgs0[0].nt;

    free(mytgs);

    return 0; // Succesful execution

fail:
    if (fp) fclose(fp);
    free(mytgs);
    free_mux2_headers(hd);
    return rc;
}


float** _read_mux2(int32_t numSrc, 
                   char **muxFileNameArray, 
                   float *weights, 
                   double *params, 
                   int32_t *number_of_stations,
                   anuga_int *permutation,
                   int32_t verbose)
{
    FILE *fp;
    int32_t total_number_of_stations, i, isrc, ista, k;
    char *muxFileName;
    int32_t istart=-1, istop=-1;
    int32_t number_of_selected_stations;
    float *muxData=NULL; // Suppress warning
    anuga_int numData;
    anuga_int *perm = NULL;
    anuga_int *permutation_temp = NULL;

    int32_t len_sts_data, error_code;
    float **sts_data = NULL;
    float *temp_sts_data = NULL;

    anuga_int offset;

    int32_t number_of_time_steps, N;
    double delta_t;

    size_t elements_read;

    // Shorthands pointing to memory blocks for each source
    int32_t *fros_per_source=NULL;     
    int32_t *lros_per_source=NULL;         


    struct mux2_headers headers = {NULL, NULL, NULL, 0};
    struct mux2_headers *hd = &headers;

    error_code = _read_mux2_headers(numSrc, 
        muxFileNameArray, 
        &total_number_of_stations,
        &number_of_time_steps,
        &delta_t,
        hd,
        verbose);
    if (error_code != 0) {
        printf("urs.c: Internal function _read_mux2_headers failed: Error code = %d\n", 
            error_code);
        free_mux2_headers(hd);   // already released on that path; harmless and explicit
        return NULL;
    }


    // Apply rule that an empty permutation file means 'take all stations'
    // We could change this later by passing in None instead of the empty 
    // permutation.
    number_of_selected_stations = *number_of_stations;  
    if (number_of_selected_stations == 0)
    {
        number_of_selected_stations = total_number_of_stations;  

        // Return possibly updated number of stations
        *number_of_stations = total_number_of_stations;     

        // Create the Identity permutation vector
        permutation_temp = (anuga_int *) malloc(number_of_selected_stations*sizeof(anuga_int));
        if (permutation_temp == NULL)
        {
            printf("ERROR: Memory for permutation_temp could not be allocated.\n");
            goto fail;
        }
        
        for (i = 0; i < number_of_selected_stations; i++)
        {
            permutation_temp[i] = (anuga_int) i;  
        }

        perm = permutation_temp;
    }
    else
    {
        perm = permutation;
    }

    // The params array is used only for passing data back to Python.
    params[0] = (double) number_of_selected_stations;
    params[1] = (double) delta_t;
    params[2] = (double) number_of_time_steps;

    // Make array(s) to hold demuxed data for stations given in the 
    // permutation file 
    sts_data = (float**) calloc(number_of_selected_stations, sizeof(float*));
    if (sts_data == NULL)
    {
        printf("ERROR: Memory for sts_data could not be allocated.\n");
        goto fail;
    }

    // For each selected station, allocate space for its data
    len_sts_data = number_of_time_steps + POFFSET; // Max length of each timeseries?
    for (i = 0; i < number_of_selected_stations; i++)
    {
        // Initialise sts_data to zero
        sts_data[i] = (float*) calloc(len_sts_data, sizeof(float));
        if (sts_data[i] == NULL)
        {
            printf("ERROR: Memory for sts_data could not be allocated.\n");
            goto fail;
        }
    }

    temp_sts_data = (float*) calloc(len_sts_data, sizeof(float));
    if (temp_sts_data == NULL)
    {
        printf("ERROR: Memory for temp_sts_data could not be allocated.\n");
        goto fail;
    }

    muxData = (float*) calloc(hd->numDataMax, sizeof(float));
    if (temp_sts_data == NULL)
    {
        printf("ERROR: Memory for muxData could not be allocated.\n");
        goto fail;
    }

    // Loop over all sources
    for (isrc = 0; isrc < numSrc; isrc++)
    {

        // Shorthands to local memory
        fros_per_source = (int32_t*) hd->fros + isrc*total_number_of_stations; 
        lros_per_source = (int32_t*) hd->lros + isrc*total_number_of_stations; 	    


        // Read in data block from mux2 file
        muxFileName = muxFileNameArray[isrc];
        if((fp = fopen(muxFileName, "rb")) == NULL)
        {
            fprintf(stderr, "cannot open file %s\n", muxFileName);
            goto fail;                    
        }

        if (verbose){
            printf("Reading mux file %s\n", muxFileName);
        }

        offset = (anuga_int) sizeof(int32_t) + total_number_of_stations*(sizeof(struct tgsrwg) + 2*sizeof(int32_t));
        //printf("\n offset %i ", (anuga_int int32_t)offset);
        fseek(fp, offset, 0);

        numData = getNumData(fros_per_source, 
            lros_per_source, 
            total_number_of_stations);
        // Note numData is larger than what it has to be.		     
        //elements_read = fread(muxData, ((int32_t) numData)*sizeof(float), 1, fp); 
        elements_read = fread(muxData, (size_t) sizeof(float), (size_t) numData, fp); 
        //printf("\n elements_read  %d, ", (int32_t)elements_read);
        //printf("\n ferror(fp)  %d, ", (int32_t)ferror(fp));
        if ((int32_t) elements_read == 0 && ferror(fp)) {
            fprintf(stderr, "Error reading mux data: %s", strerror(errno));
            if (errno == EFAULT)        // error 14 in /usr/include/asm-generic/errno-base.h
            {
                fprintf(stderr, "NOTE: This error has been seen before in low memory systems with no swap.\n");
            }

            fclose(fp);
            goto fail;
        }	

        fclose(fp);  

        // loop over stations present in the permutation array 
        //     use ista with mux data
        //     use i with the processed data to be returned         
        for(i = 0; i < number_of_selected_stations; i++)
        {               

            ista = (int32_t) perm[i]; // Get global index into mux data  

            // fill the data0 array from the mux file, and weight it
            fillDataArray(ista, 
                total_number_of_stations, 
                number_of_time_steps,
                hd->mytgs0[ista].ig, // Grid number (if -1 fill with zeros)
                fros_per_source, 
                lros_per_source, 
                temp_sts_data, 
                &istart, 
                &istop, 
                muxData);

            // Weight appropriately and add
            for(k = 0; k < hd->mytgs0[ista].nt; k++)
            {
                if((isdata(sts_data[i][k])) && isdata(temp_sts_data[k]))
                {
                    sts_data[i][k] += temp_sts_data[k] * weights[isrc];
                }
                else
                {
                    sts_data[i][k] = NODATA;
                }
                //printf("%d: temp_sts_data[%d]=%f\n", i, k, temp_sts_data[k]);	    

            }


            // Update metadata (e.g. start time and end time)
            N = number_of_time_steps;

            if (isrc == 0) {
                // Assign values for first source
                sts_data[i][N] = (float)hd->mytgs0[ista].geolat;
                sts_data[i][N+1] = (float)hd->mytgs0[ista].geolon;
                sts_data[i][N+2] = (float)hd->mytgs0[ista].z;
                sts_data[i][N+3] = (float)fros_per_source[ista];
                sts_data[i][N+4] = (float)lros_per_source[ista];
            } else {
                // Update first and last timesteps for subsequent sources
                if (sts_data[i][N+3] > (float)fros_per_source[ista]) {		
                    if (verbose) {
                        printf("Adjusting start time for station %d and source %d",
                            ista, isrc);
                        printf(" from %f to %f\n", 
                            sts_data[i][N+3], 
                            (float) fros_per_source[ista]);  
                    }
                    sts_data[i][N+3] = (float) fros_per_source[ista];
                }

                if (sts_data[i][N+4] < (float) lros_per_source[ista]) {		
                    if (verbose) {
                        printf("Adjusting end time for station %d and source %d",
                            ista, isrc);
                        printf(" from %f to %f\n", 
                            sts_data[i][N+4], 
                            (float) lros_per_source[ista]);  
                    }
                    sts_data[i][N+4] = (float) lros_per_source[ista];
                }		
            }
        }
    }

    //printf("sts_data[1,8]=%f\n", sts_data[1][8]);

    free(muxData);
    free(temp_sts_data);
    free_mux2_headers(hd);
    free(permutation_temp);

    return sts_data;

fail:
    free(muxData);
    free(temp_sts_data);
    if (sts_data)
    {
        for (i = 0; i < number_of_selected_stations; i++) free(sts_data[i]);
        free(sts_data);
    }
    free(permutation_temp);
    free_mux2_headers(hd);
    return NULL;
}
