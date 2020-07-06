/*
gcc -fPIC -c urs_ext.c -I/usr/include/python2.5 -o urs_ext.o -Wall -O
gcc -shared urs_ext.o  -o urs_ext.so
*/

#include "structure.h"
#include "math.h"
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <float.h>
#include <time.h>

#define MAX_FILE_NAME_LENGTH 128
#define NODATA 99.0
#define EPSILON  0.00001

#define DEBUG 0

#define POFFSET 5 //Number of site_params

static int *fros=NULL;  // First recorded output step 
static int *lros=NULL;  // Last recorded output step 
static struct tgsrwg* mytgs0=NULL;

static long numDataMax=0;


/*The MUX file format 


*/



/////////////////////////////////////////////////////////////////////////
//Auxiliary functions
void fillDataArray(int ista, int total_number_of_stations, int nt, int ig, int *nst, 
                   int *nft, float *data, int *istart_p, 
                   int *istop_p, float *muxData)
{
    int it, last_it, jsta;
    long int offset=0;


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


long getNumData(const int *fros, const int *lros, const int total_number_of_stations)
/* calculates the number of data in the data block of a mux file */
/* based on the first and last recorded output steps for each gauge */ 
{
    int ista, last_output_step;
    long numData = 0;

    last_output_step = 0;   
    for(ista = 0; ista < total_number_of_stations; ista++)
        if(*(fros + ista) != -1)
        {
            numData += *(lros + ista) - *(fros + ista) + 1;
            last_output_step = (last_output_step < *(lros+ista) ? 
                *(lros+ista):last_output_step);
        }   
        numData += last_output_step*total_number_of_stations; /* these are the t records */
        return numData;
}

/////////////////////////////////////////////////////////////////////////
//Internal Functions
int _read_mux2_headers(int numSrc, 
                       char **muxFileNameArray, 
                       int* total_number_of_stations,
                       int* number_of_time_steps,
                       double* delta_t,
                       //long* numDataMax,
                       int verbose)
{
    FILE *fp;
    int numsta, i, j;
    struct tgsrwg *mytgs=0;
    char *muxFileName;                                                                  
    char susMuxFileName;
    long numData;
    size_t elements_read; // fread return value
    int block_size;

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

    // Loop over all sources, read headers and check compatibility
    for (i = 0; i < numSrc; i++)
    {
        muxFileName = muxFileNameArray[i];

        // Open the mux file
        if((fp = fopen(muxFileName, "rb")) == NULL)
        {
            char *err_msg = strerror(errno);

            fprintf(stderr, "cannot open file '%s': %s\n", muxFileName, err_msg);
            return -1;  
        }

        if (!i)
        {
            elements_read = fread(total_number_of_stations, sizeof(int), 1, fp);
            if ((int) elements_read == 0 && ferror(fp)){
                fprintf(stderr, "Error reading total number of stations\n");
                fclose(fp);
                return -2;
            }

            fros = (int*) malloc(*total_number_of_stations*numSrc*sizeof(int));
            lros = (int*) malloc(*total_number_of_stations*numSrc*sizeof(int));

            mytgs0 = (struct tgsrwg*) malloc(*total_number_of_stations*sizeof(struct tgsrwg));
            mytgs = (struct tgsrwg*) malloc(*total_number_of_stations*sizeof(struct tgsrwg));

            block_size = *total_number_of_stations*sizeof(struct tgsrwg);
            elements_read = fread(mytgs0, block_size , 1, fp);
            if ((int) elements_read == 0 && ferror(fp)){
                fprintf(stderr, "Error reading mytgs0\n");
                fclose(fp);
                return -2;
            }
        }
        else
        {
            // Check that the mux files are compatible
            elements_read = fread(&numsta, sizeof(int), 1, fp);
            if ((int) elements_read == 0 && ferror(fp)){
                fprintf(stderr, "Error reading numsta\n");
                fclose(fp);
                return -2;
            }

            if(numsta != *total_number_of_stations)
            {
                fprintf(stderr,"%s has different number of stations to %s\n", 
                    muxFileName, 
                    muxFileNameArray[0]);
                fclose(fp);
                return -1;   
            }

            block_size = numsta*sizeof(struct tgsrwg);
            elements_read = fread(mytgs, block_size, 1, fp); 
            if ((int) elements_read == 0 && ferror(fp)){
                fprintf(stderr, "Error reading mgtgs\n");
                fclose(fp);
                return -2;
            }	    


            for (j = 0; j < numsta; j++)
            {
                if (mytgs[j].dt != mytgs0[j].dt)
                {
                    fprintf(stderr, "%s has different sampling rate to %s\n", 
                        muxFileName, 
                        muxFileNameArray[0]);
                    fclose(fp);
                    return -1;            
                }   
                if (mytgs[j].nt != mytgs0[j].nt)
                {
                    fprintf(stderr, "%s has different series length to %s\n", 
                        muxFileName, 
                        muxFileNameArray[0]);
                    fclose(fp);
                    return -1;            
                }

                if (mytgs[j].nt != mytgs0[0].nt)
                {
                    printf("Station 0 has different series length to Station %d\n", j); 
                }
            }
        }

        /* Read the start and stop times for this source */
        elements_read = fread(fros + i*(*total_number_of_stations), 
            *total_number_of_stations*sizeof(int), 1, fp);
        if ((int) elements_read == 0 && ferror(fp)){
            fprintf(stderr, "Error reading start times\n");
            fclose(fp);
            return -3;
        }	    


        elements_read = fread(lros + i*(*total_number_of_stations), 
            *total_number_of_stations*sizeof(int), 1, fp);
        if ((int) elements_read == 0 && ferror(fp)){
            fprintf(stderr, "Error reading stop times\n");
            fclose(fp);
            return -3;
        }	    	      

        /* Compute the size of the data block for this source */
        numData = getNumData(fros + i*(*total_number_of_stations), 
            lros + i*(*total_number_of_stations), 
            (*total_number_of_stations));

        /* Sanity check */
        if (numData < 0)
        {
            fprintf(stderr,"Size of data block appears to be negative!\n");
            return -1;        
        }

        if (numDataMax < numData)
        {
            numDataMax = numData;
        }

        fclose(fp);          
    }


    // Store time resolution and number of timesteps    
    // These are the same for all stations as tested above, so 
    // we take the first one.
    *delta_t = (double)mytgs0[0].dt;
    *number_of_time_steps = mytgs0[0].nt;

    free(mytgs);

    return 0; // Succesful execution
}


float** _read_mux2(int numSrc, 
                   char **muxFileNameArray, 
                   float *weights, 
                   double *params, 
                   int *number_of_stations,
                   long *permutation,
                   int verbose)
{
    FILE *fp;
    int total_number_of_stations, i, isrc, ista, k;
    char *muxFileName;
    int istart=-1, istop=-1;
    int number_of_selected_stations;
    float *muxData=NULL; // Suppress warning
    long numData;
    long *perm = NULL;
    long *permutation_temp = NULL;

    int len_sts_data, error_code;
    float **sts_data;
    float *temp_sts_data;

    long int offset;

    int number_of_time_steps, N;
    double delta_t;

    size_t elements_read;

    // Shorthands pointing to memory blocks for each source
    int *fros_per_source=NULL;     
    int *lros_per_source=NULL;         


    error_code = _read_mux2_headers(numSrc, 
        muxFileNameArray, 
        &total_number_of_stations,
        &number_of_time_steps,
        &delta_t,
        verbose);
    if (error_code != 0) {
        printf("urs.c: Internal function _read_mux2_headers failed: Error code = %d\n", 
            error_code);

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
        permutation_temp = (long *) malloc(number_of_selected_stations*sizeof(long));
        if (permutation_temp == NULL)
        {
            printf("ERROR: Memory for permutation_temp could not be allocated.\n");
            return NULL;
        }
        
        for (i = 0; i < number_of_selected_stations; i++)
        {
            permutation_temp[i] = (long) i;  
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
    sts_data = (float**) malloc(number_of_selected_stations*sizeof(float*));
    if (sts_data == NULL)
    {
        printf("ERROR: Memory for sts_data could not be allocated.\n");
        return NULL;
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
            return NULL;
        }
    }

    temp_sts_data = (float*) calloc(len_sts_data, sizeof(float));
    if (temp_sts_data == NULL)
    {
        printf("ERROR: Memory for temp_sts_data could not be allocated.\n");
        return NULL;
    }

    muxData = (float*) calloc(numDataMax, sizeof(float));
    if (temp_sts_data == NULL)
    {
        printf("ERROR: Memory for muxData could not be allocated.\n");
        return NULL;
    }

    // Loop over all sources
    for (isrc = 0; isrc < numSrc; isrc++)
    {

        // Shorthands to local memory
        fros_per_source = (int*) fros + isrc*total_number_of_stations; 
        lros_per_source = (int*) lros + isrc*total_number_of_stations; 	    


        // Read in data block from mux2 file
        muxFileName = muxFileNameArray[isrc];
        if((fp = fopen(muxFileName, "rb")) == NULL)
        {
            fprintf(stderr, "cannot open file %s\n", muxFileName);
            free(muxData);
            free(temp_sts_data);
            free(muxData);

            return NULL;                    
        }

        if (verbose){
            printf("Reading mux file %s\n", muxFileName);
        }

        offset = (long int)sizeof(int) + total_number_of_stations*(sizeof(struct tgsrwg) + 2*sizeof(int));
        //printf("\n offset %i ", (long int)offset);
        fseek(fp, offset, 0);

        numData = getNumData(fros_per_source, 
            lros_per_source, 
            total_number_of_stations);
        // Note numData is larger than what it has to be.		     
        //elements_read = fread(muxData, ((int) numData)*sizeof(float), 1, fp); 
        elements_read = fread(muxData, (size_t) sizeof(float), (size_t) numData, fp); 
        //printf("\n elements_read  %d, ", (int)elements_read);
        //printf("\n ferror(fp)  %d, ", (int)ferror(fp));
        if ((int) elements_read == 0 && ferror(fp)) {
            fprintf(stderr, "Error reading mux data: %s", strerror(errno));
            if (errno == EFAULT)        // error 14 in /usr/include/asm-generic/errno-base.h
            {
                fprintf(stderr, "NOTE: This error has been seen before in low memory systems with no swap.\n");
            }

            fclose(fp);
            free(muxData);
            free(temp_sts_data);
            free(muxData);

            return NULL;
        }	

        fclose(fp);  

        // loop over stations present in the permutation array 
        //     use ista with mux data
        //     use i with the processed data to be returned         
        for(i = 0; i < number_of_selected_stations; i++)
        {               

            ista = (int) perm[i]; // Get global index into mux data  

            // fill the data0 array from the mux file, and weight it
            fillDataArray(ista, 
                total_number_of_stations, 
                number_of_time_steps,
                mytgs0[ista].ig, // Grid number (if -1 fill with zeros)
                fros_per_source, 
                lros_per_source, 
                temp_sts_data, 
                &istart, 
                &istop, 
                muxData);

            // Weight appropriately and add
            for(k = 0; k < mytgs0[ista].nt; k++)
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
                sts_data[i][N] = (float)mytgs0[ista].geolat;
                sts_data[i][N+1] = (float)mytgs0[ista].geolon;
                sts_data[i][N+2] = (float)mytgs0[ista].z;
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
    free(fros);
    free(lros);
    free(mytgs0);

    if (permutation_temp)
    {
        free(permutation_temp);
    }

    return sts_data;
}
