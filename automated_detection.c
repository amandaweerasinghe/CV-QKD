
/*code for automated detection of header,phase correction and excess noise calculation
Author: Amanda Weerasinghe 
05/2022*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#define PI 3.14159265

//function declaration for phase correction
int phase_correction(float *bob_phase_sig, float *bob_phase_ref, float *bob_amp_sig, float *corrected_phase, float *corrx, float *corrp);

//function declaration for excess noise calculation 
float excess_noise (float *alice_x, float *corrx, float *shot_x);

// start of main function 
int main () {
    FILE * pFile, * pFile1;
    long lSize,lSize1;
    char * buffer, *buffer1;
    size_t result, result1;

    /*x quadrature file*/
    pFile = fopen ( "ADC0-11-5.bin" , "rb" );
    if (pFile==NULL) {fputs ("File error",stderr); exit (1);}

    // obtain file size:
    fseek (pFile , 0 , SEEK_END);
    lSize = ftell (pFile);
    rewind (pFile);
    printf("size:%ld\n",lSize);
    // allocate memory to contain the whole file:
    buffer = (char*) malloc (sizeof(char)*lSize);
    if (buffer == NULL) {fputs ("Memory error",stderr); exit (2);}

    // copy the file into the buffer:
    result = fread (buffer,1,lSize,pFile);
    if (result != lSize) {fputs ("Reading error",stderr); exit (3);}
 
    /* the whole file is now loaded in the memory buffer. */
    printf("buffer:%d\n",buffer[3]);
    // terminate
    fclose (pFile);

    /*p quadrature file*/
    pFile1 = fopen ( "ADC1-11-5.bin" , "rb" );
    if (pFile1==NULL) {fputs ("File error",stderr); exit (4);}
    // obtain file size:
    fseek (pFile1 , 0 , SEEK_END);
    lSize1 = ftell (pFile1);
    rewind (pFile1);
    printf("size:%ld\n",lSize1);
    // allocate memory to contain the whole file:
    buffer1 = (char*) malloc (sizeof(char)*lSize1);
    if (buffer1 == NULL) {fputs ("Memory error",stderr); exit (5);}

    // copy the file into the buffer:
    result1 = fread (buffer1,1,lSize1,pFile1);
    if (result1 != lSize1) {fputs ("Reading error",stderr); exit (6);}
 
    /* the whole file is now loaded in the memory buffer. */
    printf("buffer1:%d\n",buffer1[3]);
    // terminate
    fclose (pFile1);

    //conversion of a and p values into phase and amplitude values 
    int a=0; //counter for for loop for cart2pol conversion
    float *theta, *rho; 
    theta = (float*) malloc (sizeof(float)*lSize);
    if (theta == NULL) {fputs ("Memory error",stderr); exit (2);}
    rho = (float*) malloc (sizeof(float)*lSize);
    if (rho == NULL) {fputs ("Memory error",stderr); exit (2);}
    while(a<lSize){
        rho[a]=hypot(buffer1[a],buffer[a]);
        theta[a]=atan2f(buffer1[a],buffer[a]);
        if(theta[a]<0){theta[a]=theta[a]+2*PI;}
        if(a>0&&a<100){printf("rho[%d]:%f\n",a,rho[a]);}
        a=a+1;
    }

    /*reading phase and pulse files sent from alice*/
    char *token;
    int j=0; //counter for outer while loop that reads line by line
    float *phase;
    FILE *fp =fopen("phase_gaussian_header_10000copy.txt","r");
    if(fp==NULL){
        perror("Unable to open file!");
        exit(1);
    }
    //float phase[600000]={0};
     phase= (float*) malloc (sizeof(float)*600000);
     if (phase == NULL) {fputs ("Memory error",stderr); exit (7);}
    char chunk[128]; //string for reading each line
    while (fgets(chunk,sizeof(chunk),fp)!= NULL){
        token=strtok(chunk,"\n"); //token to separate columns 
        //while loop to separate columns in each line and save the value in second column to float array 
        while(token !=NULL){
            phase[j]=strtod(token,NULL);//convert string to float
            phase[j]=(phase[j]+1)*PI;
            if(j>199016 && j<199078){
                printf("phase[%d]= %f\n",j,phase[j]); //store in float array
                }
            token=strtok(NULL,","); //call for next token
            }
        j+=1;//increment outer loop counter
        }
    fclose(fp);

    char *token1;
    int k=0; //counter for outer while loop that reads line by line
    float *pulse;
    FILE *fp1 =fopen("pulse_gaussian_header_10000copy.txt","r");
    if(fp1==NULL){
        perror("Unable to open file!");
        exit(1);
    }
    pulse= (float*) malloc (sizeof(float)*600000);
    if (pulse == NULL) {fputs ("Memory error",stderr); exit (8);}
    char chunk1[128]; //string for reading each line
    while (fgets(chunk1,sizeof(chunk1),fp1)!= NULL){
        token1=strtok(chunk1,"\n"); //token to separate columns 
        //while loop to separate columns in each line and save the value in second column to float array 
        while(token1 !=NULL){
            pulse[k]=strtod(token1,NULL);//convert string to float
            pulse[k]=pulse[k]+1;
            if(k>500000 && k<500100){
                printf("pulse[%d]= %f\n",k,pulse[k]); //store in float array
                }
            token1=strtok(NULL,","); //call for next token
            }
        k+=1;//increment outer loop counter
        }
    fclose(fp1);
   
    /*choose correct pulse and phase values from the text files*/
    int first_alice=100012, spacing_alice=50, block_length=10000,counter_alice=0;
    float phase_alice[10000]={0},amplitude_alice[10000]={0}; //memory in the stack for phase and amplitude values of alice
    float alice_x[10000]={0},alice_p[10000]={0};
    //assigning the phase and amplitude values of alice to the heap one element at a time
    while(counter_alice<block_length){
        phase_alice[counter_alice]=phase[first_alice+spacing_alice*(counter_alice)]; //alice's phase
        amplitude_alice[counter_alice]=pulse[first_alice+spacing_alice*(counter_alice)]; //alice's amplitude

        //converting phase and amplitude to x and p quadrature values
        alice_x[counter_alice]=amplitude_alice[counter_alice]*cosf(phase_alice[counter_alice]);
        alice_p[counter_alice]=amplitude_alice[counter_alice]*sinf(phase_alice[counter_alice]);

        if(counter_alice<11){
        printf("phase_alice[%d]=%f\n",counter_alice,phase_alice[counter_alice]);
        printf("amplitude_alice[%d]:%f\n",counter_alice,amplitude_alice[counter_alice]);
        printf("alicex[%d]:%f\n",counter_alice,alice_x[counter_alice]);
        printf("alicep[%d]:%f\n",counter_alice,alice_p[counter_alice]);
        }
        counter_alice=counter_alice+1; //incrementing counter for the while loop
    }
    free(phase); //free heap for phase read from txt file
    free(pulse); //free heap for pulse read from txt file

    int ref_spacing=160, sig_spacing=160,first_ref_yes=0,counter_ref=0;

    //locating the first reference before the first header
    while(first_ref_yes==0){
        if(rho[counter_ref]>119 && counter_ref>100){first_ref_yes=1;}
        counter_ref=counter_ref+1;
    }
    printf("counter_ref:%d\n",counter_ref);

    int ref_index=counter_ref, header_yes=0;

    //locating the start of the header
    while(header_yes==0){
        if(rho[ref_index]>118){ref_index=ref_index+ref_spacing;}
        if(rho[ref_index]>98 && rho[ref_index]<116 && rho[ref_index+ref_spacing]>98 && rho[ref_index+ref_spacing]<116){
            header_yes=1;
        }
    }

    int header_index=ref_index, ref_yes=0;
    printf("counter_ref:%d\n",header_index);

    //locating the end of the header and the start of the data block 
    while(ref_yes==0){
        if(rho[header_index]>98 && rho[header_index]<116){header_index=header_index+ref_spacing;}
        if(rho[header_index]>119 && rho[header_index+ref_spacing]>119){ref_yes=1;}
    }
    printf("counter_ref:%d\n",header_index);
    int counter_first_ref=1, new_header_index;

    // .....................NEW COMMENT...............................
    //locating the start to middle of the first reference if the automatic detection has located on to the last portiion of the ref pulse 
    while(counter_first_ref<11){
        if(rho[header_index-counter_first_ref]>119){
            new_header_index=header_index-counter_first_ref;
            
        }
        counter_first_ref=counter_first_ref+1;
    }
   

    int start_sig=new_header_index-ref_spacing/2, start_ref=new_header_index; //starting points of ref and signal pulses
    int counter_data=0; //counter for loop for ref and sig pulse value storage 
    //array declarations for x and p values at Bob
    float bob_x_sig[10000]={0}, bob_p_sig[10000]={0}, bob_amp_sig[10000]={0},bob_phase_sig[10000]={0}; //signals
    float bob_x_ref[10000]={0}, bob_p_ref[10000]={0}, bob_amp_ref[10000]={0},bob_phase_ref[10000]={0}; //references
    float shot_x[10000]={0}, shot_p[10000]={0};
    float corrected_phase[10000]={0},corrx[10000]={0},corrp[10000]={0}; //arrays for phase correction
    float phase_error[10000]={0}; //array for storing phase error
    float phase_error_mean=0; //phase error mean 

    //separation of x and p values at bob and data values for shot noise calculation
    while(counter_data<block_length){
        //signal pulses
        bob_x_sig[counter_data]=buffer[start_sig+counter_data*sig_spacing];
        bob_p_sig[counter_data]=buffer1[start_sig+counter_data*sig_spacing];
        bob_amp_sig[counter_data]=rho[start_sig+counter_data*sig_spacing];
        bob_phase_sig[counter_data]=theta[start_sig+counter_data*sig_spacing];

        //reference pulses
        bob_x_ref[counter_data]=buffer[start_ref+counter_data*ref_spacing];
        bob_p_ref[counter_data]=buffer1[start_ref+counter_data*ref_spacing];
        bob_amp_ref[counter_data]=rho[start_ref+counter_data*ref_spacing];
        bob_phase_ref[counter_data]=theta[start_ref+counter_data*ref_spacing];

        //shot noise data
        shot_x[counter_data]=buffer[start_sig+40+counter_data*sig_spacing];
        shot_p[counter_data]=buffer1[start_sig+40+counter_data*sig_spacing];

        //correct phase using phase of adjacent ref pulses
        phase_correction(&bob_phase_sig[counter_data], &bob_phase_ref[counter_data],&bob_amp_sig[counter_data], &corrected_phase[counter_data],&corrx[counter_data],&corrp[counter_data]);

        //..........calculation of phase error:: NEW COMMENT.............
        phase_error[counter_data]=corrected_phase[counter_data]-phase_alice[counter_data];
        //wrap around to 2Pi
        if((phase_error[counter_data]>PI)|| (phase_error[counter_data]< (-1*PI))){
            phase_error[counter_data]=phase_error[counter_data]-2*PI;
        }
        //........phase error mean calculation:: NEW COMMENT..............
        phase_error_mean=phase_error_mean+phase_error[counter_data];

        if(counter_data<30){
            printf("corrx[%d]=%f\n", counter_data,corrx[counter_data]);
            printf("corrp[%d]=%f\n", counter_data,corrp[counter_data]);
            printf("phase_error[%d]=%f\n",counter_data,phase_error[counter_data]);
        }
        counter_data=counter_data+1;
    }
    //.....NEW COMMENT:phase_error_mean calculation......
    phase_error_mean=fabs(phase_error_mean/block_length);
    printf("\n********* final parameters for this block *******\n\n");
    printf("phase error mean(radians):%f\n",phase_error_mean);
    
    //call function for excess noise calcualtion 
    float excessnoise= excess_noise(&alice_x[0], &corrx[0], &shot_x[0]);
    //free heap 
    free (buffer);
    free(buffer1);
    free(theta);
    free(rho);

    return 0;
}

//definition of phase correction function 

int phase_correction(float *bob_phase_sig, float *bob_phase_ref, float *bob_amp_sig, float *corrected_phase, float *corrx, float *corrp){

    //phase correction using ref phase
    *corrected_phase= *bob_phase_sig-*bob_phase_ref;
    //wrap phase to 2Pi
    if(*corrected_phase<0){*corrected_phase= *corrected_phase+2*PI;} 
    //convert to x and p quadratures using corrected phase
    *corrx= (*bob_amp_sig)*cosf(*corrected_phase);
    *corrp= (*bob_amp_sig)*sinf(*corrected_phase);
    return 1;
}

//definition of function for covariance, transmittance, variances and excess noise calculation 

float excess_noise (float *alice_x, float *corrx, float *shot_x){
    float sum_alice=0, sum_bob=0,product_sum=0, cov_alice_bob=0, V_alice=0, Va=2.12, shot_noise=0, shot_noise_sum=0, V_bob=0;
    float det_eff=0.62;//detector efficiency
    float T=0, T_1=0,T_snu=0, excess_noise=0, excess_noise_snu=0; //transmittance and excess noise values 
    float distance=0; //equaivalent distance

    int counter=0;
    while(counter<10000){
        sum_alice=sum_alice+ *(alice_x+counter);
        sum_bob=sum_bob+ *(corrx+counter);
        product_sum=product_sum+ (*(alice_x+counter))*(*(corrx+counter));
        shot_noise_sum=shot_noise_sum+ *(shot_x+counter); //for shot noise variance calculation

        counter=counter+1;
    }
    sum_alice=sum_alice/10000;
    sum_bob=sum_bob/10000;
    product_sum=product_sum/10000;
    cov_alice_bob=product_sum-sum_alice*sum_bob; //covariance beween alice and bob
    shot_noise_sum=shot_noise_sum/10000; //mean of shot noise

    counter=0; //reset counter 
    //shot noise variance, alice variance, and bob variance calculations 
    while(counter<10000){
        shot_noise= shot_noise+ (*(shot_x+counter)-shot_noise_sum)*(*(shot_x+counter)-shot_noise_sum);
        V_alice=V_alice+ (*(alice_x+counter)-sum_alice)*(*(alice_x+counter)-sum_alice);
        V_bob=V_bob+ (*(corrx+counter)-sum_bob)*(*(corrx+counter)-sum_bob);
        counter=counter+1;
    }

    shot_noise=shot_noise/10000; //shoit noise variance 
    V_alice=V_alice/10000;// alice variance from values read from the file
    V_bob=V_bob/10000; //bob variance from measured values at bob
    printf("V_alice=%f\n",V_alice);
    printf("cov=%f\n",cov_alice_bob); 

    //transmittance calculation
    T=cov_alice_bob*cov_alice_bob/(Va*Va*det_eff);
    T_snu=T/(shot_noise);
    T_1=cov_alice_bob*cov_alice_bob/(V_alice*V_alice*det_eff);
    printf("T=%f , T_snu=%f , T_1=%f\n",T,T_snu,T_1);
    //...............NEW COMMENT:transmittance distance caculation..............
    distance=log10f(T_snu)*10/(-0.2);
    printf("distance(km)=%f\n",distance);


    //excess noise calculation 
    excess_noise=(V_bob-det_eff*T_1*V_alice-shot_noise-0.1*shot_noise)/(det_eff*T_1);
    excess_noise_snu=sqrt(excess_noise*excess_noise)/sqrt(shot_noise);
    printf("excess_noise=%f, excess_noise_snu=%f\n",excess_noise,excess_noise_snu);

    return excess_noise_snu;
}