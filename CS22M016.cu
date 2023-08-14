#include <iostream>
#include <stdio.h>
#include <cuda.h>

#define max_N 100000
#define max_P 30
#define BLOCKSIZE 1024

using namespace std;

//*******************************************

// Write down the kernels here

// First I am calculating the prefix sum of facids capacity so that i can map centre no and facility ids of the processing request.
// Here i am calculating prefix sum using Reduction method.
// After that I am sending threads which is equal to total no of facility.
// Each thread represing each facility.
// After that I have made worklist  to push all the request which are going to use that facility and slot array  of 25 hours which is intially zero.
// Then I will check whether the facility his free for the that request if it is allowed it will process further otherwise it wont update anything.

__global__ void  exclusivesum(int n, int *garr , int* ccopy,int *csum)
{
    
    if(threadIdx.x+ blockDim.x*blockIdx.x<n){
     garr[threadIdx.x+ blockDim.x*blockIdx.x]-=ccopy[threadIdx.x+ blockDim.x*blockIdx.x];
     if((threadIdx.x+ blockDim.x*blockIdx.x) == n-1)
      csum[0] = ccopy[n-1] + garr[n-1];

    }
}


__global__ void  dkernel(int n,int *a,int start, int end)
{     
     int id=start+threadIdx.x;
  
     if(id<n)
     {
          int  tmp,off=1;
          while(off<n)
        {
            if(threadIdx.x>=off){ 
                tmp=a[id-off];
            __syncthreads();
              a[id]+=tmp;
            __syncthreads();
            }
            int d= off*2;
            off=d;
        } 
    }
}


__global__ void  dke(int L , int K , int C,int *sum, int *gcentre,int *gfacility,int *gcapacity,int *gfac_ids,int *gsucc_reqs,int *gtot_reqs,int *greq_id,int *greq_cen,int *greq_fac,int *greq_start,int *greq_slots,int *gfacps,int R,int N)
{
    int iterations=sum[0], id=threadIdx.x + blockIdx.x * blockDim.x;
    
    if(id < iterations)
    {

    int  counter  = 0;
    int capacity= gcapacity[id];
    int slot[25];
    int centerno=-1;
    int facilityno=-1; 
    int worklist[10000];
    
   
    int ii=0;
    while(ii<25)
    {
      slot[ii] = 0;
      ii++;
    }
    
    int i=0;
    while(i<N)
    { 
      int val=gfacps[i];
      if(id < val)
        {
          centerno = i - 1;
          facilityno =  id - gfacps[i-1];
          break;
        }
        i++;
    }

    if(centerno==-1)
     {
        centerno=N-1;
        facilityno = threadIdx.x + blockIdx.x * blockDim.x - gfacps[centerno];
     }

    int j=0;
    while(j<R)
    {
      if(greq_cen[j] == centerno )
      { 
         if( greq_fac[j] == facilityno){
          worklist[counter] = greq_id[j];
          counter++;
         }
      }
      j++;
    }
    
    
    for(int i=0;i<counter;i++)
    {

        int req=worklist[i], start_slot= greq_start[req], slot_count= greq_slots[req], check=0;
        for(int i=start_slot;i<start_slot+slot_count;i++)
        {
           if(slot[i]>=capacity)
            { 
              check=1;
              break;
            } 
        }

        if(check==0)
           for(int i=start_slot;i<start_slot+slot_count;i++)
             slot[i]=slot[i]+1;

        if(check==0)
           atomicAdd(&gsucc_reqs[centerno],1);
        
    }

    }

    
}
__global__ void  exclusivesum(int n, int *garr , int* ccopy,int *csum)
{
    
    if(threadIdx.x+ blockDim.x*blockIdx.x<n){
     garr[threadIdx.x+ blockDim.x*blockIdx.x]-=ccopy[threadIdx.x+ blockDim.x*blockIdx.x];
     if((threadIdx.x+ blockDim.x*blockIdx.x) == n-1)
      csum[0] = ccopy[n-1] + garr[n-1];

    }
}


__global__ void  dkernel(int n,int *a,int start, int end)
{     
     int id=start+threadIdx.x;
  
     if(id<n)
     {
          int  tmp,off=1;
          while(off<n)
        {
            if(threadIdx.x>=off){ 
                tmp=a[id-off];
            __syncthreads();
              a[id]+=tmp;
            __syncthreads();
            }
            int d= off*2;
            off=d;
        } 
    }
}





//***********************************************


int main(int argc,char **argv)
{
	// variable declarations...
    int N,*centre,*facility,*capacity,*fac_ids, *succ_reqs, *tot_reqs;
    

    FILE *inputfilepointer;
    
    //File Opening for read
    char *inputfilename = argv[1];
    inputfilepointer    = fopen( inputfilename , "r");

    if ( inputfilepointer == NULL )  {
        printf( "input.txt file failed to open." );
        return 0; 
    }

    fscanf( inputfilepointer, "%d", &N ); // N is number of centres
	
    // Allocate memory on cpu
    centre=(int*)malloc(N * sizeof (int));  // Computer  centre numbers
    facility=(int*)malloc(N * sizeof (int));  // Number of facilities in each computer centre
    fac_ids=(int*)malloc(max_P * N  * sizeof (int));  // Facility room numbers of each computer centre
    capacity=(int*)malloc(max_P * N * sizeof (int));  // stores capacities of each facility for every computer centre 


    int success=0;  // total successful requests
    int fail = 0;   // total failed requests
    tot_reqs = (int *)malloc(N*sizeof(int));   // total requests for each centre
    succ_reqs = (int *)malloc(N*sizeof(int)); // total successful requests for each centre

    // Input the computer centres data
    int k1=0 , k2 = 0;
    for(int i=0;i<N;i++)
    {
      fscanf( inputfilepointer, "%d", &centre[i] );
      fscanf( inputfilepointer, "%d", &facility[i] );
      
      for(int j=0;j<facility[i];j++)
      {
        fscanf( inputfilepointer, "%d", &fac_ids[k1] );
        k1++;
      }
      for(int j=0;j<facility[i];j++)
      {
        fscanf( inputfilepointer, "%d", &capacity[k2]);
        k2++;     
      }
    }

    // variable declarations
    int *req_id, *req_cen, *req_fac, *req_start, *req_slots;   // Number of slots requested for every request
    
    // Allocate memory on CPU 
	int R;
	fscanf( inputfilepointer, "%d", &R); // Total requests
    req_id = (int *) malloc ( (R) * sizeof (int) );  // Request ids
    req_cen = (int *) malloc ( (R) * sizeof (int) );  // Requested computer centre
    req_fac = (int *) malloc ( (R) * sizeof (int) );  // Requested facility
    req_start = (int *) malloc ( (R) * sizeof (int) );  // Start slot of every request
    req_slots = (int *) malloc ( (R) * sizeof (int) );   // Number of slots requested for every request
    
    // Input the user request data
    for(int j = 0; j < R; j++)
    {
       fscanf( inputfilepointer, "%d", &req_id[j]);
       fscanf( inputfilepointer, "%d", &req_cen[j]);
       fscanf( inputfilepointer, "%d", &req_fac[j]);
       fscanf( inputfilepointer, "%d", &req_start[j]);
       fscanf( inputfilepointer, "%d", &req_slots[j]);
       tot_reqs[req_cen[j]]+=1;  
    }
		


    //*********************************
    // Call the kernels here
    int *gcentre,*gfacility,*gcapacity,*gfac_ids, *gsucc_reqs, *gtot_reqs,*greq_id, *greq_cen, *greq_fac, *greq_start, *greq_slots;
    cudaMalloc(&gcapacity,sizeof(int)*max_P * N );
    int *csum;
    cudaMalloc(&csum,sizeof(int));
    cudaMalloc(&gfac_ids,sizeof(int)*max_P * N );
    cudaMalloc(&gsucc_reqs,sizeof(int)*N);
    cudaMalloc(&gtot_reqs,sizeof(int)*N);
    cudaMalloc(&greq_id,sizeof(int)*R);
    cudaMalloc(&greq_cen,sizeof(int)*R);
    cudaMalloc(&greq_fac,sizeof(int)*R);
    cudaMalloc(&greq_start,sizeof(int)*R);
    cudaMalloc(&greq_slots,sizeof(int)*R);
    cudaMalloc(&gcentre,sizeof(int)*N);
    cudaMalloc(&gfacility,sizeof(int)*N);
   
    
   
    
    
    cudaMemcpy(greq_id, req_id, R * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(greq_cen, req_cen, R * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(greq_fac, req_fac, R * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(greq_start, req_start, R * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(greq_slots, req_slots, R * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(gcentre,centre , N* sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(gfacility, facility, N * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(gcapacity,capacity , max_P * N * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(gfac_ids, fac_ids, max_P * N * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(gsucc_reqs,succ_reqs ,  N * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(gtot_reqs,tot_reqs, N * sizeof(int), cudaMemcpyHostToDevice);

    cudaMemset(gsucc_reqs, 0, N * sizeof(int));
    int *gfacps,*ccopy;
    cudaMalloc(&gfacps,sizeof(int)*(N));
    cudaMalloc(&ccopy,sizeof(int)*(N));
    cudaMemcpy(gfacps,facility, N * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(ccopy,facility, N * sizeof(int), cudaMemcpyHostToDevice);
    
    int blocksize=1023;
    int noofblocks= ceil((float)N/blocksize);
    

    for(int i=0;i<noofblocks;i++)
    {  
        int start=i*blocksize, end;
        if(i!=noofblocks-1)
           end= (i+1)*blocksize-1;
        else
            end=N-1;

       
       if(i==0)
       dkernel<<<1,blocksize>>>(N,gfacps,start,end);
       else
       dkernel<<<1,blocksize+1>>>(N,gfacps,start-1,end);
       cudaDeviceSynchronize();
    }

    exclusivesum<<<noofblocks,1024>>>(N,gfacps,ccopy,csum);
    int sum[1];
    cudaMemcpy(sum,csum,sizeof(int),cudaMemcpyDeviceToHost);
    int threadsPerBlock = 1024;
    int numBlocks = (sum[0] + threadsPerBlock - 1) / threadsPerBlock;
    
    dke<<<numBlocks, threadsPerBlock>>>(0,0,0,csum, gcentre, gfacility, gcapacity, gfac_ids, gsucc_reqs, gtot_reqs, greq_id, greq_cen, greq_fac, greq_start, greq_slots, gfacps, R, N);
    cudaDeviceSynchronize();

    cudaMemcpy(succ_reqs ,gsucc_reqs,  N * sizeof(int), cudaMemcpyDeviceToHost);
    for(int i=0;i<N;i++)
     success+=succ_reqs[i];

    fail=R-success; 
    
     

  


    //********************************





    // Output
    char *outputfilename = argv[2]; 
    FILE *outputfilepointer;
    outputfilepointer = fopen(outputfilename,"w");

    fprintf( outputfilepointer, "%d %d\n", success, fail);
    //printf("%d %d\n", success, fail);
    for(int j = 0; j < N; j++)
    {
        fprintf( outputfilepointer, "%d %d\n", succ_reqs[j], tot_reqs[j]-succ_reqs[j]);
        //printf("%d %d\n", succ_reqs[j], tot_reqs[j]-succ_reqs[j]);
    }
    fclose( inputfilepointer );
    fclose( outputfilepointer );
    cudaDeviceSynchronize();
	return 0;
}
