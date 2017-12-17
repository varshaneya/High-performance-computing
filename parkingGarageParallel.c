#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<mpi.h>
#include<omp.h>

#define BLOCK_LOW(id,p,n)((id)*(n)/(p))
#define BLOCK_HIGH(id,p,n)(BLOCK_LOW(id+1,p,n) - 1)
#define BLOCK_SIZE(id,p,n)(BLOCK_LOW(id+1,p,n) - BLOCK_LOW(id,p,n))

typedef struct garage
{
    int occupied;
    double TTL;
}garage;

int doubleComparator(const void *a,const void *b)
{
    double *c = (double*)a;
    double *d = (double*)b;
        
    if(*c < *d)
        return -1;
    else if(*c > *d)
        return 1;
    else
        return 0;
}

double getNormalVariate(double mu,double sigma)
{
    static int hasNormalVariate=0;
    static double nv2;
    
    if(hasNormalVariate)
    {
        hasNormalVariate=0;
        return nv2;
    }
    
    double v1,v2,r,f,nv1;

    do
    {
        v1 = 2.0 * ((double)rand()/RAND_MAX) -1.0;
        v2 = 2.0 * ((double)rand()/RAND_MAX) -1.0;
        r = v1*v1 + v2*v2;
    }while(r<=0.0 ||  r>=1.0);
    
    f = sqrt(-2.0*log(r)/r);
    nv1 = fabs((sigma*f*v1)+mu);
    nv2 = fabs((sigma*f*v2)+mu);
    hasNormalVariate=1;
    
    return nv1;
}

void generateInterArrivalTimes(double interArrivalTimes[],double A,int size,int argc,char*argv[])
{

    int myrank,numProcs;
    double *sendBuf;
    int i;
    int sendCnt;
    int *recvCnt=NULL,*recvDisp=NULL;

    MPI_Init(&argc,&argv);        
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    MPI_Comm_size(MPI_COMM_WORLD,&numProcs);
    
    srand(time(NULL)+(unsigned int)myrank);
    
    sendCnt = BLOCK_SIZE(myrank,numProcs,size);
    sendBuf = (double*)malloc(sizeof(double)*sendCnt);

    if(myrank == 0)
    {
        recvCnt = malloc(sizeof(int)*numProcs);
        recvDisp = malloc(sizeof(int)*numProcs);
        int count;
        
        for(i=0,count=0;i<numProcs;i++)
        {
            recvDisp[i] = count;
            recvCnt[i]=BLOCK_SIZE(i,numProcs,size);
            count += recvCnt[i];
        }
            
        sendBuf[0]=0.0;
        
        omp_set_num_threads(omp_get_num_procs());
        #pragma parallel for
        for(i=1;i<sendCnt;i++)
            sendBuf[i] = -1*A*log((double)rand()/RAND_MAX);
    }
    
    else
    {
        omp_set_num_threads(omp_get_num_procs());
        #pragma parallel for
        for(i=0;i<sendCnt;i++)
            sendBuf[i] = -1*A*log((double)rand()/RAND_MAX);
    }
    
    MPI_Gatherv(sendBuf,sendCnt,MPI_DOUBLE,interArrivalTimes,recvCnt,recvDisp,MPI_DOUBLE,0,MPI_COMM_WORLD);
    
    free(sendBuf);

    if(myrank == 0)
    {
        free(recvCnt);
        free(recvDisp);
    }
    
    MPI_Finalize();
    
    if(myrank == 0)
        return;
    else
        exit(0);
}

void initialize(garage stalls[],int size) 
{
    int i;
    omp_set_num_threads(omp_get_num_procs());
    #pragma omp parallel for
    for(i=0;i<size;i++)
        stalls[i].occupied = 0;
}

int checkStatus(garage stalls[],int size,double currentTime)
{
    int i;
    int found = -1;
    
    for(i=0;i<size;i++)
        if(stalls[i].occupied == 0)
            found = i;
        else if(stalls[i].occupied == 1 && stalls[i].TTL <= currentTime)
        {
            stalls[i].occupied = 0;
            stalls[i].TTL=0.0;
            if(found == -1)
                found = i;
        }
    return found;
}

int assignCar(garage stalls[],int size,double currentTime,double mu,double sigma)
{
    int status = checkStatus(stalls,size,currentTime);
    
    if(status != -1)    //when an empty space is found
    {
        stalls[status].occupied = 1;
        stalls[status].TTL = currentTime + getNormalVariate(mu,sigma);
    }
    
    return status;
}

int calNumStallsOccupied(garage stalls[],int size)
{
    int i;
    int occupied=0;
    
    for(i=0;i<size;i++)
        if(stalls[i].occupied == 1)
            ++occupied;
    
   return occupied;
        
}

int main(int argc,char* argv[])
{
    if(argc != 5)
    {
        printf("usage ./parkingGarage <simulation size> <no of stalls in the garage> <mean time between car arrivals> \
        <mean time of stay in the garage>");
        return -1;
    }
    
    int IATsize = atoi(argv[1]);  //simulation size
    int numStalls = atoi(argv[2]);    //no of stalls in the garage
    double A = atoi(argv[3]);   //mean time between car arrivals
    double M = atoi(argv[4]);  //mean time of stay in the garage
    
    
    int i,status,carsTurnedAway=0;
    long numStallsOccupied=0;
    double *interArrivalTimes = (double*)malloc(IATsize*sizeof(double));
    
    if(!interArrivalTimes)
    {
        printf("No memory for simulation. Quitting\n");
        return -1;
    }
    
    garage* stalls = (garage*)malloc(numStalls*sizeof(garage));
    
    if(!stalls)
    {
        printf("No memory for simulation. Quitting\n");
        return -1;
    }
    
    omp_set_num_threads(omp_get_num_procs());
    
    #pragma omp parallel
    {
        #pragma parallel sections
        {
            #pragma parallel section
            generateInterArrivalTimes(interArrivalTimes,A,IATsize,argc,argv);
            #pragma parallel section    
            initialize(stalls,numStalls);
        }
    }
    double currentTime=0.0;
    
    for(i=0;i<IATsize;i++)
    {
        currentTime += interArrivalTimes[i];
        if((status=assignCar(stalls,numStalls,currentTime,M,M/4))== -1)
            ++carsTurnedAway;
        numStallsOccupied+=calNumStallsOccupied(stalls,numStalls);
    }

   printf("Average no of stalls occupied: %lf\n",(double)numStallsOccupied/IATsize);
   printf("probabiltiy of a car being turned away due to full parking space: %lf\n",(double)carsTurnedAway/IATsize);

    return 0;
}
