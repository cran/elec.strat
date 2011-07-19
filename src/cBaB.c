/*
 *  lerntobab.h
 *  
 *
 *  Created by Michael Higgins on 5/10/11.
*/

#include <R.h>
#include <Rmath.h>
#include <R_ext/Utils.h>

//Computes the LKPBound
//Requires: 
//m: current number of fixed components
//*cost: the cost vector pointer
//*value: the value vector pointer
//*numUnits: the number of units pointer
//dynThresh: the threshold subtracted by the values of the fixed components
//currmin: the sum of costs of just the fixed components
static double LKPBound(int m, double *cost, double *value, int *numUnits, double dynThresh, double currmin){
		
	//Initialize the lkpBound
	//Other variables:
	//loop: loop = 1, keep looping
	//count: keeps track of number of free components to sum over
	double lkpBound = currmin;
	int loop = 1, count = 0;
	
	//While looping
	while(loop == 1){
		
		//If we have summed over all the units, stop looping
		if(m + count >= *numUnits){
			loop = 0;
		}
		
		else{
			
			//If adding the next unit will push us over the threshold
			if(dynThresh - *(value + m + count) < 0){
				
				//Add fraction of last free cost to LKP Bound, and stop looping
				lkpBound += dynThresh / (*(value + m + count)) * (*(cost + m + count));
				loop = 0;
			}
			
			//If adding next unit doesn't push us over the threshold
			else{
				
				//Add next free cost to LKP bound.
				//Reduce threshold by next free value.
				lkpBound += *(cost + m + count);
				dynThresh -=  *(value + m + count);
			}
		}

		//Increment count
		count++;
	}
	return lkpBound;
}

//minBab is the function, called recursively, used for branching.

//m: number of fixed components
//lambda: stores the exact value of max or min problem
//solution: stores the solution
//value: vector of values 
//cost: vector of costs
//StratID: vector of strat ID's.  Same length as solution, value, and cost.
//numStrat: number of stratum
//numUnits: number of units
//threshold: Threshold on sum of values
//indStrat: indicator to decide whether to branch on a stratum
//worksol: current vector of fixed components that is being evaluated

static void minBab(int m, double *lambda, int *solution, double *value, double *cost, 
				   int *stratID, int *numStrat, int *numUnits, double *threshold, 
				   int *indStrat, int *worksol)
{
	
	//dynThresh: the threshold subtracted by the values of the fixed components
	//currmin: sum of costs of just the fixed components
	//sumIndStrat: sum of indStrat, used to stop branching
	//branch: if branch = 1, keep branching
	//stratIndex: holds value of stratum ID for current fixed component m
	double dynThresh = *threshold;
	double currmin = 0.0;
	int sumIndStrat = 0;
	int branch = 1;
	int stratIndex = 0;
	
	//Check user interrupt
	R_CheckUserInterrupt();
	
	//If there is at least one fixed component
	if(m != 0){
		
		//Subtract threshold by sum of values of fixed components
		//Take sum of costs of fixed components
		for(int i = 0; i < m; i++){
			dynThresh -= (*(value+i)) * (*(worksol+i));
			currmin += (*(cost+i))*(*(worksol+i));
		}
		
		//Take sum of indicator variables indStrat
		for(int i = 0; i < *numStrat; i++){
			sumIndStrat += *(indStrat + i);
		}
		
		//If sum of fixed values greater than the threshold
		if(dynThresh <= 0){
			
			//If sum of fixed costs less than lambda,
			//replace lambda with currmin... 
			//currmin is now the best candidate for the minimum
			//
			//(If sum of fixed costs is greater than lambda,
			//then current vector of fixed components is not optimal)
			if(currmin < *lambda){
				
				//Replace lambda with currmin
				*lambda = currmin;
				
				//Replace current best solution with worksol
				for(int i = 0; i < *numUnits; i++){
					*(solution + i) = *(worksol + i);
				}
			}
			
			//Stop branching from this vector
			branch = 0;
		}
		
		//If sum of the indicator variables is zero
		else if(sumIndStrat == 0){
			
			//Stop branching
			branch = 0;
		}
		
		//If the LKPBound given the current fixed components is greater than lambda
		else if(LKPBound(m, cost, value, numUnits, dynThresh, currmin) >= *lambda){
			
			//No branch from this current fixed vector will be optimal.
			//Stop branching.
			branch = 0;
		}
		
	}
	
	//If we are still branching from this fixed vector
	if(branch == 1){
		
		//Find stratumID corresponding to current fixed component
		stratIndex = *(stratID + m) - 1;
		
		//If the indicator variable indStrat corresponding to
		//the stratumID is 1
		if(*(indStrat + stratIndex) == 1){
			
			//Branch to vector with (m+1)st fixed component set to 1
			worksol[m] = 1;
			minBab(m+1, lambda, solution, 
				   value, cost, stratID, numStrat, numUnits, threshold, 
				   indStrat, worksol);
			
			//Branch to vector with (m+1)st fixed component set to 0
			//Set component of indStrat corresponding to the strata of mth fixed component to 0
			worksol[m] = 0;
			indStrat[stratIndex] = 0;
			minBab(m+1, lambda, solution, 
				   value, cost, stratID, numStrat, numUnits, threshold, 
				   indStrat, worksol);
			
			//Restore value of indStrat
			indStrat[stratIndex] = 1;
		}
		
		//If the indicator variable indStrat corresponding to
		//the stratumID is 1
		else{
			
			//By initialization, worksol[m] = 0
			//Branch to vector with (m+1)st fixed component set to 0
			minBab(m+1, lambda, solution, 
				   value, cost, stratID, numStrat, numUnits, threshold, 
				   indStrat, worksol);
		}
	}
}


//Initializes the BaB algorithm
//Calls the BaB recursive function
//
//lambda: stores the exact value of max or min problem
//solution: stores the solution
//value: vector of values 
//cost: vector of costs
//StratID: vector of strat ID's.  Same length as solution, value, and cost.
//numStrat: number of stratum
//numUnits: number of units
//threshold: Threshold on sum of values
//minimize Not implemented at this time: 1 = minimize, 0 = maximize: 

void BaB(double *lambda, int *solution, 
			   double *value, double *cost, int *stratID, int *numStrat, int *numUnits, double *threshold) 
{	
	//Sanity check removed
	//int sane = 1;
	
	//indStrat: indicator to decide whether to branch on a stratum
	//worksol: current vector of fixed components that is being evaluated
	//m: number of fixed components, will be 0 unless some strata are not sampled.
	int indStrat [*numStrat];
	int worksol [*numUnits]; 
	int m = 0;
	
	//Initialize branch on a stratum indicator vector
	for(int i = 0; i < *numStrat; i++){
		indStrat[i]  = 1;
	}
	
	//Initialize current fixed component vector
	for(int i = 0; i < *numUnits; i++){
		worksol[i]  =  0;
	}
	
	//Are there units with 0 cost?  If so, throw as much error as possible into them.
	//Note that 0 cost stratum have to occur for smallest values of r
	for(int i = 0; i < *numUnits; i++){
		
		//Only 0 samples if the cost is 0.
		if(*(cost + i) == 0.0){
			
			//Include unit in solution
			worksol[i] = 1;
			
			//Set indStrat = 0.  
			//If eventually non-zero cost in that stratum, will make that 1.
			indStrat[*(stratID + i) - 1] = 0;
			
			//Increment m
			m++;	
		}
		
		//If eventually non-zero cost in that stratum, set indStrat to 1
		else{
			indStrat[*(stratID + i) - 1] = 1;
		}
	}
	
	//Sanity check removed
	//
	//Must have minimize either 0 or 1
	//if(*minimize != 0 && *minimize != 1){
	//	Rprintf("Invalid value of minimize.\n");
	//	Rprintf("Set minimize = 1 for minimization or set minimize = 0 for maximization.\n");
	//	sane = 0;
	//}
	//if(*minimize == 1 && sane == 1){
	
	
	//Run BaB function
	minBab(m, lambda, solution, 
			   value, cost, stratID, numStrat, numUnits, threshold, 
			   indStrat, worksol);
	
	//Not implemented at this time
	//if(*minimize == 0 && sane == 1){
	//	Rprintf("maximize it!");
	//	maxBab;
	//}
} 

