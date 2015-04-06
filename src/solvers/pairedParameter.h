#include <vector>
#include <array>
#include "myFun.h"
class PairedParameters
{
public: 
	double c1;
	double c2;
	double c12;
	double c21;
	double c1_plus_c12;
	double c2_plus_c21;
	double a[2];
	double at[2];
	double v_plus[2];
	double v_minus[2];
	double vB_plus[2];
	double vB_minus[2];
	double lambda[2][4];
	double int_lambda[2];
	double d_int_labda[2];//derivative of int_lambda
	double dlambda[2][4];
	double dp[2][4];
	double c1_plus_c12_plus_c2_plus_c21;
	double temp;
	double lambda_plus;
	double lambda_minus;
	double c_plus;
	double c_minus;
	double c_plus_over_at;
	double c_minus_over_at;
	double d_plus;
	double d_minus;
	double detV;
	double v0_1;
	double v1_0;
	double detA;
	double lambdaB_plus;
	double lambdaB_minus;
	double vB0_1;
	double vB0_0;
	double vB1_0;
	double vB1_1;
	double e_lambdaB_plus_t;
	double e_lambdaB_minus_t;
	double lambdaB_e_lambdaB_plus_t;
	double lambdaB_e_lambdaB_minus_t;
	double detVB;
	double cB_plus;
	double cB_minus;
	double v0_1a0;
	double v1_0a0;
	double v1_1a1;
	double v0_0a0;
	double v1_0a1;
	double v0_1a1;
	double vB0_1c1;
	double vB1_0c1;
	double vB1_1c2;
	double vB0_0c1;
	double vB1_0c2;
	double vB0_1c2;
	double e_lambda_plus_t;
	double e_lambda_minus_t;
	double c00;
	double c01;
	double c10;
	double c11;
	double p01_plus_p03;
	double p10_plus_p12;
	double mean;
	double dmean;
	double var;
	double dvar;
	double n0;
	double n1;
	double vB0_1n0;
	double vB1_0n0;
	double vB1_1n0;
	double vB0_0n1;
	double vB0_1n1;
	double vB1_0n1;
	double EX_T[2];
	double varX;
	double p[2][4];
	double lambda2p[2][4];
	static std::size_t NumberOfReactions;
#ifdef DEBUG
	std::vector<double> Propensities;//actually it is jacobian rather than the propensity
	std::vector<int> directions;//-1:backward, 0:not involved, 1:forward
	std::vector<double> forwardPropensities;//actually it is jacobian rather than the propensity
	std::vector<int> forwardReactions;
	std::vector<double> backwardPropensities;//actualy it is jacobian rather than the propensity
	std::vector<int> backwardReactions;
	std::vector<std::vector<int> > inputReactions;//input reactions to the paired species, no catalyst reactions
	std::vector<std::vector<int> > outputReactions;//output reactions from the paired species, no catalyst reactions
	std::vector<std::vector<double> > outputPropensities;//the jacobian of output reactions
#else
	double* Propensities;//actually it is jacobian rather than the propensity
	int* directions;//-1:backward, 0:not involved, 1:forward
	double* forwardPropensities;//actually it is jacobian rather than the propensity
	int* forwardReactions;
	double* backwardPropensities;//actualy it is jacobian rather than the propensity
	int* backwardReactions;
	int** inputReactions;//input reactions to the paired species, no catalyst reactions
	int** outputReactions;//output reactions from the paired species, no catalyst reactions
	double** outputPropensities;//the jacobian of output reactions
#endif
	std::size_t numberOfOutputReactions[2];//outputReactions[i].size()
	std::size_t numberOfInputReactions[2];//inputReactions[i].size()
	std::size_t numberOfForwardReactions;//forwardReactions.size()
	std::size_t numberOfBackwardReactions;//backwardReactions.size()
	//constructor
/*	PairedParameters(int NumberOfReactions):
		p(2, std::vector<double> (4, 0) ),
		lambda2p(2, std::vector<double> (4, 0) ),
		Propensities(std::vector<double> (NumberOfReactions, 0)),
		directions(std::vector<int> (NumberOfReactions, 0)),
		inputReactions(2, std::vector<int> (0)),
		outputReactions(2, std::vector<int> (0)),
		outputPropensities(2, std::vector<double> (0)){}
*/
	
	PairedParameters():numberOfForwardReactions(0), numberOfBackwardReactions(0)
	{
		myNew(inputReactions, 2);
		myNew(outputReactions, 2);
		myNew(outputPropensities, 2);
		myNew(Propensities, NumberOfReactions, 0.0);
		myNew(directions, NumberOfReactions, 0);
	}

	virtual ~PairedParameters()
	{
#ifndef DEBUG
		delete[] inputReactions;
		delete[] outputReactions;
		delete[] outputPropensities;
		delete[] Propensities;
		delete[] directions;
#endif
	}

	void update_full_time_invariant()
	{
		c1_plus_c12=c1+c12;
		c2_plus_c21=c2+c21;
		c1_plus_c12_plus_c2_plus_c21=c1_plus_c12+c2_plus_c21;
		if(c12==0 || c21==0) //just for decrease the number of calculation, no other purpose
		{
			//for matrix A
			if(c1_plus_c12>=c2_plus_c21)
			{
				lambda_plus=c1_plus_c12;
				lambda_minus=c2_plus_c21;
				lambdaB_plus=-c2_plus_c21;
				lambdaB_minus=-c1_plus_c12;
			}
			else
			{
				lambda_plus=c2_plus_c21;
				lambda_minus=c1_plus_c12;
				lambdaB_plus=-c1_plus_c12;
				lambdaB_minus=-c2_plus_c21;
			}
		}
		else
		{
			temp=sqrt(pow(c1_plus_c12-c2_plus_c21, 2)+4*c12*c21);
			//for matrix A
			lambda_plus=(c1_plus_c12_plus_c2_plus_c21+temp)/2.0;
			lambda_minus=(c1_plus_c12_plus_c2_plus_c21-temp)/2.0;
			//for matrix B
			lambdaB_plus=(temp-c1_plus_c12_plus_c2_plus_c21)/2;
			lambdaB_minus=(-c1_plus_c12_plus_c2_plus_c21-temp)/2;
		}
		if(c12!=0)
		{
			v_plus[0]=c12;
			v_plus[1]=c1_plus_c12-lambda_plus;
			v_minus[0]=c12;
			v_minus[1]=c1_plus_c12-lambda_minus;
			vB_plus[0]=c2_plus_c21+lambdaB_plus;
			vB_plus[1]=c12;
			vB_minus[0]=c2_plus_c21+lambdaB_minus;
			vB_minus[1]=c12;
		}
		else if(c21!=0)
		{
			v_plus[0]=c2_plus_c21-lambda_plus;
			v_plus[1]=c21;
			v_minus[0]=c2_plus_c21-lambda_minus;
			v_minus[1]=c21;
			vB_plus[0]=c21;
			vB_plus[1]=c1_plus_c12+lambdaB_plus;
			vB_minus[0]=c21;
			vB_minus[1]=c1_plus_c12+lambdaB_minus;
		}
		else 
		{
			if(c1>=c2)
			{
				v_plus[0]=1;
				v_plus[1]=0;
				v_minus[0]=0;
				v_minus[1]=1;
				vB_plus[0]=0;
				vB_plus[1]=1;
				vB_minus[0]=1;
				vB_minus[1]=0;
			}
			else
			{
				v_plus[0]=0;
				v_plus[1]=1;
				v_minus[0]=1;
				v_minus[1]=0;
				vB_plus[0]=1;
				vB_plus[1]=0;
				vB_minus[0]=0;
				vB_minus[1]=1;
			}
		}
		detV=v_plus[0]*v_minus[1]-v_plus[1]*v_minus[0];
		v0_1=v_plus[0]*v_minus[1]/detV;
		v1_0=v_plus[1]*v_minus[0]/detV;
		v0_0a0=a[0]*v_plus[0]*v_minus[0]/detV;
		v1_1a1=a[1]*v_plus[1]*v_minus[1]/detV;
		v0_1a0=v0_1*a[0];
		v1_0a0=v1_0*a[0];
		v1_0a1=v1_0*a[1];
		v0_1a1=v0_1*a[1];
		detA=c12*c2+c1*c2+c1*c21;
		c00=c1*c2_plus_c21/detA;
		c01=c1*c21/detA;
		c10=c2*c12/detA;
		c11=c2*c1_plus_c12/detA;
		//for matrix B
		//if(c21!=0)
		//{
		//	vB_plus[0]=c21;
		//	vB_plus[1]=c1_plus_c12+lambdaB_plus;
		//	vB_minus[0]=c21;
		//	vB_minus[1]=c1_plus_c12+lambdaB_minus;
		//}
		//else
		//{
		//	vB_plus[0]=c2_plus_c21+lambdaB_plus;
		//	vB_plus[1]=c12;
		//	vB_minus[0]=c2_plus_c21+lambdaB_minus;
		//	vB_minus[1]=c12;
		//}
		detVB=vB_plus[0]*vB_minus[1]-vB_plus[1]*vB_minus[0];
		vB0_1=vB_plus[0]*vB_minus[1]/detVB;
		vB0_0=vB_plus[0]*vB_minus[0]/detVB;
		vB1_0=vB_plus[1]*vB_minus[0]/detVB;
		vB1_1=vB_plus[1]*vB_minus[1]/detVB;
		vB0_1c1=vB0_1*c1;
		vB1_0c1=vB1_0*c1;
		vB1_1c2=vB1_1*c2;
		vB0_0c1=vB0_0*c1;
		vB1_0c2=vB1_0*c2;
		vB0_1c2=vB0_1*c2;
		vB0_1n0=vB0_1*n0;
		vB1_0n0=vB1_0*n0;
		vB1_1n0=vB1_1*n0;
		vB0_0n1=vB0_0*n1;
		vB0_1n1=vB0_1*n1;
		vB1_0n1=vB1_0*n1;

		//int i, n;
		//n=forwardReactions.size();
		//for(i=0; i<n; i++)
		//	Propensities[forwardReactions[i]]=forwardPropensities[i]/c12;
		//n=backwardReactions.size();
		//for(i=0; i<n; i++)
		//	Propensities[backwardReactions[i]]=backwardPropensities[i]/c21;
	}

	void getMean01(double t)//the other species work with 0 and 1 has the opposite mean for their change
	{
		at[0]=a[0]*t;
		at[1]=a[1]*t;
		e_lambda_plus_t=exp(-lambda_plus*t);
		e_lambda_minus_t=exp(-lambda_minus*t);
		c_plus=(1-e_lambda_plus_t)/lambda_plus;//entry of the diagnal matrix
		c_minus=(1-e_lambda_minus_t)/lambda_minus;
		lambda[0][0]=c_plus*v0_1a0-c_minus*v1_0a0;
		lambda[0][1]=(c_minus-c_plus)*v0_0a0;
		lambda[1][0]=(c_plus-c_minus)*v1_1a1;
		lambda[1][1]=c_minus*v0_1a1-c_plus*v1_0a1;
		lambda[0][3]=(at[0]-lambda[0][0])*c10-lambda[0][1]*c11;
		lambda[1][2]=-lambda[1][0]*c00+(at[1]-lambda[1][1])*c01;
		dlambda[0][0]=e_lambda_plus_t*v0_1a0-e_lambda_minus_t*v1_0a0;
		dlambda[0][1]=(e_lambda_minus_t-e_lambda_plus_t)*v0_0a0;
		dlambda[1][0]=(e_lambda_plus_t-e_lambda_minus_t)*v1_1a1;
		dlambda[1][1]=e_lambda_minus_t*v0_1a1-e_lambda_plus_t*v1_0a1;
		dlambda[0][3]=(a[0]-dlambda[0][0])*c10-dlambda[0][1]*c11;
		dlambda[1][2]=-dlambda[1][0]*c00+(a[1]-dlambda[1][1])*c01;

		//for matrix B
		e_lambdaB_plus_t=exp(lambdaB_plus*t);
		e_lambdaB_minus_t=exp(lambdaB_minus*t);
		lambdaB_e_lambdaB_plus_t=lambdaB_plus*e_lambdaB_plus_t;
		lambdaB_e_lambdaB_minus_t=lambdaB_minus*e_lambdaB_minus_t;
		p[0][1]=(e_lambdaB_plus_t-e_lambdaB_minus_t)*vB1_1;
		p[1][0]=(e_lambdaB_minus_t-e_lambdaB_plus_t)*vB0_0;
		cB_plus=(e_lambdaB_plus_t-1)/lambdaB_plus;
		cB_minus=(e_lambdaB_minus_t-1)/lambdaB_minus;
		p[0][3]=(cB_plus-cB_minus)*vB1_1c2;
		p[1][2]=(cB_minus-cB_plus)*vB0_0c1;
		dp[0][1]=(lambdaB_e_lambdaB_plus_t-lambdaB_e_lambdaB_minus_t)*vB1_1;
		dp[1][0]=(lambdaB_e_lambdaB_minus_t-lambdaB_e_lambdaB_plus_t)*vB0_0;
		dp[0][3]=c2*p[0][1];
		dp[1][2]=c1*p[1][0];

		//return
		mean=lambda[1][0]+lambda[1][2]-lambda[0][1]-lambda[0][3]-n0*(p[0][1]+p[0][3])+n1*(p[1][0]+p[1][2]);
		dmean=dlambda[1][0]+dlambda[1][2]-dlambda[0][1]-dlambda[0][3]-n0*(dp[0][1]+dp[0][3])+n1*(dp[1][0]+dp[1][2]);
	}

	void getVar01(double t)//the other species work with 0 and 1 has the same variance for their change
	{
		at[0]=a[0]*t;
		at[1]=a[1]*t;
		e_lambda_plus_t=exp(-lambda_plus*t);
		e_lambda_minus_t=exp(-lambda_minus*t);
		c_plus=(1-e_lambda_plus_t)/lambda_plus;//entry of the diagnal matrix
		c_minus=(1-e_lambda_minus_t)/lambda_minus;
		lambda[0][0]=c_plus*v0_1a0-c_minus*v1_0a0;
		lambda[0][1]=(c_minus-c_plus)*v0_0a0;
		lambda[1][0]=(c_plus-c_minus)*v1_1a1;
		lambda[1][1]=c_minus*v0_1a1-c_plus*v1_0a1;
		lambda[0][3]=(at[0]-lambda[0][0])*c10-lambda[0][1]*c11;
		lambda[1][2]=-lambda[1][0]*c00+(at[1]-lambda[1][1])*c01;
		dlambda[0][0]=e_lambda_plus_t*v0_1a0-e_lambda_minus_t*v1_0a0;
		dlambda[0][1]=(e_lambda_minus_t-e_lambda_plus_t)*v0_0a0;
		dlambda[1][0]=(e_lambda_plus_t-e_lambda_minus_t)*v1_1a1;
		dlambda[1][1]=e_lambda_minus_t*v0_1a1-e_lambda_plus_t*v1_0a1;
		dlambda[0][3]=(a[0]-dlambda[0][0])*c10-dlambda[0][1]*c11;
		dlambda[1][2]=-dlambda[1][0]*c00+(a[1]-dlambda[1][1])*c01;

		//for matrix B
		e_lambdaB_plus_t=exp(lambdaB_plus*t);
		e_lambdaB_minus_t=exp(lambdaB_minus*t);
		lambdaB_e_lambdaB_plus_t=lambdaB_plus*e_lambdaB_plus_t;
		lambdaB_e_lambdaB_minus_t=lambdaB_minus*e_lambdaB_minus_t;
		p[0][1]=(e_lambdaB_plus_t-e_lambdaB_minus_t)*vB1_1;
		p[1][0]=(e_lambdaB_minus_t-e_lambdaB_plus_t)*vB0_0;
		cB_plus=(e_lambdaB_plus_t-1)/lambdaB_plus;
		cB_minus=(e_lambdaB_minus_t-1)/lambdaB_minus;
		p[0][3]=(cB_plus-cB_minus)*vB1_1c2;
		p[1][2]=(cB_minus-cB_plus)*vB0_0c1;
		dp[0][1]=(lambdaB_e_lambdaB_plus_t-lambdaB_e_lambdaB_minus_t)*vB1_1;
		dp[1][0]=(lambdaB_e_lambdaB_minus_t-lambdaB_e_lambdaB_plus_t)*vB0_0;
		dp[0][3]=c2*p[0][1];
		dp[1][2]=c1*p[1][0];

		//return
		p01_plus_p03=p[0][1]+p[0][3];
		p10_plus_p12=p[1][0]+p[1][2];
		var=lambda[1][0]+lambda[1][2]+lambda[0][1]+lambda[0][3]+n0*p01_plus_p03*(1-p01_plus_p03)+n1*p10_plus_p12*(1-p10_plus_p12);
		dvar=dlambda[1][0]+dlambda[1][2]+dlambda[0][1]+dlambda[0][3]+n0*(dp[0][1]+dp[0][3])*(1-2*p01_plus_p03)+n1*(dp[1][0]+dp[1][2])*(1-2*p10_plus_p12);
	}

	void getMean2(double t)
	{
		at[0]=a[0]*t;
		at[1]=a[1]*t;
		e_lambda_plus_t=exp(-lambda_plus*t);
		e_lambda_minus_t=exp(-lambda_minus*t);
		c_plus=(1-e_lambda_plus_t)/lambda_plus;//entry of the diagnal matrix
		c_minus=(1-e_lambda_minus_t)/lambda_minus;
		lambda[0][0]=c_plus*v0_1a0-c_minus*v1_0a0;
		lambda[0][1]=(c_minus-c_plus)*v0_0a0;
		lambda[1][0]=(c_plus-c_minus)*v1_1a1;
		lambda[1][1]=c_minus*v0_1a1-c_plus*v1_0a1;
		lambda[0][2]=(at[0]-lambda[0][0])*c00-lambda[0][1]*c01;
		lambda[1][2]=-lambda[1][0]*c00+(at[1]-lambda[1][1])*c01;
		dlambda[0][0]=e_lambda_plus_t*v0_1a0-e_lambda_minus_t*v1_0a0;
		dlambda[0][1]=(e_lambda_minus_t-e_lambda_plus_t)*v0_0a0;
		dlambda[1][0]=(e_lambda_plus_t-e_lambda_minus_t)*v1_1a1;
		dlambda[1][1]=e_lambda_minus_t*v0_1a1-e_lambda_plus_t*v1_0a1;
		dlambda[0][2]=(a[0]-dlambda[0][0])*c00-dlambda[0][1]*c01;
		dlambda[1][2]=-dlambda[1][0]*c00+(a[1]-dlambda[1][1])*c01;

		//for matrix B
		e_lambdaB_plus_t=exp(lambdaB_plus*t);
		e_lambdaB_minus_t=exp(lambdaB_minus*t);
		p[0][0]=e_lambdaB_plus_t*vB0_1-e_lambdaB_minus_t*vB1_0;
		p[1][0]=(e_lambdaB_minus_t-e_lambdaB_plus_t)*vB0_0;
		cB_plus=(e_lambdaB_plus_t-1)/lambdaB_plus;
		cB_minus=(e_lambdaB_minus_t-1)/lambdaB_minus;
		p[0][2]=cB_plus*vB0_1c1-cB_minus*vB1_0c1;
		p[1][2]=(cB_minus-cB_plus)*vB0_0c1;
		dp[0][2]=c1*p[0][0];
		dp[1][2]=c1*p[1][0];

		//return
		mean=lambda[0][2]+lambda[1][2]+n0*p[0][2]+n1*p[1][2];
		dmean=dlambda[0][2]+dlambda[1][2]+n0*dp[0][2]+n1*dp[1][2];
	}

	void getVar2(double t)
	{
		at[0]=a[0]*t;
		at[1]=a[1]*t;
		e_lambda_plus_t=exp(-lambda_plus*t);
		e_lambda_minus_t=exp(-lambda_minus*t);
		c_plus=(1-e_lambda_plus_t)/lambda_plus;//entry of the diagnal matrix
		c_minus=(1-e_lambda_minus_t)/lambda_minus;
		lambda[0][0]=c_plus*v0_1a0-c_minus*v1_0a0;
		lambda[0][1]=(c_minus-c_plus)*v0_0a0;
		lambda[1][0]=(c_plus-c_minus)*v1_1a1;
		lambda[1][1]=c_minus*v0_1a1-c_plus*v1_0a1;
		lambda[0][2]=(at[0]-lambda[0][0])*c00-lambda[0][1]*c01;
		lambda[1][2]=-lambda[1][0]*c00+(at[1]-lambda[1][1])*c01;
		dlambda[0][0]=e_lambda_plus_t*v0_1a0-e_lambda_minus_t*v1_0a0;
		dlambda[0][1]=(e_lambda_minus_t-e_lambda_plus_t)*v0_0a0;
		dlambda[1][0]=(e_lambda_plus_t-e_lambda_minus_t)*v1_1a1;
		dlambda[1][1]=e_lambda_minus_t*v0_1a1-e_lambda_plus_t*v1_0a1;
		dlambda[0][2]=(a[0]-dlambda[0][0])*c00-dlambda[0][1]*c01;
		dlambda[1][2]=-dlambda[1][0]*c00+(a[1]-dlambda[1][1])*c01;

		//for matrix B
		e_lambdaB_plus_t=exp(lambdaB_plus*t);
		e_lambdaB_minus_t=exp(lambdaB_minus*t);
		p[0][0]=e_lambdaB_plus_t*vB0_1-e_lambdaB_minus_t*vB1_0;
		p[1][0]=(e_lambdaB_minus_t-e_lambdaB_plus_t)*vB0_0;
		cB_plus=(e_lambdaB_plus_t-1)/lambdaB_plus;
		cB_minus=(e_lambdaB_minus_t-1)/lambdaB_minus;
		p[0][2]=cB_plus*vB0_1c1-cB_minus*vB1_0c1;
		p[1][2]=(cB_minus-cB_plus)*vB0_0c1;
		dp[0][2]=c1*p[0][0];
		dp[1][2]=c1*p[1][0];

		//return
		var=lambda[0][2]+lambda[1][2]+n0*p[0][2]*(1-p[0][2])+n1*p[1][2]*(1-p[1][2]);
		dvar=dlambda[0][2]+dlambda[1][2]+n0*(1-2*p[0][2])*dp[0][2]+n1*(1-2*p[1][2])*dp[1][2];
	}

	void getMean3(double t)
	{
		at[0]=a[0]*t;
		at[1]=a[1]*t;
		e_lambda_plus_t=exp(-lambda_plus*t);
		e_lambda_minus_t=exp(-lambda_minus*t);
		c_plus=(1-e_lambda_plus_t)/lambda_plus;//entry of the diagnal matrix
		c_minus=(1-e_lambda_minus_t)/lambda_minus;
		lambda[0][0]=c_plus*v0_1a0-c_minus*v1_0a0;
		lambda[0][1]=(c_minus-c_plus)*v0_0a0;
		lambda[1][0]=(c_plus-c_minus)*v1_1a1;
		lambda[1][1]=c_minus*v0_1a1-c_plus*v1_0a1;
		lambda[0][3]=(at[0]-lambda[0][0])*c10-lambda[0][1]*c11;
		lambda[1][3]=-lambda[1][0]*c10+(at[1]-lambda[1][1])*c11;
		dlambda[0][0]=e_lambda_plus_t*v0_1a0-e_lambda_minus_t*v1_0a0;
		dlambda[0][1]=(e_lambda_minus_t-e_lambda_plus_t)*v0_0a0;
		dlambda[1][0]=(e_lambda_plus_t-e_lambda_minus_t)*v1_1a1;
		dlambda[1][1]=e_lambda_minus_t*v0_1a1-e_lambda_plus_t*v1_0a1;
		dlambda[0][3]=(a[0]-dlambda[0][0])*c10-dlambda[0][1]*c11;
		dlambda[1][3]=-dlambda[1][0]*c10+(a[1]-dlambda[1][1])*c11;

		//for matrix B
		e_lambdaB_plus_t=exp(lambdaB_plus*t);
		e_lambdaB_minus_t=exp(lambdaB_minus*t);
		p[0][1]=(e_lambdaB_plus_t-e_lambdaB_minus_t)*vB1_1;
		p[1][1]=e_lambdaB_minus_t*vB0_1-e_lambdaB_plus_t*vB1_0;
		cB_plus=(e_lambdaB_plus_t-1)/lambdaB_plus;
		cB_minus=(e_lambdaB_minus_t-1)/lambdaB_minus;
		p[0][3]=(cB_plus-cB_minus)*vB1_1c2;
		p[1][3]=cB_minus*vB0_1c2-cB_plus*vB1_0c2;
		dp[0][3]=c2*p[0][1];
		dp[1][3]=c2*p[1][1];

		//return
		mean=lambda[0][3]+lambda[1][3]+n0*p[0][3]+n1*p[1][3];
		dmean=dlambda[0][3]+dlambda[1][3]+n0*dp[0][3]+n1*dp[1][3];
	}

	void getVar3(double t)
	{
		at[0]=a[0]*t;
		at[1]=a[1]*t;
		e_lambda_plus_t=exp(-lambda_plus*t);
		e_lambda_minus_t=exp(-lambda_minus*t);
		c_plus=(1-e_lambda_plus_t)/lambda_plus;//entry of the diagnal matrix
		c_minus=(1-e_lambda_minus_t)/lambda_minus;
		lambda[0][0]=c_plus*v0_1a0-c_minus*v1_0a0;
		lambda[0][1]=(c_minus-c_plus)*v0_0a0;
		lambda[1][0]=(c_plus-c_minus)*v1_1a1;
		lambda[1][1]=c_minus*v0_1a1-c_plus*v1_0a1;
		lambda[0][3]=(at[0]-lambda[0][0])*c10-lambda[0][1]*c11;
		lambda[1][3]=-lambda[1][0]*c10+(at[1]-lambda[1][1])*c11;
		dlambda[0][0]=e_lambda_plus_t*v0_1a0-e_lambda_minus_t*v1_0a0;
		dlambda[0][1]=(e_lambda_minus_t-e_lambda_plus_t)*v0_0a0;
		dlambda[1][0]=(e_lambda_plus_t-e_lambda_minus_t)*v1_1a1;
		dlambda[1][1]=e_lambda_minus_t*v0_1a1-e_lambda_plus_t*v1_0a1;
		dlambda[0][3]=(a[0]-dlambda[0][0])*c10-dlambda[0][1]*c11;
		dlambda[1][3]=-dlambda[1][0]*c10+(a[1]-dlambda[1][1])*c11;

		//for matrix B
		e_lambdaB_plus_t=exp(lambdaB_plus*t);
		e_lambdaB_minus_t=exp(lambdaB_minus*t);
		p[0][1]=(e_lambdaB_plus_t-e_lambdaB_minus_t)*vB1_1;
		p[1][1]=e_lambdaB_minus_t*vB0_1-e_lambdaB_plus_t*vB1_0;
		cB_plus=(e_lambdaB_plus_t-1)/lambdaB_plus;
		cB_minus=(e_lambdaB_minus_t-1)/lambdaB_minus;
		p[0][3]=(cB_plus-cB_minus)*vB1_1c2;
		p[1][3]=cB_minus*vB0_1c2-cB_plus*vB1_0c2;
		dp[0][3]=c2*p[0][1];
		dp[1][3]=c2*p[1][1];

		//return
		var=lambda[0][3]+lambda[1][3]+n0*p[0][3]*(1-p[0][3])+n1*p[1][3]*(1-p[1][3]);
		dvar=dlambda[0][3]+dlambda[1][3]+n0*(1-2*p[0][3])*dp[0][3]+n1*(1-2*p[1][3])*dp[1][3];
	}

	void getMean(double t, int reaction)//the expected number of firings for a reaction
	{
		e_lambda_plus_t=exp(-lambda_plus*t);
		e_lambda_minus_t=exp(-lambda_minus*t);
		c_plus=(1-e_lambda_plus_t)/lambda_plus;//entry of the diagnal matrix
		c_minus=(1-e_lambda_minus_t)/lambda_minus;
		d_plus=(t-c_plus)/lambda_plus;//diagnal entry for the integral of lambda
		d_minus=(t-c_minus)/lambda_minus;

		//for matrix B
		e_lambdaB_plus_t=exp(lambdaB_plus*t);
		e_lambdaB_minus_t=exp(lambdaB_minus*t);
		lambdaB_e_lambdaB_plus_t=lambdaB_plus*e_lambdaB_plus_t;
		lambdaB_e_lambdaB_minus_t=lambdaB_minus*e_lambdaB_minus_t;
		cB_plus=(e_lambdaB_plus_t-1)/lambdaB_plus;
		cB_minus=(e_lambdaB_minus_t-1)/lambdaB_minus;

		if(directions[reaction]==1)//forward reaction only s0 is interested
		{
			int_lambda[0]=d_plus*v0_1a0-d_minus*v1_0a0;
			int_lambda[0]+=(d_plus-d_minus)*v1_1a1;
			d_int_labda[0]=c_plus*v0_1a0-c_minus*v1_0a0;
			d_int_labda[0]+=(c_plus-c_minus)*v1_1a1;

			int_lambda[0]+=cB_plus*vB0_1n0-cB_minus*vB1_0n0;
			int_lambda[0]+=(cB_minus-cB_plus)*vB0_0n1;
			d_int_labda[0]+=e_lambdaB_plus_t*vB0_1n0-e_lambdaB_minus_t*vB1_0n0;
			d_int_labda[0]+=(e_lambdaB_minus_t-e_lambdaB_plus_t)*vB0_0n1;

			//return
			mean=int_lambda[0]*Propensities[reaction];
			dmean=d_int_labda[0]*Propensities[reaction];
			
			lambda[0][0]=c_plus*v0_1a0-c_minus*v1_0a0;
			lambda[1][0]=(c_plus-c_minus)*v1_1a1;
			p[0][0]=e_lambdaB_plus_t*vB0_1-e_lambdaB_minus_t*vB1_0;
			p[1][0]=(e_lambdaB_minus_t-e_lambdaB_plus_t)*vB0_0;
			dlambda[0][0]=e_lambda_plus_t*v0_1a0-e_lambda_minus_t*v1_0a0;
			dlambda[1][0]=(e_lambda_plus_t-e_lambda_minus_t)*v1_1a1;
			dp[0][0]=lambdaB_e_lambdaB_plus_t*vB0_1-lambdaB_e_lambdaB_minus_t*vB1_0;
			dp[1][0]=(lambdaB_e_lambdaB_minus_t-lambdaB_e_lambdaB_plus_t)*vB0_0;
			varX=lambda[0][0]+lambda[1][0]+n0*p[0][0]*(1-p[0][0])+n1*p[1][0]*(1-p[1][0]);
			temp=pow(t*Propensities[reaction]/2, 2);
			var=mean+temp*varX;
			dvar=dmean+temp*2/t*varX+temp*(dlambda[0][0]+dlambda[1][0]+n0*dp[0][0]*(1-2*p[0][0])+n1*dp[1][0]*(1-2*p[1][0]));
		}
		else
		{
			int_lambda[1]=(d_minus-d_plus)*v0_0a0;
			int_lambda[1]+=d_minus*v0_1a1-d_plus*v1_0a1;
			d_int_labda[1]=(c_minus-c_plus)*v0_0a0;
			d_int_labda[1]+=c_minus*v0_1a1-c_plus*v1_0a1;

			int_lambda[1]+=(cB_plus-cB_minus)*vB1_1n0;
			int_lambda[1]+=cB_minus*vB0_1n1-cB_plus*vB1_0n1;
			d_int_labda[1]+=(e_lambdaB_plus_t-e_lambdaB_minus_t)*vB1_1n0;
			d_int_labda[1]+=e_lambdaB_minus_t*vB0_1n1-e_lambdaB_plus_t*vB1_0n1;

			//return
			mean=int_lambda[1]*Propensities[reaction];
			dmean=d_int_labda[1]*Propensities[reaction];

			lambda[0][1]=(c_minus-c_plus)*v0_0a0;
			lambda[1][1]=c_minus*v0_1a1-c_plus*v1_0a1;
			p[0][1]=(e_lambdaB_plus_t-e_lambdaB_minus_t)*vB1_1;
			p[1][1]=e_lambdaB_minus_t*vB0_1-e_lambdaB_plus_t*vB1_0;
			dlambda[0][1]=(e_lambda_minus_t-e_lambda_plus_t)*v0_0a0;
			dlambda[1][1]=e_lambda_minus_t*v0_1a1-e_lambda_plus_t*v1_0a1;
			dp[0][1]=(lambdaB_e_lambdaB_plus_t-lambdaB_e_lambdaB_minus_t)*vB1_1;
			dp[1][1]=lambdaB_e_lambdaB_minus_t*vB0_1-lambdaB_e_lambdaB_plus_t*vB1_0;
			varX=lambda[0][1]+lambda[1][1]+n0*p[0][1]*(1-p[0][1])+n1*p[1][1]*(1-p[1][1]);
			temp=pow(t*Propensities[reaction]/2, 2);
			var=mean+temp*varX;
			dvar=dmean+temp*2/t*varX+temp*(dlambda[0][1]+dlambda[1][1]+n0*dp[0][1]*(1-2*p[0][1])+n1*dp[1][1]*(1-2*p[1][1]));
		}
	}

	void update_full_time_dependent(double t)
	{
		int i;
		//at[0]=a[0]*t;
		//at[1]=a[1]*t;
		e_lambda_plus_t=exp(-lambda_plus*t);
		e_lambda_minus_t=exp(-lambda_minus*t);
		c_plus=(1-e_lambda_plus_t)/lambda_plus;//entry of the diagnal matrix
		c_minus=(1-e_lambda_minus_t)/lambda_minus;
		//lambda[0][0]=c_plus*v0_1a0-c_minus*v1_0a0;
		//lambda[0][1]=(c_minus-c_plus)*v0_0a0;
		//lambda[1][0]=(c_plus-c_minus)*v1_1a1;
		//lambda[1][1]=c_minus*v0_1a1-c_plus*v1_0a1;
		//lambda[0][2]=(at[0]-lambda[0][0])*c00-lambda[0][1]*c01;
		//lambda[0][3]=(at[0]-lambda[0][0])*c10-lambda[0][1]*c11;
		//lambda[1][2]=-lambda[1][0]*c00+(at[1]-lambda[1][1])*c01;
		//lambda[1][3]=-lambda[1][0]*c10+(at[1]-lambda[1][1])*c11;

		////lambda to p
		//for(int i=0;i<2;i++)
		//{
		//	if(at[i]!=0)
		//	{
		//		for(int j=0;j<4;j++)
		//			lambda2p[i][j]=lambda[i][j]/at[i];
		//	}
		//}
		if(a[0]==0)
		{
			for(i=0; i<4; i++)
			{
				lambda[0][i]=0;
				lambda2p[0][i]=0;
			}
		}
		else
		{
			at[0]=a[0]*t;
			c_plus_over_at=c_plus/at[0];
			c_minus_over_at=c_minus/at[0];
			lambda2p[0][0]=c_plus_over_at*v0_1a0-c_minus_over_at*v1_0a0;
			lambda2p[0][1]=(c_minus_over_at-c_plus_over_at)*v0_0a0;
			lambda2p[0][2]=(1-lambda2p[0][0])*c00-lambda2p[0][1]*c01;
			lambda2p[0][3]=(1-lambda2p[0][0])*c10-lambda2p[0][1]*c11;
			for(int j=0;j<4;j++)
				lambda[0][j]=lambda2p[0][j]*at[0];
		}

		if(a[1]==0)
		{
			for(i=0; i<4; i++)
			{
				lambda[1][i]=0;
				lambda2p[1][i]=0;
			}
		}
		else
		{
			at[1]=a[1]*t;
			c_plus_over_at=c_plus/at[1];
			c_minus_over_at=c_minus/at[1];
			lambda2p[1][0]=(c_plus_over_at-c_minus_over_at)*v1_1a1;
			lambda2p[1][1]=c_minus_over_at*v0_1a1-c_plus_over_at*v1_0a1;
			lambda2p[1][2]=-lambda2p[1][0]*c00+(1-lambda2p[1][1])*c01;
			lambda2p[1][3]=-lambda2p[1][0]*c10+(1-lambda2p[1][1])*c11;
			for(int j=0;j<4;j++)
				lambda[1][j]=lambda2p[1][j]*at[1];
		}

		//for matrix B
		e_lambdaB_plus_t=exp(lambdaB_plus*t);
		e_lambdaB_minus_t=exp(lambdaB_minus*t);
		lambdaB_e_lambdaB_plus_t=lambdaB_plus*e_lambdaB_plus_t;
		lambdaB_e_lambdaB_minus_t=lambdaB_minus*e_lambdaB_minus_t;
		p[0][0]=e_lambdaB_plus_t*vB0_1-e_lambdaB_minus_t*vB1_0;
		p[0][1]=(e_lambdaB_plus_t-e_lambdaB_minus_t)*vB1_1;
		p[1][0]=(e_lambdaB_minus_t-e_lambdaB_plus_t)*vB0_0;
		p[1][1]=e_lambdaB_minus_t*vB0_1-e_lambdaB_plus_t*vB1_0;
		cB_plus=(e_lambdaB_plus_t-1)/lambdaB_plus;
		cB_minus=(e_lambdaB_minus_t-1)/lambdaB_minus;
		p[0][2]=cB_plus*vB0_1c1-cB_minus*vB1_0c1;
		p[0][3]=(cB_plus-cB_minus)*vB1_1c2;
		p[1][2]=(cB_minus-cB_plus)*vB0_0c1;
		p[1][3]=cB_minus*vB0_1c2-cB_plus*vB1_0c2;
		
		d_plus=(t-c_plus)/lambda_plus;//diagnal entry for the integral of lambda
		d_minus=(t-c_minus)/lambda_minus;
		int_lambda[0]=d_plus*v0_1a0-d_minus*v1_0a0;
		int_lambda[0]+=(d_plus-d_minus)*v1_1a1;
		int_lambda[0]+=cB_plus*vB0_1n0-cB_minus*vB1_0n0;
		int_lambda[0]+=(cB_minus-cB_plus)*vB0_0n1;
		int_lambda[1]=(d_minus-d_plus)*v0_0a0;
		int_lambda[1]+=d_minus*v0_1a1-d_plus*v1_0a1;
		int_lambda[1]+=(cB_plus-cB_minus)*vB1_1n0;
		int_lambda[1]+=cB_minus*vB0_1n1-cB_plus*vB1_0n1;
		EX_T[0]=n0*p[0][0]+n1*p[1][0]+lambda[0][0]+lambda[1][0];
		EX_T[1]=n0*p[0][1]+n1*p[1][1]+lambda[0][1]+lambda[1][1];
	}
/*	void update_full_time_dependent(double t)
	{
		at[0]=a[0]*t;
		at[1]=a[1]*t;
		e_lambda_plus_t=exp(-lambda_plus*t);
		e_lambda_minus_t=exp(-lambda_minus*t);
		c_plus=(1-e_lambda_plus_t)/lambda_plus;//entry of the diagnal matrix
		c_minus=(1-e_lambda_minus_t)/lambda_minus;
		lambda[0][0]=c_plus*v0_1a0-c_minus*v1_0a0;
		lambda[0][1]=(c_minus-c_plus)*v0_0a0;
		lambda[1][0]=(c_plus-c_minus)*v1_1a1;
		lambda[1][1]=c_minus*v0_1a1-c_plus*v1_0a1;
		lambda[0][2]=(at[0]-lambda[0][0])*c00-lambda[0][1]*c01;
		lambda[0][3]=(at[0]-lambda[0][0])*c10-lambda[0][1]*c11;
		lambda[1][2]=-lambda[1][0]*c00+(at[1]-lambda[1][1])*c01;
		lambda[1][3]=-lambda[1][0]*c10+(at[1]-lambda[1][1])*c11;
		dlambda[0][0]=e_lambda_plus_t*v0_1a0-e_lambda_minus_t*v1_0a0;
		dlambda[0][1]=(e_lambda_minus_t-e_lambda_plus_t)*v0_0a0;
		dlambda[1][0]=(e_lambda_plus_t-e_lambda_minus_t)*v1_1a1;
		dlambda[1][1]=e_lambda_minus_t*v0_1a1-e_lambda_plus_t*v1_0a1;
		dlambda[0][2]=(a[0]-dlambda[0][0])*c00-dlambda[0][1]*c01;
		dlambda[0][3]=(a[0]-dlambda[0][0])*c10-dlambda[0][1]*c11;
		dlambda[1][2]=-dlambda[1][0]*c00+(a[1]-dlambda[1][1])*c01;
		dlambda[1][3]=-dlambda[1][0]*c10+(a[1]-dlambda[1][1])*c11;

		//for matrix B
		e_lambdaB_plus_t=exp(lambdaB_plus*t);
		e_lambdaB_minus_t=exp(lambdaB_minus*t);
		lambdaB_e_lambdaB_plus_t=lambdaB_plus*e_lambdaB_plus_t;
		lambdaB_e_lambdaB_minus_t=lambdaB_minus*e_lambdaB_minus_t;
		p[0][0]=e_lambdaB_plus_t*vB0_1-e_lambdaB_minus_t*vB1_0;
		p[0][1]=(e_lambdaB_plus_t-e_lambdaB_minus_t)*vB1_1;
		p[1][0]=(e_lambdaB_minus_t-e_lambdaB_plus_t)*vB0_0;
		p[1][1]=e_lambdaB_minus_t*vB0_1-e_lambdaB_plus_t*vB1_0;
		cB_plus=(e_lambdaB_plus_t-1)/lambdaB_plus;
		cB_minus=(e_lambdaB_minus_t-1)/lambdaB_minus;
		p[0][2]=cB_plus*vB0_1c1-cB_minus*vB1_0c2;
		p[0][3]=(cB_plus-cB_minus)*vB1_1c2;
		p[1][2]=(cB_minus-cB_plus)*vB0_0c1;
		p[1][3]=cB_minus*vB0_1c2-cB_plus*vB1_0c2;
		dp[0][0]=lambdaB_e_lambdaB_plus_t*vB0_1-lambdaB_e_lambdaB_minus_t*vB1_0;
		dp[0][1]=(lambdaB_e_lambdaB_plus_t-lambdaB_e_lambdaB_minus_t)*vB1_1;
		dp[1][0]=(lambdaB_e_lambdaB_minus_t-lambdaB_e_lambdaB_plus_t)*vB0_0;
		dp[1][1]=lambdaB_e_lambdaB_minus_t*vB0_1-lambdaB_e_lambdaB_plus_t*vB1_0;
		dp[0][2]=c1*p[0][0];
		dp[0][3]=c2*p[0][1];
		dp[1][2]=c1*p[1][0];
		dp[1][3]=c2*p[1][1];
	}*/
};

std::size_t PairedParameters::NumberOfReactions=0;