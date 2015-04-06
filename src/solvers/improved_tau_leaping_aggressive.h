/*!
\brief the explicit tau-leaping with dynamic step size selection that switches to SSA when step size is small
*/

#ifndef _IMPROVED_TAU_LEAPING_AGGRESSIVE_H_
#define _IMPROVED_TAU_LEAPING_AGGRESSIVE_H_

#include <iostream>
#include <vector>
#include <deque>
#include <list>
#include <limits>
#include <algorithm>
#include <cmath>
#include <functional>
#include <boost/bind/bind.hpp>
#include "Random.h"
#include "StandardDriverTypes.h"
#include "SSA_Direct.h"
#include "TauLeapingExplicitAdaptive.h"
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include "pairedParameter.h"
#include "myFun.h"
/*! 
\param _denseVectorType the population vector type, should be dense
\param _matrixType
\param _propensitiesFunctorType functor takes reaction index and _denseVectorType population
and returns the propensity for that reaction
\param _dependencyGraphType
*/
/*class Jacobian
{
public:
Jacobian(int the_index, const std::vector<int> &the_reactant, double (*the_function)(int, int)):index(the_index), reactant(the_reactant), function(the_function){}

double getJacobian(int i)//i:species index
{
if(i==reactant[0])
return function(reactant[1], index);
else if(i==reactant[1])
return function(reactant[0], index);
else
return 0;
}
private:
int index;
std::vector<int> reactant;
double (*function)(int, int);
};*/
namespace STOCHKIT
{
	template<typename _denseVectorType, 
		typename _matrixType,
		typename _propensitiesFunctorType,
		typename _dependencyGraphType>
	class improvedTauLeaping_aggressive : public TauLeapingExplicitAdaptive<_denseVectorType, _matrixType, _propensitiesFunctorType, _dependencyGraphType>
	{	
	public:
		typedef TauLeapingExplicitAdaptive<_denseVectorType, _matrixType, _propensitiesFunctorType, _dependencyGraphType> TauLeaping;
		typedef TauLeaping::SSA_Direct<_denseVectorType, _matrixType, _propensitiesFunctorType, _dependencyGraphType> SSA;
		using SSA::currentTime;
		using SSA::currentPopulation;
		using SSA::currentPropensities;
		using SSA::calculateAllPropensities;
		using SSA::NumberOfReactions;
		using SSA::NumberOfSpecies;
		using SSA::stoichiometry;
		using SSA::propensities;
		using SSA::randomGenerator;
		using SSA::initialPopulation;
		using TauLeaping::criticalPropensitySum;
		using TauLeaping::g;
		using TauLeaping::speciesToReaction;
		using TauLeaping::affectedReactions;
		using TauLeaping::NumberOfReactants;
		using TauLeaping::SSASteps;
		using TauLeaping::epsilon;
		using TauLeaping::threshold;
		using TauLeaping::trimed_list;
		//	using SSA::initialize;//for neg

#ifdef MATRIX_STOICHIOMETRY
		typedef StandardDriverTypes::stoichiometryRow matrixrow;
#endif

		//	boost::numeric::ublas::vector<double> previousReactionCounts;
		std::list<std::size_t> criticalSpecies;//for neg
		std::list<std::size_t> noncriticalSpecies;//for neg
		std::list<std::size_t> revisedSpecies;//for simplification
		std::list<std::vector<int> > revisedPair;//reversible pair
		std::list<std::size_t> mergedSpecies;//noncritical and critical species
		//	std::size_t criticalThreshold;//for neg
		//	double criticalPropensitySum;//for neg
		boost::numeric::ublas::vector<std::size_t> tagList;//for neg and simplification
		//	std::vector<std::vector<std::size_t> >speciesToReaction;//for neg
		//	std::vector<std::vector<std::size_t> >reactionToSpecies;//for neg
		//	boost::numeric::ublas::vector<std::size_t> affectedReactions;//for neg
		//	boost::numeric::ublas::vector<std::size_t> affectedSpecies;//for neg
		std::vector<std::vector<std::size_t> >relatedSpecies;//for simplification
		boost::numeric::ublas::vector<std::size_t> relatedSpeciesSize;//for simplification
		std::vector<std::vector<bool> > connectedSpecies;
		boost::numeric::ublas::vector<std::size_t> connectedSpeciesSize;

	//	std::vector<std::vector<std::size_t> >inputReactions;//for simplification
	//	std::vector<std::vector<std::size_t> >outputReactions;//for simplification
	//	std::vector<std::vector<std::size_t> >inputStoi;//for simplification, not used
	//	std::vector<std::vector<std::size_t> >outputStoi;//for simplification, not used
	//	boost::numeric::ublas::vector<bool>revisable;//for simplification
	//	boost::numeric::ublas::vector<double>stepsize;//for simplification
		//	boost::numeric::ublas::vector<double>last_stepsize;//stepsize for the last step
	//	std::deque<std::vector<double> >possible_size;//for simplification, not used

		std::size_t NumberOfRevisedSpecies;
		//	std::list<int> relaxedSpecies//for simplification
		boost::numeric::ublas::vector<std::size_t> stepsize_order;//for simplification

		inline static bool cmp(std::size_t i, std::size_t j, improvedTauLeaping_aggressive<_denseVectorType, _matrixType, _propensitiesFunctorType, _dependencyGraphType> *pt)
		{
			return pt->stepsize[i]<pt->stepsize[j];
		}

		//! the minimum number of reactions one tauleaping step should have
		//	std::size_t threshold;//set to zero to prevent switching to ssa
		//! the number of ssa steps taken when switching from tauleaping to ssa.
		//	std::size_t SSASteps;
		//! Epsilon used by tau-leaping to determine the tau
		//	double epsilon;

		/*! \brief the G vector described in "Efficient step size selection for the tau-leaping simulation method" */
		//_denseVectorType G;
		/* use instead of a vector--just set it to 2.0 which is conservative
		except if species can dimerize with itself and has small population */
		//	double g;

		//! squared elements stoichiometric matrix
		//	_matrixType squaredVj;
		//	_denseVectorType mu;
		//	_denseVectorType sigmaSquared;
		std::size_t tauleapingSteps;

		//reevaluate the stepsize
		std::vector<std::vector<std::size_t> >influencedSpecies;//for simplification
	//	boost::numeric::ublas::vector<std::size_t> influencedSpeciesSize;//for simplification

		double accumulator;
		boost::numeric::ublas::vector<double>reevaluatedSpecies;//reevaluated species are given value accumulator
		//	boost::numeric::ublas::vector<double>reevaluatedReaction;//reevaluated reactions are given value accumulator
	//	boost::numeric::ublas::vector<std::size_t>responsiblespeices;
	//	std::vector<int> reverseReaction;//not used
		void getReevaluated();//assign accumulator to the reevaluated species and reactions

		void getParameters();//calculate c, c_sum, a_sum;

		boost::numeric::ublas::matrix<double>stoi;
		boost::numeric::ublas::matrix<double>squarestoi;

		boost::numeric::ublas::matrix<int>revisable_2nd;

		inline double fmean(std::size_t i, double t);//i: the reaction index
		inline double dmean(std::size_t i, double t);
		inline double fvar(std::size_t i, double t);
		inline double dvar(std::size_t i, double t);
		inline double fcatalyst(std::size_t i, double t, double &dmean);//i: the reaction index
		inline double varcatalyst(std::size_t i, double t, double &dvar);//i: the reaction index

		void form_formular(std::size_t i); //creat the normalReaction list and revisedReaction list for species i
		double meanf(std::size_t i, double t, bool &ifpositive);
		double meandf(std::size_t i, double t, bool ifpositive);
		double varf(std::size_t i, double t);
		double vardf(std::size_t i, double t);
		double solve_tau(std::size_t i, double tau);
		double reevaluate_tau(double criticalStepsize, double endtime);
#ifdef PRINT_RELAXED
		double reevaluate_tau(double criticalStepsize, double endtime, int &confineSpecies);
#endif

		//for Jacobian
		inline double singleMolecule(int i, int j)//i:species; j:reactions.
		{
			return rateList[j];
		}

		inline double differentMolecule(int i, int j)//i:species; j:reactions
		{
			return rateList[j]*currentPopulation[i];
		}

		inline double sameMolecule(int i, int j)//i:species; j:reactions
		{
			return rateList[j]*(currentPopulation[i]-0.5);
		}

		std::vector<double (improvedTauLeaping_aggressive<_denseVectorType,_matrixType,_propensitiesFunctorType,_dependencyGraphType>::*)(int, int)> Jacobian;//the list of Jacobian
		inline double getJacobian(int i, int j);//i:species index

	private:
		//! default constructor not implemented
		improvedTauLeaping_aggressive();

		//for calculation speed
	//	std::vector<std::vector<int> > torevisedreaction;
		void update_time_invariant();
		void update_time_dependent(double tau);
#ifdef DEBUG
		std::vector<std::vector<std::size_t> >inputReactions;//for simplification
		std::vector<std::size_t> NumberOfInputReactions;//NumberOfInputReactions[i]=inputReactions[i].size()
		std::vector<std::vector<std::size_t> >outputReactions;//for simplification
		std::vector<std::size_t> NumberOfOutputReactions;//NumberOfOutputReactions[i]=onputReactions[i].size()

		boost::numeric::ublas::vector<std::size_t> influencedSpeciesSize;//for simplification
		boost::numeric::ublas::vector<bool>revisable;//for simplification
		boost::numeric::ublas::vector<double>stepsize;//for simplification
		boost::numeric::ublas::vector<std::size_t>responsiblespeices;
		std::vector<std::vector<int> >responsiblepair;
		boost::numeric::ublas::vector<double>c;//rate
		boost::numeric::ublas::vector<double>c_sum;
		boost::numeric::ublas::vector<double>a_sum;
		boost::numeric::ublas::vector<double>c_in;//c21 for s1
		boost::numeric::ublas::vector<double>c_out;//c12 for s1
		std::vector<std::vector<std::vector<int> > >revisable_2nd_reaction;
		std::vector<std::size_t> number_of_revisable_catalyze_reaction;//revisable_2nd_reaction[i][j].size()
		std::vector<std::vector<std::size_t> >number_of_revisable_2nd_reaction;
		std::vector<std::vector<int> > revisable_catalyze_reaction;
		std::vector<int> normalReaction;//used by the function form_formular and solve_tau
		std::vector<int> revisedReaction;//used by the function form_formular and solve_tau
		std::vector<int> pairProduct;//used by the function form_formular and solve_tau
		std::vector<int> pairReactant;//used by the function form_formular and solve_tau
		std::vector<int> oneDirectionReaction;//used by the function form_formular and solve_tau
		std::vector<int> catalyzedReaction;//used by the function form_formular and solve_tau
		std::vector<std::vector<int> > reactantList;
		std::vector<double> rateList;
		std::vector<double> c_over_csum;
		std::vector<double> asum_over_csum;
		std::vector<double> cx;
		std::vector<double> exp_ct;
		std::vector<double> one_exp_ct;
		std::vector<double> c2;
		std::vector<std::vector<PairedParameters> >pairedparameter;
#else
		std::size_t **inputReactions;//for simplification
		std::size_t *NumberOfInputReactions;//NumberOfInputReactions[i]=inputReactions[i].size()
		std::size_t **outputReactions;//for simplification
		std::size_t *NumberOfOutputReactions;//NumberOfOutputReactions[i]=onputReactions[i].size()
		
		std::size_t *influencedSpeciesSize;//for simplification
		bool *revisable;//for simplification
		double *stepsize;//for simplification
		std::size_t *responsiblespeices;
		int **responsiblepair;
		double *c;//rate
		double *c_sum;
		double *a_sum;
		double *c_in;//c21 for s1
		double *c_out;//c12 for s1
		int ***revisable_2nd_reaction;
		std::size_t **number_of_revisable_2nd_reaction;//revisable_2nd_reaction[i][j].size()
		int **revisable_catalyze_reaction;
		std::size_t *number_of_revisable_catalyze_reaction;
		int *normalReaction;//used by the function form_formular and solve_tau
		int *revisedReaction;//used by the function form_formular and solve_tau
		int *pairProduct;//used by the function form_formular and solve_tau
		int *pairReactant;//used by the function form_formular and solve_tau
		int *oneDirectionReaction;//used by the function form_formular and solve_tau
		int *catalyzedReaction;//used by the function form_formular and solve_tau

		int **reactantList;
		double *rateList;
		double *c_over_csum;
		double *asum_over_csum;
		double *cx;
		double *exp_ct;
		double *one_exp_ct;
		double *c2;
		PairedParameters **pairedparameter;
		void freeMemory();
#endif
		std::size_t NumberOfNormalReactions;
		std::size_t NumberOfRevisedReactions;
		std::size_t NumberOfPairProducts;
		std::size_t NumberOfPairReactants;
		std::size_t NumberOfOneDirectionReactions;
		std::size_t NumberOfCatalyzedReactions;
		double savedmean, savedvar;//save the dmean and dvar for paired reactions to save calculation

	public:

		//! Constructor
#ifdef JACOBIAN
		improvedTauLeaping_aggressive(const _denseVectorType& initialPop,
			const _matrixType& stoich,
			const _propensitiesFunctorType& propensitiesFunctor,
			const _dependencyGraphType& depGraph,
			const std::vector<std::vector<int> > &reactant,
			const std::vector<double> &rate,
			int seed=time(NULL));
#else
		improvedTauLeaping_aggressive(const _denseVectorType& initialPop,
			const _matrixType& stoich,
			const _propensitiesFunctorType& propensitiesFunctor,
			const _dependencyGraphType& depGraph,
			int seed=time(NULL));
#endif

		//! compiler-generated copy constructor OK
		//! compiler-generated assignment operator OK

		//! destructor
		virtual ~improvedTauLeaping_aggressive(){}

		void initialize(double startTime, double endTime);

		void creatPairMatrix();//creat revisable_2nd

		void findPair(int speciesId, int reactionId, std::vector<std::vector<std::vector<int> > > &Revisable_2nd_reaction);//find the reverble pair of speciesId in reactionId

		/*	void setThreshold(std::size_t threshold) {
		this->threshold=threshold;
		}*/

		/*	void setSSASteps(std::size_t ssaSteps) {
		SSASteps=ssaSteps;
		}*/

		void settauleapingSteps(std::size_t tauSteps) {
			tauleapingSteps=tauSteps;
		}

		//	void setEpsilon(double epsilon);

		void selectTau(double &noncriticalStepsize, double &criticalStepsize);

		//should this return (a reference to) the vector of reaction counts?
		int selectReactions(double leapSize, bool runCritical);

		//should this take (a reference to) the vector of reaction counts?
		//	bool fireReactions(int criticalIndex);

		//methods for critical reactions
		//	double critical_selectStepSize();
		//	int critical_selectReaction();
		//	bool critical_fireReaction(int reactionIndex);
		//	void critical_rollBack(int reactionIndex);
		void updateTagLists();

		/*!
		\brief run an ensemble simulation with output recorded at fixed time intervals

		output must have a conforming initialize, getOutputTimes, and record method
		outputTimes should be set in output prior to calling simulate
		if doValidate=true (the default) calls validate before ensemble
		calls initialize before each realization

		\param realizations number of simulations in the ensemble
		\param startTime the initial value of currentTime for each realization
		\param endTime the end time of each realization
		\param Output the class that handles storing the output for the simulation
		*/
		template<typename IntervalOutputType>
		void simulate(std::size_t realizations, double startTime, double endTime, IntervalOutputType& output, bool doValidate=true);

#ifdef PROFILE_SIMULATE_TAU_LEAPING_EXPLICIT_ADAPTIVE
		double numberOfIteration, numberOfEquation;
#endif

		int lastIterationSpecies;
		double lastIterationStepsize;
		double try_solve_tau(std::size_t i, double tau);

	};//end TauLeapingExplicitAdaptive class
}

#define _IMPROVED_TAU_LEAPING_AGGRESSIVE_IPP_
#include "improved_tau_leaping_aggressive.ipp"
#undef _IMPROVED_TAU_LEAPING_AGGRESSIVE_IPP_

#endif
