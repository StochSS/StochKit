#if !defined(__TestModelModel_h__)
#define __TestModelModel_h__

//model files will have a consistent name like "ProblemDefinition"
//notice NumberOfReactions and NumberOfSpecies appear several times
//create a parameters list that can be modified from driver (not preprocessor #DEFINEs)
//for example, VOLUME might show up in propensities and initial population, only want to define/change in one place
//events?

//_matrixType should be a dense vector of (dense or sparse) stoichiometry vectors
//the stoichiometry vectors should be compatible with the initial population's _populationVectorType via + operator
template<typename _matrixType>
_matrixType TestModelStoichiometry() {
  typedef typename _matrixType::value_type vectorType;
  const int NumberOfReactions=4;
  const int NumberOfSpecies=3;
  _matrixType nu(NumberOfReactions, vectorType(NumberOfSpecies));
  nu[0][0]=-1;//reaction 0: S0->0
  nu[1][0]=-2; nu[1][1]=1;//reaction 1: 2S0->S1
  nu[2][0]=2; nu[2][1]=-1;//reaction 2: S1->2S0
  nu[3][1]=-1; nu[3][2]=1;//reaction 3: S1->S2

  return nu;
};

//dependency graph is dense vector (length=NumberOfReactions) of variable-length vectors
//the ith variable-length vector stores the rxn indices of propensities affected by rxn i
//not as efficient as storing in contiguous memory?
//could make more efficient by specifying a capacity for the variable-length vectors?
template<typename _graphType>
_graphType TestModelDependencyGraph() {
  const int NumberOfReactions=4;
  _graphType dg(NumberOfReactions);
  dg[0].push_back(0); dg[0].push_back(1);//rxn 0 affects reactions 0 and 1 (compare to rxn 3)
  dg[1].push_back(0); dg[1].push_back(1); dg[1].push_back(2); dg[1].push_back(3);//rxn 1 affects all rxns (0,1,2 and 3)
  dg[2].push_back(0); dg[2].push_back(1); dg[2].push_back(2); dg[2].push_back(3);//rxn 2 affects all rxns (0,1,2 and 3)
  dg[3].push_back(2); dg[3].push_back(3);//rxn 3 affects rxns 2 and 3
  
  return dg;
};

template<typename _populationVectorType>
_populationVectorType TestModelInitialPopulations() {
  const int NumberOfSpecies=3;
  _populationVectorType X(NumberOfSpecies);//will initialize values to 0
  X[0]=10000;
  return X;
};

template<typename _populationVectorType>
class TestModelPropensities {
 public:
  static const int NumberOfReactions = 4;

  //! A pointer to a member function that computes a single propensity.
  typedef double (TestModelPropensities::* PropensityMember) (_populationVectorType&);

  PropensityMember _propensityFunctions[NumberOfReactions];

  //demonstrate how to make "global" (model-level) variables
  //note they are initialized in the constructor (when it calls init)
  double VOLUME;
  double global_c;

  void init() {
    VOLUME=1.0;
    global_c=1.0;
    _propensityFunctions[0] = &TestModelPropensities::f0;
    _propensityFunctions[1] = &TestModelPropensities::f1;
    _propensityFunctions[2] = &TestModelPropensities::f2;
    _propensityFunctions[3] = &TestModelPropensities::f3;
  }

  //! Assignment operator not implemented.
  TestModelPropensities& operator=(const TestModelPropensities&);

  //! Default constructor.
  TestModelPropensities() {
    init();
  }

  //! Copy constructor.
  TestModelPropensities
    (const TestModelPropensities& other) {
    init();
  }

  //! Destructor.
  ~TestModelPropensities() 
  {}

  //! Return the specified propensity function.
  double
  operator()(const int n, _populationVectorType& populations) {
    return (this->*_propensityFunctions[n])(populations);
  }


  //@}
  //--------------------------------------------------------------------------
  // Compute propensities.
private:

  double
  f0(_populationVectorType& x) {
    //show how to use global (model-level) parameters
    return VOLUME*global_c*x[0];
  }
  
  double
  f1(_populationVectorType& x) {
    //show how to use local (reaction-level) constants
    //these can't be changed by other code
    double local_c=.002;
    double local_scale_factor=2.0;
    return (local_c / local_scale_factor) * x[0] * (x[0] - 1);
  }
  
  double
  f2(_populationVectorType& x) {
    return 0.5 * x[1];
  }
  
  double
  f3(_populationVectorType& x) {
    return 0.04 * x[1];
  }
  
};

#endif
