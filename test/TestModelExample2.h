#if !defined(__TestModelModelExample2_h__)
#define __TestModelModelExample2_h__

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
  const int NumberOfReactions=3;
  const int NumberOfSpecies=2;
  _matrixType nu(NumberOfReactions, vectorType(NumberOfSpecies));
  nu[0][0]=-1;//reaction 0: S0->S2
  nu[1][0]=1;//reaction 1: S2->S0
  nu[2][1]=-1;//reaction 2: S0+S1->S0+S3

  return nu;
};

//dependency graph is dense vector (length=NumberOfReactions) of variable-length vectors
//the ith variable-length vector stores the rxn indices of propensities affected by rxn i
//not as efficient as storing in contiguous memory?
//could make more efficient by specifying a capacity for the variable-length vectors?
//Correctness????
template<typename _graphType>
_graphType TestModelDependencyGraph() {
  const int NumberOfReactions=3;
  _graphType dg(NumberOfReactions);
  dg[0].push_back(0); dg[0].push_back(1); dg[0].push_back(2);//rxn 0 affects reactions 0 and 1 (compare to rxn 3)
  dg[1].push_back(0); dg[1].push_back(1); dg[1].push_back(2);//rxn 1 affects all rxns (0,1,2 and 3)
  dg[2].push_back(2); //rxn 2 affects all rxns (0,1,2 and 3)
  
  return dg;
};

template<typename _populationVectorType>
_populationVectorType TestModelInitialPopulations() {
  const int NumberOfSpecies=2;
  _populationVectorType X(NumberOfSpecies);//will initialize values to 0
  X[0]=10000;
  X[1]=100;
  return X;
};

template<typename _populationVectorType>
class TestModelPropensities {
 public:
  static const int NumberOfReactions = 3;
  static const int NumberOfSpecies = 2;

  //! A pointer to a member function that computes a single propensity.
  typedef double (TestModelPropensities::* PropensityMember) (const _populationVectorType&) const;

  PropensityMember _propensityFunctions[NumberOfReactions];
  PropensityMember _propensityFunctionsJacobian[NumberOfReactions][NumberOfSpecies];

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
    _propensityFunctionsJacobian[0][0] = &TestModelPropensities::f0_ds0;
    _propensityFunctionsJacobian[0][1] = &TestModelPropensities::f0_ds1;
    _propensityFunctionsJacobian[1][0] = &TestModelPropensities::f1_ds0;
    _propensityFunctionsJacobian[1][1] = &TestModelPropensities::f1_ds1;
    _propensityFunctionsJacobian[2][0] = &TestModelPropensities::f2_ds0;
    _propensityFunctionsJacobian[2][1] = &TestModelPropensities::f2_ds1;
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
  operator()(const int n, const _populationVectorType& populations) const {
    return (this->*_propensityFunctions[n])(populations);
  }

  //! Return the specified propensity jacobian function.
  double
  operator()(const int n, const int s, const _populationVectorType& populations) const {
    return (this->*_propensityFunctionsJacobian[n][s])(populations);
  }

  //@}
  //--------------------------------------------------------------------------
  // Compute propensities.
private:

  double
  f0(const _populationVectorType& x) const {
    //show how to use global (model-level) parameters
    const double c1 = 100000.0;
    return c1*x[0];
  }
  
  double
  f1(const _populationVectorType& x) const {
    //show how to use local (reaction-level) constants
    //these can't be changed by other code
    const double c2 = 100000.0;
    const double xT = 20000.0;
    return c2 * (xT - x[0]);
  }
  
  double
  f2(const _populationVectorType& x) const {
     const double c3 = 0.0005;
     return c3 * x[0] * x[1];
  }
  
//compute Jacobians
  double
  f0_ds0(const _populationVectorType& x) const {
    //show how to use global (model-level) parameters
    const double c1 = 100000.0;
    return c1;
  }

  double
  f0_ds1(const _populationVectorType& x) const {
    return 0.0;
  }

  double
  f1_ds0(const _populationVectorType& x) const {
    //show how to use local (reaction-level) constants
    //these can't be changed by other code
    const double c2 = 100000.0;
    //const double xT = 20000.0;
    return -c2;
  }

  double
  f1_ds1(const _populationVectorType& x) const {
    return 0.0;
  }

  double
  f2_ds0(const _populationVectorType& x) const {
     const double c3 = 0.0005;
     return c3 * x[1];
  }

  double
  f2_ds1(const _populationVectorType& x) const {
     const double c3 = 0.0005;
     return c3 * x[0];
  }
};

#endif
