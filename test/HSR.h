#if !defined(__TestModelModel_h__)
#define __TestModelModel_h__


//model files will have a consistent name like "ProblemDefinition"
//notice NumberOfReactions and NumberOfSpecies appear several times
//create a parameters list that can be modified from driver (not preprocessor #DEFINEs)
//for example, VOLUME might show up in propensities and initial population, only want to define/change in one place
//events?

//_matrixType should be a dense vector of [dense or sparse) stoichiometry vectors
//the stoichiometry vectors should be compatible with the initial population's _populationVectorType via + operator
template<typename _matrixType>
_matrixType TestModelStoichiometry() {
  typedef typename _matrixType::value_type vectorType;
  const int NumberOfReactions=61;
  const int NumberOfSpecies=28;
  _matrixType nu(NumberOfReactions, vectorType(NumberOfSpecies));
/*  nu[0][0]=-1;//reaction 0: S0->0
  nu[1][0]=-2; nu[1][1]=1;//reaction 1: 2S0->S1
  nu[2][0]=2; nu[2][1]=-1;//reaction 2: S1->2S0
  nu[3][1]=-1; nu[3][2]=1;//reaction 3: S1->S2
*/
   nu[0][0] = -1;
   nu[0][1] = -1;
   nu[0][2] = 1;

   nu[1][0] = 1;
   nu[1][1] = 1;
   nu[1][2] = -1;

   nu[2][0] = -1;
   nu[2][3] = -1;
   nu[2][4] = 1;

   nu[3][0] = 1;
   nu[3][3] = 1;
   nu[3][4] = -1;

   nu[4][0] = -1;
   nu[4][5] = -1;
   nu[4][6] = 1;

   nu[5][0] = 1;
   nu[5][5] = 1;
   nu[5][6] = -1;

   nu[6][3] = -1;
   nu[6][13] = -1;
   nu[6][14] = 1;

   nu[7][3] = 1;
   nu[7][13] = 1;
   nu[7][14] = -1;

   nu[8][13] = -1;
   nu[8][15] = -1;
   nu[8][16] = 1;

   nu[9][13] = 1;
   nu[9][15] = 1;
   nu[9][16] = -1;

   nu[10][2] = -1;
   nu[10][5] = -1;
   nu[10][7] = 1;

   nu[11][2] = 1;
   nu[11][5] = 1;
   nu[11][7] = -1;

   nu[12][4] = -1;
   nu[12][5] = -1;
   nu[12][8] = 1;

   nu[13][4] = 1;
   nu[13][5] = 1;
   nu[13][8] = -1;

   nu[14][2] = -1;
   nu[14][9] = -1;
   nu[14][11] = 1;

   nu[15][2] = 1;
   nu[15][9] = 1;
   nu[15][11] = -1;

   nu[16][4] = -1;
   nu[16][10] = -1;
   nu[16][12] = 1;

   nu[17][4] = 1;
   nu[17][10] = 1;
   nu[17][12] = -1;

   nu[18][14] = -1;
   nu[18][17] = -1;
   nu[18][18] = 1;

   nu[19][14] = 1;
   nu[19][17] = 1;
   nu[19][18] = -1;

   nu[20][21] = 1;

   nu[21][21] = -1;

   nu[22][13] = 1;

   nu[23][13] = -1;

   nu[24][15] = 1;
   nu[24][16] = -1;

   nu[25][3] = 1;
   nu[25][14] = -1;

   nu[26][3] = 1;
   nu[26][17] = 1;
   nu[26][18] = -1;

   nu[27][3] = 1;
   nu[27][26] = 1;
   nu[27][27] = -1;

   nu[28][22] = 1;

   nu[29][22] = -1;

   nu[30][17] = 1;

   nu[31][17] = -1;

   nu[32][14] = 1;
   nu[32][18] = -1;

   nu[33][24] = 1;

   nu[34][24] = -1;

   nu[35][3] = 1;

   nu[36][3] = -1;

   nu[37][13] = 1;
   nu[37][17] = 1;
   nu[37][18] = -1;

   nu[38][19] = 1;
   nu[38][20] = -1;

   nu[39][13] = 1;
   nu[39][26] = 1;
   nu[39][27] = -1;

   nu[40][23] = 1;

   nu[41][23] = -1;

   nu[42][19] = 1;

   nu[43][19] = -1;

   nu[44][3] = 1;
   nu[44][20] = -1;

   nu[45][3] = -1;
   nu[45][19] = -1;
   nu[45][20] = 1;

   nu[46][3] = 1;
   nu[46][19] = 1;
   nu[46][20] = -1;

   nu[47][25] = 1;

   nu[48][25] = -1;

   nu[49][26] = 1;

   nu[50][26] = -1;

   nu[51][14] = 1;
   nu[51][27] = -1;

   nu[52][14] = -1;
   nu[52][26] = -1;
   nu[52][27] = 1;

   nu[53][14] = 1;
   nu[53][26] = 1;
   nu[53][27] = -1;

   nu[54][0] = 1;
   nu[54][4] = -1;

   nu[55][0] = 1;
   nu[55][10] = 1;
   nu[55][12] = -1;

   nu[56][6] = 1;
   nu[56][8] = -1;

   nu[57][13] = 1;
   nu[57][14] = -1;

   nu[58][13] = 1;
   nu[58][17] = 1;
   nu[58][18] = -1;

   nu[59][14] = 1;
   nu[59][26] = 1;
   nu[59][27] = -1;

   nu[60][19] = 1;
   nu[60][20] = -1;
  
  return nu;
};

//dependency graph is dense vector (length=NumberOfReactions) of variable-length vectors
//the ith variable-length vector stores the rxn indices of propensities affected by rxn i
//not as efficient as storing in contiguous memory?
//could make more efficient by specifying a capacity for the variable-length vectors?
template<typename _graphType>
_graphType TestModelDependencyGraph() {
  const int NumberOfReactions=61;
  _graphType dg(NumberOfReactions);
//  dg[0].push_back(0); dg[0].push_back(1);//rxn 0 affects reactions 0 and 1 (compare to rxn 3)
//  dg[1].push_back(0); dg[1].push_back(1); dg[1].push_back(2); dg[1].push_back(3);//rxn 1 affects all rxns (0,1,2 and 3)
//  dg[2].push_back(0); dg[2].push_back(1); dg[2].push_back(2); dg[2].push_back(3);//rxn 2 affects all rxns (0,1,2 and 3)
//  dg[3].push_back(2); dg[3].push_back(3);//rxn 3 affects rxns 2 and 3
  dg[0].push_back(0);   dg[0].push_back(1);   dg[0].push_back(2);   dg[0].push_back(4);   dg[0].push_back(10);   dg[0].push_back(14); 
  dg[1].push_back(0);   dg[1].push_back(1);   dg[1].push_back(2);   dg[1].push_back(4);   dg[1].push_back(10);   dg[1].push_back(14); 
  dg[2].push_back(0);   dg[2].push_back(2);   dg[2].push_back(3);   dg[2].push_back(4);   dg[2].push_back(6);   dg[2].push_back(12);   dg[2].push_back(16);   dg[2].push_back(36);   dg[2].push_back(45);   dg[2].push_back(54); 
  dg[3].push_back(0);   dg[3].push_back(2);   dg[3].push_back(3);   dg[3].push_back(4);   dg[3].push_back(6);   dg[3].push_back(12);   dg[3].push_back(16);   dg[3].push_back(36);   dg[3].push_back(45);   dg[3].push_back(54); 
  dg[4].push_back(0);   dg[4].push_back(2);   dg[4].push_back(4);   dg[4].push_back(5);   dg[4].push_back(10);   dg[4].push_back(12); 
  dg[5].push_back(0);   dg[5].push_back(2);   dg[5].push_back(4);   dg[5].push_back(5);   dg[5].push_back(10);   dg[5].push_back(12); 
  dg[6].push_back(2);   dg[6].push_back(6);   dg[6].push_back(7);   dg[6].push_back(8);   dg[6].push_back(18);   dg[6].push_back(23);   dg[6].push_back(25);   dg[6].push_back(36);   dg[6].push_back(45);   dg[6].push_back(52);   dg[6].push_back(57); 
  dg[7].push_back(2);   dg[7].push_back(6);   dg[7].push_back(7);   dg[7].push_back(8);   dg[7].push_back(18);   dg[7].push_back(23);   dg[7].push_back(25);   dg[7].push_back(36);   dg[7].push_back(45);   dg[7].push_back(52);   dg[7].push_back(57); 
  dg[8].push_back(6);   dg[8].push_back(8);   dg[8].push_back(9);   dg[8].push_back(23);   dg[8].push_back(24); 
  dg[9].push_back(6);   dg[9].push_back(8);   dg[9].push_back(9);   dg[9].push_back(23);   dg[9].push_back(24); 
  dg[10].push_back(1);   dg[10].push_back(4);   dg[10].push_back(10);   dg[10].push_back(11);   dg[10].push_back(12);   dg[10].push_back(14); 
  dg[11].push_back(1);   dg[11].push_back(4);   dg[11].push_back(10);   dg[11].push_back(11);   dg[11].push_back(12);   dg[11].push_back(14); 
  dg[12].push_back(3);   dg[12].push_back(4);   dg[12].push_back(10);   dg[12].push_back(12);   dg[12].push_back(13);   dg[12].push_back(16);   dg[12].push_back(54);   dg[12].push_back(56); 
  dg[13].push_back(3);   dg[13].push_back(4);   dg[13].push_back(10);   dg[13].push_back(12);   dg[13].push_back(13);   dg[13].push_back(16);   dg[13].push_back(54);   dg[13].push_back(56); 
  dg[14].push_back(1);   dg[14].push_back(10);   dg[14].push_back(14);   dg[14].push_back(15);   dg[14].push_back(33); 
  dg[15].push_back(1);   dg[15].push_back(10);   dg[15].push_back(14);   dg[15].push_back(15);   dg[15].push_back(33); 
  dg[16].push_back(3);   dg[16].push_back(12);   dg[16].push_back(16);   dg[16].push_back(17);   dg[16].push_back(20);   dg[16].push_back(28);   dg[16].push_back(40);   dg[16].push_back(47);   dg[16].push_back(54);   dg[16].push_back(55); 
  dg[17].push_back(3);   dg[17].push_back(12);   dg[17].push_back(16);   dg[17].push_back(17);   dg[17].push_back(20);   dg[17].push_back(28);   dg[17].push_back(40);   dg[17].push_back(47);   dg[17].push_back(54);   dg[17].push_back(55); 
  dg[18].push_back(7);   dg[18].push_back(18);   dg[18].push_back(19);   dg[18].push_back(25);   dg[18].push_back(26);   dg[18].push_back(31);   dg[18].push_back(32);   dg[18].push_back(37);   dg[18].push_back(52);   dg[18].push_back(57);   dg[18].push_back(58); 
  dg[19].push_back(7);   dg[19].push_back(18);   dg[19].push_back(19);   dg[19].push_back(25);   dg[19].push_back(26);   dg[19].push_back(31);   dg[19].push_back(32);   dg[19].push_back(37);   dg[19].push_back(52);   dg[19].push_back(57);   dg[19].push_back(58); 
  dg[20].push_back(17);   dg[20].push_back(20);   dg[20].push_back(21);   dg[20].push_back(22);   dg[20].push_back(28);   dg[20].push_back(40);   dg[20].push_back(47);   dg[20].push_back(55); 
  dg[21].push_back(21);   dg[21].push_back(22); 
  dg[22].push_back(6);   dg[22].push_back(8);   dg[22].push_back(21);   dg[22].push_back(22);   dg[22].push_back(23); 
  dg[23].push_back(6);   dg[23].push_back(8);   dg[23].push_back(23); 
  dg[24].push_back(8);   dg[24].push_back(9);   dg[24].push_back(24); 
  dg[25].push_back(2);   dg[25].push_back(6);   dg[25].push_back(7);   dg[25].push_back(18);   dg[25].push_back(25);   dg[25].push_back(36);   dg[25].push_back(45);   dg[25].push_back(52);   dg[25].push_back(57); 
  dg[26].push_back(2);   dg[26].push_back(6);   dg[26].push_back(18);   dg[26].push_back(19);   dg[26].push_back(26);   dg[26].push_back(31);   dg[26].push_back(32);   dg[26].push_back(36);   dg[26].push_back(37);   dg[26].push_back(45);   dg[26].push_back(58); 
  dg[27].push_back(2);   dg[27].push_back(6);   dg[27].push_back(27);   dg[27].push_back(36);   dg[27].push_back(39);   dg[27].push_back(45);   dg[27].push_back(50);   dg[27].push_back(51);   dg[27].push_back(52);   dg[27].push_back(53);   dg[27].push_back(59); 
  dg[28].push_back(17);   dg[28].push_back(20);   dg[28].push_back(28);   dg[28].push_back(29);   dg[28].push_back(30);   dg[28].push_back(40);   dg[28].push_back(47);   dg[28].push_back(55); 
  dg[29].push_back(29);   dg[29].push_back(30); 
  dg[30].push_back(18);   dg[30].push_back(29);   dg[30].push_back(30);   dg[30].push_back(31); 
  dg[31].push_back(18);   dg[31].push_back(31); 
  dg[32].push_back(7);   dg[32].push_back(18);   dg[32].push_back(19);   dg[32].push_back(25);   dg[32].push_back(26);   dg[32].push_back(32);   dg[32].push_back(37);   dg[32].push_back(52);   dg[32].push_back(57);   dg[32].push_back(58); 
  dg[33].push_back(15);   dg[33].push_back(33);   dg[33].push_back(34);   dg[33].push_back(35); 
  dg[34].push_back(34);   dg[34].push_back(35); 
  dg[35].push_back(2);   dg[35].push_back(6);   dg[35].push_back(34);   dg[35].push_back(35);   dg[35].push_back(36);   dg[35].push_back(45); 
  dg[36].push_back(2);   dg[36].push_back(6);   dg[36].push_back(36);   dg[36].push_back(45); 
  dg[37].push_back(6);   dg[37].push_back(8);   dg[37].push_back(18);   dg[37].push_back(19);   dg[37].push_back(23);   dg[37].push_back(26);   dg[37].push_back(31);   dg[37].push_back(32);   dg[37].push_back(37);   dg[37].push_back(58); 
  dg[38].push_back(38);   dg[38].push_back(43);   dg[38].push_back(44);   dg[38].push_back(45);   dg[38].push_back(46);   dg[38].push_back(60); 
  dg[39].push_back(6);   dg[39].push_back(8);   dg[39].push_back(23);   dg[39].push_back(27);   dg[39].push_back(39);   dg[39].push_back(50);   dg[39].push_back(51);   dg[39].push_back(52);   dg[39].push_back(53);   dg[39].push_back(59); 
  dg[40].push_back(17);   dg[40].push_back(20);   dg[40].push_back(28);   dg[40].push_back(40);   dg[40].push_back(41);   dg[40].push_back(42);   dg[40].push_back(47);   dg[40].push_back(55); 
  dg[41].push_back(41);   dg[41].push_back(42); 
  dg[42].push_back(41);   dg[42].push_back(42);   dg[42].push_back(43);   dg[42].push_back(45); 
  dg[43].push_back(43);   dg[43].push_back(45); 
  dg[44].push_back(2);   dg[44].push_back(6);   dg[44].push_back(36);   dg[44].push_back(38);   dg[44].push_back(44);   dg[44].push_back(45);   dg[44].push_back(46);   dg[44].push_back(60); 
  dg[45].push_back(2);   dg[45].push_back(6);   dg[45].push_back(36);   dg[45].push_back(38);   dg[45].push_back(43);   dg[45].push_back(44);   dg[45].push_back(45);   dg[45].push_back(46);   dg[45].push_back(60); 
  dg[46].push_back(2);   dg[46].push_back(6);   dg[46].push_back(36);   dg[46].push_back(38);   dg[46].push_back(43);   dg[46].push_back(44);   dg[46].push_back(45);   dg[46].push_back(46);   dg[46].push_back(60); 
  dg[47].push_back(17);   dg[47].push_back(20);   dg[47].push_back(28);   dg[47].push_back(40);   dg[47].push_back(47);   dg[47].push_back(48);   dg[47].push_back(49);   dg[47].push_back(55); 
  dg[48].push_back(48);   dg[48].push_back(49); 
  dg[49].push_back(48);   dg[49].push_back(49);   dg[49].push_back(50);   dg[49].push_back(52); 
  dg[50].push_back(50);   dg[50].push_back(52); 
  dg[51].push_back(7);   dg[51].push_back(18);   dg[51].push_back(25);   dg[51].push_back(27);   dg[51].push_back(39);   dg[51].push_back(51);   dg[51].push_back(52);   dg[51].push_back(53);   dg[51].push_back(57);   dg[51].push_back(59); 
  dg[52].push_back(7);   dg[52].push_back(18);   dg[52].push_back(25);   dg[52].push_back(27);   dg[52].push_back(39);   dg[52].push_back(50);   dg[52].push_back(51);   dg[52].push_back(52);   dg[52].push_back(53);   dg[52].push_back(57);   dg[52].push_back(59); 
  dg[53].push_back(7);   dg[53].push_back(18);   dg[53].push_back(25);   dg[53].push_back(27);   dg[53].push_back(39);   dg[53].push_back(50);   dg[53].push_back(51);   dg[53].push_back(52);   dg[53].push_back(53);   dg[53].push_back(57);   dg[53].push_back(59); 
  dg[54].push_back(0);   dg[54].push_back(2);   dg[54].push_back(3);   dg[54].push_back(4);   dg[54].push_back(12);   dg[54].push_back(16);   dg[54].push_back(54); 
  dg[55].push_back(0);   dg[55].push_back(2);   dg[55].push_back(4);   dg[55].push_back(16);   dg[55].push_back(17);   dg[55].push_back(20);   dg[55].push_back(28);   dg[55].push_back(40);   dg[55].push_back(47);   dg[55].push_back(55); 
  dg[56].push_back(5);   dg[56].push_back(13);   dg[56].push_back(56); 
  dg[57].push_back(6);   dg[57].push_back(7);   dg[57].push_back(8);   dg[57].push_back(18);   dg[57].push_back(23);   dg[57].push_back(25);   dg[57].push_back(52);   dg[57].push_back(57); 
  dg[58].push_back(6);   dg[58].push_back(8);   dg[58].push_back(18);   dg[58].push_back(19);   dg[58].push_back(23);   dg[58].push_back(26);   dg[58].push_back(31);   dg[58].push_back(32);   dg[58].push_back(37);   dg[58].push_back(58); 
  dg[59].push_back(6);   dg[59].push_back(8);   dg[59].push_back(23);   dg[59].push_back(27);   dg[59].push_back(39);   dg[59].push_back(50);   dg[59].push_back(51);   dg[59].push_back(52);   dg[59].push_back(53);   dg[59].push_back(59); 
  dg[60].push_back(38);   dg[60].push_back(43);   dg[60].push_back(44);   dg[60].push_back(45);   dg[60].push_back(46);   dg[60].push_back(60);   
  return dg;
};

template<typename _populationVectorType>
_populationVectorType TestModelInitialPopulations() {
  const int NumberOfSpecies=28;
  _populationVectorType X(NumberOfSpecies);//will initialize values to 0

   X[4] = 1;
   X[5] = 4645669;
   X[6] = 1324;
   X[7] = 80;
   X[8] = 16;
   X[9] = 3413;
   X[10] = 29;
   X[11] = 584;
   X[12] = 1;
   X[13] = 22;
   X[15] = 171440;
   X[16] = 9150;
   X[17] = 2280;
   X[18] = 6;
   X[19] = 596;
   X[20] = 0;
   X[21] = 13;
   X[22] = 3;
   X[23] = 3;
   X[24] = 7;
   X[26] = 260;

  return X;
};

template<typename _populationVectorType>
class TestModelPropensities {
 public:
  static const int NumberOfReactions = 61;

  //! A pointer to a member function that computes a single propensity.
  typedef double (TestModelPropensities::* PropensityMember) (_populationVectorType&);

  std::vector<PropensityMember> _propensityFunctions;
//  PropensityMember _propensityFunctions[NumberOfReactions];

  //demonstrate how to make "global" (model-level) variables
  //note they are initialized in the constructor (when it calls init)
  double VOLUME;
  double global_c;

  void init() {
    VOLUME=1.0;
    global_c=1.0;
	_propensityFunctions.resize(NumberOfReactions);
    _propensityFunctions[0] = &TestModelPropensities::f0;
    _propensityFunctions[1] = &TestModelPropensities::f1;
    _propensityFunctions[2] = &TestModelPropensities::f2;
    _propensityFunctions[3] = &TestModelPropensities::f3;
    _propensityFunctions[4] = &TestModelPropensities::f4;
    _propensityFunctions[5] = &TestModelPropensities::f5;
    _propensityFunctions[6] = &TestModelPropensities::f6;
    _propensityFunctions[7] = &TestModelPropensities::f7;
    _propensityFunctions[8] = &TestModelPropensities::f8;
    _propensityFunctions[9] = &TestModelPropensities::f9;
    _propensityFunctions[10] = &TestModelPropensities::f10;
    _propensityFunctions[11] = &TestModelPropensities::f11;
    _propensityFunctions[12] = &TestModelPropensities::f12;
    _propensityFunctions[13] = &TestModelPropensities::f13;
    _propensityFunctions[14] = &TestModelPropensities::f14;
    _propensityFunctions[15] = &TestModelPropensities::f15;
    _propensityFunctions[16] = &TestModelPropensities::f16;
    _propensityFunctions[17] = &TestModelPropensities::f17;
    _propensityFunctions[18] = &TestModelPropensities::f18;
    _propensityFunctions[19] = &TestModelPropensities::f19;
    _propensityFunctions[20] = &TestModelPropensities::f20;
    _propensityFunctions[21] = &TestModelPropensities::f21;
    _propensityFunctions[22] = &TestModelPropensities::f22;
    _propensityFunctions[23] = &TestModelPropensities::f23;
    _propensityFunctions[24] = &TestModelPropensities::f24;
    _propensityFunctions[25] = &TestModelPropensities::f25;
    _propensityFunctions[26] = &TestModelPropensities::f26;
    _propensityFunctions[27] = &TestModelPropensities::f27;
    _propensityFunctions[28] = &TestModelPropensities::f28;
    _propensityFunctions[29] = &TestModelPropensities::f29;
    _propensityFunctions[30] = &TestModelPropensities::f30;
    _propensityFunctions[31] = &TestModelPropensities::f31;
    _propensityFunctions[32] = &TestModelPropensities::f32;
    _propensityFunctions[33] = &TestModelPropensities::f33;
    _propensityFunctions[34] = &TestModelPropensities::f34;
    _propensityFunctions[35] = &TestModelPropensities::f35;
    _propensityFunctions[36] = &TestModelPropensities::f36;
    _propensityFunctions[37] = &TestModelPropensities::f37;
    _propensityFunctions[38] = &TestModelPropensities::f38;
    _propensityFunctions[39] = &TestModelPropensities::f39;
    _propensityFunctions[40] = &TestModelPropensities::f40;
    _propensityFunctions[41] = &TestModelPropensities::f41;
    _propensityFunctions[42] = &TestModelPropensities::f42;
    _propensityFunctions[43] = &TestModelPropensities::f43;
    _propensityFunctions[44] = &TestModelPropensities::f44;
    _propensityFunctions[45] = &TestModelPropensities::f45;
    _propensityFunctions[46] = &TestModelPropensities::f46;
    _propensityFunctions[47] = &TestModelPropensities::f47;
    _propensityFunctions[48] = &TestModelPropensities::f48;
    _propensityFunctions[49] = &TestModelPropensities::f49;
    _propensityFunctions[50] = &TestModelPropensities::f50;
    _propensityFunctions[51] = &TestModelPropensities::f51;
    _propensityFunctions[52] = &TestModelPropensities::f52;
    _propensityFunctions[53] = &TestModelPropensities::f53;
    _propensityFunctions[54] = &TestModelPropensities::f54;
    _propensityFunctions[55] = &TestModelPropensities::f55;
    _propensityFunctions[56] = &TestModelPropensities::f56;
    _propensityFunctions[57] = &TestModelPropensities::f57;
    _propensityFunctions[58] = &TestModelPropensities::f58;
    _propensityFunctions[59] = &TestModelPropensities::f59;
    _propensityFunctions[60] = &TestModelPropensities::f60;
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
   double C1=2.54;
   return C1*x[0]*x[1];
   }

   double
   f1(_populationVectorType& x) { 
   double C2=1;
   return C2*x[2];
   }

   double
   f2(_populationVectorType& x) { 
   double C3=0.254;
   return C3*x[0]*x[3];
   }

   double
   f3(_populationVectorType& x) { 
   double C4=1;
   return C4*x[4];
   }

   double
   f4(_populationVectorType& x) { 
   double C5=0.0254;
   return C5*x[0]*x[5];
   }

   double
   f5(_populationVectorType& x) { 
   double C6=10;
   return C6*x[6];
   }

   double
   f6(_populationVectorType& x) { 
   double C7=254;
   return C7*x[3]*x[13];
   }

   double
   f7(_populationVectorType& x) { 
   double C8=10000;
   return C8*x[14];
   }

   double
   f8(_populationVectorType& x) { 
   double C9=0.000254;
   return C9*x[13]*x[15];
   }

   double
   f9(_populationVectorType& x) { 
   double C10=0.01;
   return C10*x[16];
   }

   double
   f10(_populationVectorType& x) { 
   double C11=0.000254;
   return C11*x[2]*x[5];
   }

   double
   f11(_populationVectorType& x) { 
   double C12=1;
   return C12*x[7];
   }

   double
   f12(_populationVectorType& x) { 
   double C13=0.000254;
   return C13*x[4]*x[5];
   }

   double
   f13(_populationVectorType& x) { 
   double C14=1;
   return C14*x[8];
   }

   double
   f14(_populationVectorType& x) { 
   double C15=2.54;
   return C15*x[2]*x[9];
   }

   double
   f15(_populationVectorType& x) { 
   double C16=1;
   return C16*x[11];
   }

   double
   f16(_populationVectorType& x) { 
   double C17=2540;
   return C17*x[4]*x[10];
   }

   double
   f17(_populationVectorType& x) { 
   double C18=1000;
   return C18*x[12];
   }

   double
   f18(_populationVectorType& x) { 
   double C19=0.0254;
   return C19*x[14]*x[17];
   }

   double
   f19(_populationVectorType& x) { 
   double C20=1;
   return C20*x[18];
   }

   double
   f20(_populationVectorType& x) { 
   double C21=6.62;
   return C21*x[12];
   }

   double
   f21(_populationVectorType& x) { 
   double C22=0.5;
   return C22*x[21];
   }

   double
   f22(_populationVectorType& x) { 
   double C23=20;
   return C23*x[21];
   }

   double
   f23(_populationVectorType& x) { 
   double C24=0.03;
   return C24*x[13];
   }

   double
   f24(_populationVectorType& x) { 
   double C25=0.03;
   return C25*x[16];
   }

   double
   f25(_populationVectorType& x) { 
   double C26=0.03;
   return C26*x[14];
   }

   double
   f26(_populationVectorType& x) { 
   double C27=0.03;
   return C27*x[18];
   }

   double
   f27(_populationVectorType& x) { 
   double C28=0.03;
   return C28*x[27];
   }

   double
   f28(_populationVectorType& x) { 
   double C29=1.67;
   return C29*x[12];
   }

   double
   f29(_populationVectorType& x) { 
   double C30=0.5;
   return C30*x[22];
   }

   double
   f30(_populationVectorType& x) { 
   double C31=20;
   return C31*x[22];
   }

   double
   f31(_populationVectorType& x) { 
   double C32=0.03;
   return C32*x[17];
   }

   double
   f32(_populationVectorType& x) { 
   double C33=0.03;
   return C33*x[18];
   }

   double
   f33(_populationVectorType& x) { 
   double C34=0.00625;
   return C34*x[11];
   }

   double
   f34(_populationVectorType& x) { 
   double C35=0.5;
   return C35*x[24];
   }

   double
   f35(_populationVectorType& x) { 
   double C36=7;
   return C36*x[24];
   }

   double
   f36(_populationVectorType& x) { 
   double C37=0.03;
   return C37*x[3];
   }

   double
   f37(_populationVectorType& x) { 
   double C38=3;
   return C38*x[18];
   }

   double
   f38(_populationVectorType& x) { 
   double C39=0.7;
   return C39*x[20];
   }

   double
   f39(_populationVectorType& x) { 
   double C40=0.5;
   return C40*x[27];
   }

   double
   f40(_populationVectorType& x) { 
   double C41=1;
   return C41*x[12];
   }

   double
   f41(_populationVectorType& x) { 
   double C42=0.5;
   return C42*x[23];
   }

   double
   f42(_populationVectorType& x) { 
   double C43=20;
   return C43*x[23];
   }

   double
   f43(_populationVectorType& x) { 
   double C44=0.03;
   return C44*x[19];
   }

   double
   f44(_populationVectorType& x) { 
   double C45=0.03;
   return C45*x[20];
   }

   double
   f45(_populationVectorType& x) { 
   double C46=2.54;
   return C46*x[3]*x[19];
   }

   double
   f46(_populationVectorType& x) { 
   double C47=10000;
   return C47*x[20];
   }

   double
   f47(_populationVectorType& x) { 
   double C48=0.43333;
   return C48*x[12];
   }

   double
   f48(_populationVectorType& x) { 
   double C49=0.5;
   return C49*x[25];
   }

   double
   f49(_populationVectorType& x) { 
   double C50=20;
   return C50*x[25];
   }

   double
   f50(_populationVectorType& x) { 
   double C51=0.03;
   return C51*x[26];
   }

   double
   f51(_populationVectorType& x) { 
   double C52=0.03;
   return C52*x[27];
   }

   double
   f52(_populationVectorType& x) { 
   double C53=2.54;
   return C53*x[14]*x[26];
   }

   double
   f53(_populationVectorType& x) { 
   double C54=10000;
   return C54*x[27];
   }

   double
   f54(_populationVectorType& x) { 
   double C55=0.03;
   return C55*x[4];
   }

   double
   f55(_populationVectorType& x) { 
   double C56=0.03;
   return C56*x[12];
   }

   double
   f56(_populationVectorType& x) { 
   double C57=0.03;
   return C57*x[8];
   }

   double
   f57(_populationVectorType& x) { 
   double C58=0.03;
   return C58*x[14];
   }

   double
   f58(_populationVectorType& x) { 
   double C59=0.03;
   return C59*x[18];
   }

   double
   f59(_populationVectorType& x) { 
   double C60=0.03;
   return C60*x[27];
   }

   double
   f60(_populationVectorType& x) { 
   double C61=0.03;
   return C61*x[20];
   }

};
#endif
