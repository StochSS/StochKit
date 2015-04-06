/*!

*/

#if !defined(_IMPROVED_TAU_LEAPING_AGGRESSIVE_IPP_)
#error This file is the implementation of improved_TauLeaping_Aggressive
#endif
namespace STOCHKIT
{

#define templateList template<typename _denseVectorType, typename _matrixType, typename _propensitiesFunctorType, typename _dependencyGraphType>
#define	typeList improvedTauLeaping_aggressive<_denseVectorType, _matrixType, _propensitiesFunctorType, _dependencyGraphType>

	templateList typeList::improvedTauLeaping_aggressive
		(const _denseVectorType& initialPop,
		const _matrixType& stoich,
		const _propensitiesFunctorType& propensitiesFunctor,
		const _dependencyGraphType& depGraph,
		const std::vector<std::vector<int> >& reactant,
		const std::vector<double> &rate,
		int seed) : 
	TauLeapingExplicitAdaptive<_denseVectorType, _matrixType, _propensitiesFunctorType, _dependencyGraphType>(initialPop, stoich, propensitiesFunctor, depGraph, seed),
		tauleapingSteps(100),
		tagList(NumberOfReactions),
		relatedSpecies(NumberOfSpecies),
	//	inputReactions(NumberOfSpecies),
	//	outputReactions(NumberOfSpecies),
	//	inputStoi(NumberOfSpecies),
	//	outputStoi(NumberOfSpecies),
		influencedSpecies(NumberOfSpecies),
		connectedSpecies(NumberOfSpecies, std::vector<bool> (NumberOfSpecies, false) ),
	//	stepsize(NumberOfSpecies),
		stepsize_order(NumberOfSpecies),
	//	revisable(NumberOfSpecies, true),
		relatedSpeciesSize(NumberOfSpecies),
		reevaluatedSpecies(NumberOfSpecies),
		connectedSpeciesSize(NumberOfSpecies),
	//	responsiblespeices(NumberOfReactions),
	//	reverseReaction(NumberOfReactions, -1),
	//	responsiblepair(NumberOfReactions, std::vector<int> (2) ),
	//	influencedSpeciesSize(NumberOfSpecies),
		stoi(NumberOfReactions, NumberOfSpecies),
		revisable_2nd(NumberOfSpecies, NumberOfSpecies, -1),
		squarestoi(NumberOfReactions, NumberOfSpecies),
		Jacobian(NumberOfReactions)
	{
		std::size_t i, j, k;
		std::vector<std::vector<int> > Revisable_catalyze_reaction(NumberOfSpecies);
		std::vector<std::size_t> InputReactions;
		std::vector<std::size_t> OutputReactions;
		myNew(inputReactions, NumberOfSpecies);
		myNew(NumberOfInputReactions, NumberOfSpecies);
		myNew(outputReactions, NumberOfSpecies);
		myNew(NumberOfOutputReactions, NumberOfSpecies);

		myNew(influencedSpeciesSize, NumberOfSpecies);
		myNew(revisable, NumberOfSpecies, true),
		myNew(stepsize, NumberOfSpecies);
		myNew(responsiblespeices, NumberOfReactions);
		myNew(responsiblepair, NumberOfReactions, 2),
		myNew(c, NumberOfReactions);
		myNew(c_sum, NumberOfSpecies);
		myNew(a_sum, NumberOfSpecies);
		myNew(c_in, NumberOfSpecies, 0);
		myNew(c_out, NumberOfSpecies, 0);
		myNew(revisable_2nd_reaction, NumberOfSpecies, NumberOfSpecies);
		myNew(number_of_revisable_2nd_reaction, NumberOfSpecies, NumberOfSpecies);		
		myNew(c_over_csum, NumberOfReactions);
		myNew(asum_over_csum, NumberOfSpecies);
		myNew(cx, NumberOfSpecies);
		myNew(exp_ct, NumberOfReactions);
		myNew(one_exp_ct, NumberOfReactions);
		myNew(c2, NumberOfReactions);
		myNew(reactantList, NumberOfReactions, 2);
		myNew(rateList, NumberOfReactions);
		myNew(normalReaction, NumberOfReactions);
		myNew(revisedReaction, NumberOfReactions);
		myNew(pairProduct, NumberOfReactions);
		myNew(pairReactant, NumberOfReactions);
		myNew(oneDirectionReaction, NumberOfReactions);
		myNew(catalyzedReaction, NumberOfReactions);

		for(i=0; i<NumberOfReactions; i++)
		{
			reactantList[i][0]=reactant[i][0];
			reactantList[i][1]=reactant[i][1];
			rateList[i]=rate[i];
		}

		for(i=0; i<NumberOfReactions; i++)
			reactionToSpecies[i].clear();
		for(i=0; i<NumberOfSpecies; i++)
			speciesToReaction[i].clear();
		for(i=0; i<NumberOfReactions; i++)
		{
			if(reactantList[i][0]!=-1)
			{
				reactionToSpecies[i].push_back(reactantList[i][0]);
				speciesToReaction[reactantList[i][0]].push_back(i);
			}
			if(reactantList[i][1]!=-1 && reactantList[i][1]!=reactantList[i][0])
			{
				reactionToSpecies[i].push_back(reactantList[i][1]);
				speciesToReaction[reactantList[i][1]].push_back(i);
			}
		}
		for(i=0; i<NumberOfReactions; i++)
			affectedSpecies[i]=reactionToSpecies[i].size();
		for(i=0; i<NumberOfSpecies; i++)
			affectedReactions[i]=speciesToReaction[i].size();

		//	copy a dense version of the stoichiometry for random access
		for(i=0; i<NumberOfReactions; i++)
		{
			for(j=0; j<NumberOfSpecies; j++)
			{
				stoi(i, j)=stoichiometry(i,j);
				squarestoi(i, j)=stoi(i, j)*stoi(i, j);
			}
		}

		//creat the relatedSpecies, inputReactions, inputStoi, outputReactions, outputStoi, revisable, affectedSpecies
		for(j=0; j<NumberOfSpecies; j++)
		{
			InputReactions.clear();
			OutputReactions.clear();
			for(i=0; i<NumberOfReactions; i++)
			{
				if(stoi(i, j)>0)
				{
					//get relatedSpecies
					for(k=0; k<NumberOfSpecies; k++)
					{
						if(stoi(i,k)<0 && k!=j)
						{
							if(find(relatedSpecies[j].begin(), relatedSpecies[j].end(), k)==relatedSpecies[j].end())
							{
								relatedSpecies[j].push_back(k);
								connectedSpecies[j][k]=true;
							}
						}
					}
					//get input and output reactions and stoichiometries of the species
					InputReactions[j].push_back(i);
				//	inputStoi[j].push_back(int(stoi(i, j)));
					if(stoi(i, j)>1)
						revisable[j]=false;
				}
				//	else if(stoi(i, j)<0)
				if(j==reactantList[i][0] || j==reactantList[i][1])
				{
					//get relatedSpecies
					for(k=0; k<NumberOfSpecies; k++)
					{
						if(stoi(i, k)!=0 && k!=j)
						{
							if(find(influencedSpecies[j].begin(), influencedSpecies[j].end(), k)==influencedSpecies[j].end())
							{
								influencedSpecies[j].push_back(k);
								connectedSpecies[j][k]=true;
								connectedSpecies[k][j]=true;//in case when j is catalyst
							}
							//	if(stoi(i, k)<0)
							if(k==reactantList[i][0] || k==reactantList[i][1])
							{
								if(find(relatedSpecies[j].begin(), relatedSpecies[j].end(), k)==relatedSpecies[j].end())
									relatedSpecies[j].push_back(k);
							}
						}
					}
					//get input and output reactions and stoichiometries of the species
					if(stoi(i, j)<0)
					{
						OutputReactions[j].push_back(i);
					//	outputStoi[j].push_back(int(-stoi(i, j)));
					}
					if(stoi(i, j)<-1)
						revisable[j]=false;
					if(stoi(i, j)>0)//in case of S->2S
						revisable[j]=false;
				}
			}
			NumberOfInputReactions=InputReactions.size();
			myNew(inputReactions[j], InputReactions);
			NumberOfOutputReactions=OutputReactions.size();
			myNew(outputReactions[j], OutputReactions);
			relatedSpeciesSize[j]=relatedSpecies[j].size();
			influencedSpeciesSize[j]=influencedSpecies[j].size();
			connectedSpeciesSize[j]=connectedSpecies[j].size();
		}

		// construct the Jacobian List
		for(i=0; i<NumberOfReactions; i++)
		{
			if(reactantList[i][1]==-1)
				Jacobian[i]=&improvedTauLeaping_aggressive<_denseVectorType,_matrixType,_propensitiesFunctorType,_dependencyGraphType>::singleMolecule;
			else if(reactantList[i][0]!=reactantList[i][1])
				Jacobian[i]=&improvedTauLeaping_aggressive<_denseVectorType,_matrixType,_propensitiesFunctorType,_dependencyGraphType>::differentMolecule;
			else
				Jacobian[i]=&improvedTauLeaping_aggressive<_denseVectorType,_matrixType,_propensitiesFunctorType,_dependencyGraphType>::sameMolecule;
		}

		//reexamine the revisable vector
		for(i=0; i<NumberOfReactions; i++)
		{
			if(reactantList[i][0]==-1)
				continue;
			if(reactantList[i][1]==-1)
			{
				if(revisable[reactantList[i][0]])
				{
					if(stoi(i,reactantList[i][0])==0)
					{
						//push the catalyzed reaction
						Revisable_catalyze_reaction[reactantList[i][0]].push_back(i);

						//the procudt and the catalyze is an impossible pair
						for(j=0; j<NumberOfSpecies; j++)
						{
							if(stoi(i, j)>0)//find a product
							{
								revisable_2nd(reactantList[i][0], j)=-2;
								revisable_2nd(j, reactantList[i][0])=-2;
							}
						}
					}
					else if(stoi(i,reactantList[i][0])!=-1)
						revisable[reactantList[i][0]]=false;
				}
			}
			else if(reactantList[i][0]!=reactantList[i][1])
			{
				if(revisable[reactantList[i][0]])
				{
					if(stoi(i,reactantList[i][0])==0)
					{
						//push the catalyzed reaction
						Revisable_catalyze_reaction[reactantList[i][0]].push_back(i);

						//the procudt and the catalyze is an impossible pair
						for(j=0; j<NumberOfSpecies; j++)
						{
							if(stoi(i, j)>0)//find a product
							{
								revisable_2nd(reactantList[i][0], j)=-2;
								revisable_2nd(j, reactantList[i][0])=-2;
							}
						}
					}
					else if(stoi(i,reactantList[i][0])!=-1)
						revisable[reactantList[i][0]]=false;
				}
				if(revisable[reactantList[i][1]])
				{
					if(stoi(i,reactantList[i][1])==0)
					{
						//push the catalyzed reaction
						Revisable_catalyze_reaction[reactantList[i][1]].push_back(i);

						//the procudt and the catalyze is an impossible pair
						for(j=0; j<NumberOfSpecies; j++)
						{
							if(stoi(i, j)>0)//find a product
							{
								revisable_2nd(reactantList[i][1], j)=-2;
								revisable_2nd(j, reactantList[i][1])=-2;
							}
						}
					}
					else if(stoi(i,reactantList[i][1])!=-1)
						revisable[reactantList[i][1]]=false;
				}
			}
			else
			{
				if(stoi(i,reactantList[i][0])!=-2) //this if-statement may be unnecessary since 2 identical reactants are not included right now
					revisable[reactantList[i][0]]=false;
			}
		}

		//set the products as impossible pair
		std::vector<int> buffer;
		std::vector<int>::iterator Iterator;
		for(i=0; i<NumberOfReactions; i++)
		{
			buffer.clear();
			for(j=0; j<NumberOfSpecies; j++)
			{
				if(stoi(i, j)>0) //find a product
				{
					for(Iterator=buffer.begin(); Iterator!=buffer.end(); Iterator++)
					{
						revisable_2nd((*Iterator), j)=-2;
						revisable_2nd(j, (*Iterator))=-2;
					}
					buffer.push_back(j);
				}
			}
		}

		//allocate memory for PairedParameters
		PairedParameters::NumberOfReactions=NumberOfReactions;
#ifdef DEBUG
		PairedParameters temp;
		std::vector<PairedParameters> tempvector;
		for(i=0; i<NumberOfSpecies; i++)
			tempvector.push_back(temp);
		for(i=0; i<NumberOfSpecies; i++)
			pairedparameter.push_back(tempvector);
#else
		myNew(pairedparameter, NumberOfSpecies, NumberOfSpecies);
#endif
		myNew(revisable_catalyze_reaction, Revisable_catalyze_reaction);
		myNew(number_of_revisable_catalyze_reaction, NumberOfSpecies);
		for(i=0; i<NumberOfSpecies; i++)
			number_of_revisable_catalyze_reaction[i]=Revisable_catalyze_reaction[i].size();

#ifdef PROFILE_SIMULATE_TAU_LEAPING_EXPLICIT_ADAPTIVE
		numberOfIteration=0;
		numberOfEquation=0;
#endif
	}

	templateList void typeList::findPair(int speciesId, int reactionId, std::vector<std::vector<std::vector<int> > > &Revisable_2nd_reaction)
	{
		std::size_t j, k;
		int reversereaction;
		bool isPair;
		for(j=0; j<NumberOfSpecies; j++)//look for product
		{
			if(stoi(reactionId, j)==1)//find a product
			{
				if(revisable[j]==true)//no need to check unrevisable species
				{
					if(stoi(reactionId, speciesId)==0)//speciesId is a catalyst, it forms impossible pair with the product
					{
						revisable_2nd(speciesId, j)=-2;
						revisable_2nd(j, speciesId)=-2;
					}
					else
					{
						//push the reaction id into the revisable_2nd_reaction
						Revisable_2nd_reaction[speciesId][j].push_back(reactionId);
						if(revisable_2nd(speciesId, j)>-1)//already has a reaction that reactant[i][0]->j
						{
							/*	//	mark it as impossible pair
							revisable_2nd(speciesId, j)=-2;
							revisable_2nd(j, speciesId)=-2;*/

							//	mark it as multithread pair
							revisable_2nd(speciesId, j)=-3;
							revisable_2nd(j, speciesId)=-3;
						}
						else if(revisable_2nd(speciesId, j)==-1)//no previous reaction that reactant[i][0]->j
						{
							//check if it is a reverisble reaction
							reversereaction=revisable_2nd(j, speciesId);
							if(reversereaction==-1)
								revisable_2nd(speciesId, j)=reactionId;//assign the reaction id to the pair
							else 
							{
								isPair=true;
								for(k=0; k<NumberOfSpecies; k++)
								{
									if(stoi(reactionId, k)+stoi(reversereaction, k)!=0)//can't cancel each other
									{
										/*	revisable_2nd(speciesId, j)=-2;
										revisable_2nd(j, speciesId)=-2;*/

										revisable_2nd(speciesId, j)=-3;
										revisable_2nd(j, speciesId)=-3;
										isPair=false;
										break;
									}
								}
								if(isPair)
								{
									revisable_2nd(speciesId, j)=reactionId;//assign the reaction id to the pair

									//creat the reverseReaction vector
								//	reverseReaction[reactionId]=reversereaction;
								//	reverseReaction[reversereaction]=reactionId;
								}
							}
						}
					}
				}
				/*	else
				{
				//	mark it as impossible pair
				revisable_2nd(speciesId, j)=-2;
				revisable_2nd(j, speciesId)=-2;
				}*/ //not necessary to assign -2 since these lines will not be used
			}
			/*	else if(stoi(i, j)>1)//more than one molecules are produced in one firing
			{
			revisable_2nd(speciesId, j)=-2;//mark it as impossible pair
			revisable_2nd(j, speciesId)=-2;
			}*/ //not necessary to assign -2 since these lines will not be used
		}
	}

	templateList void typeList::creatPairMatrix()
	{
		std::size_t i, j, k, n;
		std::vector<std::vector<std::vector<int> > >Revisable_2nd_reaction(NumberOfSpecies, std::vector<std::vector<int> > (NumberOfSpecies));
		for(i=0; i<NumberOfReactions; i++)
		{
			if(reactantList[i][0]==-1)
				continue;
			if(reactantList[i][1]==-1)
			{
				if(revisable[reactantList[i][0]]==true)//no need to check unrevisable species
					findPair(reactantList[i][0], i, Revisable_2nd_reaction);
			}
			else if(reactantList[i][0]!=reactantList[i][1])
			{
				if(revisable[reactantList[i][0]]==true)//no need to check unrevisable species
				{
					revisable_2nd(reactantList[i][0], reactantList[i][1])=-2;
					findPair(reactantList[i][0], i, Revisable_2nd_reaction);
				}
				if(revisable[reactantList[i][1]]==true)//no need to check unrevisable species
				{
					revisable_2nd(reactantList[i][1], reactantList[i][0])=-2;
					findPair(reactantList[i][1], i, Revisable_2nd_reaction);
				}
			}
		}
		for(i=0; i<NumberOfSpecies; i++)
		{
			for(j=0; j<NumberOfSpecies; j++)
			{
				number_of_revisable_2nd_reaction[i][j]=Revisable_2nd_reaction[i][j].size();
				myNew(revisable_2nd_reaction[i][j], Revisable_2nd_reaction[i][j]);
			}
		}

		//allocate the forward and backward propensities for each parameter
		std::list<int> temp[2];
		for(i=0; i<NumberOfSpecies; i++)
		{
			if(revisable[i]==true)
			{
				for(j=0; j<NumberOfSpecies; j++)
				{
					if(revisable[j]==true && revisable_2nd(i, j)!=-2)//only prepare the revisable pairs
					{
						//copy the input and output reaction to the parameter instance

						temp[0].clear();
						temp[1].clear();
						std::copy(&inputReactions[j][0], &inputReactions[j][0]+NumberOfInputReactions[j], std::back_inserter(temp[0]) );
						std::copy(&outputReactions[i][0], &outputReactions[i][0]+NumberOfOutputReactions[i], std::back_inserter(temp[1]) );
					//	n=Revisable_2nd_reaction[i][j].size();
						n=number_of_revisable_2nd_reaction[i][j];

						/*	pairedparameter[i][j].forwardPropensities.resize(n);
						pairedparameter[i][j].forwardReactions=Revisable_2nd_reaction[i][j];*/
						myNew(pairedparameter[i][j].forwardPropensities, n);
						pairedparameter[i][j].numberOfForwardReactions=n;
						myNew(pairedparameter[i][j].forwardReactions, Revisable_2nd_reaction[i][j]);
						for(k=0; k<n; k++)
						{
							pairedparameter[i][j].directions[revisable_2nd_reaction[i][j][k]]=1;//forward reaction
							//remove the forward reactions from output and input reaction lists
							temp[0].remove(revisable_2nd_reaction[i][j][k]);
							temp[1].remove(revisable_2nd_reaction[i][j][k]);
						}

						/*	std::copy(temp[0].begin(), temp[0].end(), std::back_inserter( pairedparameter[i][j].inputReactions[1] ) );
						std::copy(temp[1].begin(), temp[1].end(), std::back_inserter( pairedparameter[i][j].outputReactions[0] ) );*/
						pairedparameter[i][j].numberOfInputReactions[1]=temp[0].size();
						myNew(pairedparameter[i][j].inputReactions[1], temp[0]);
						pairedparameter[i][j].numberOfOutputReactions[0]=temp[1].size();
						myNew(pairedparameter[i][j].outputReactions[0], temp[1]);

						temp[0].clear();
						temp[1].clear();
						std::copy(&inputReactions[i][0], &inputReactions[i][0]+NumberOfInputReactions[i], std::back_inserter(temp[0]) );
						std::copy(&outputReactions[j][0], &outputReactions[j][0]+NumberOfOutputReactions[j], std::back_inserter(temp[1]) );
					//	n=revisable_2nd_reaction[j][i].size();
						n=number_of_revisable_2nd_reaction[j][i];

						/*	pairedparameter[i][j].backwardPropensities.resize(n);
						pairedparameter[i][j].backwardReactions=revisable_2nd_reaction[j][i];*/
						myNew(pairedparameter[i][j].backwardPropensities, n);
						pairedparameter[i][j].numberOfBackwardReactions=n;
						myNew(pairedparameter[i][j].backwardReactions, Revisable_2nd_reaction[j][i]);

						for(k=0; k<n; k++)
						{
							pairedparameter[i][j].directions[revisable_2nd_reaction[j][i][k]]=-1;//backward reaction
							//remove the forward reactions from output and input reaction lists
							temp[0].remove(revisable_2nd_reaction[j][i][k]);
							temp[1].remove(revisable_2nd_reaction[j][i][k]);
						}
						/*	std::copy(temp[0].begin(), temp[0].end(), std::back_inserter( pairedparameter[i][j].inputReactions[0] ) );
						std::copy(temp[1].begin(), temp[1].end(), std::back_inserter( pairedparameter[i][j].outputReactions[1] ) );*/
						pairedparameter[i][j].numberOfInputReactions[0]=temp[0].size();
						myNew(pairedparameter[i][j].inputReactions[0], temp[0]);
						pairedparameter[i][j].numberOfOutputReactions[1]=temp[1].size();
						myNew(pairedparameter[i][j].outputReactions[1], temp[1]);

						//allocate the output propensities
						/*	pairedparameter[i][j].numberOfOutputReactions[0]=pairedparameter[i][j].outputReactions[0].size();
						pairedparameter[i][j].outputPropensities[0].resize(pairedparameter[i][j].numberOfOutputReactions[0]);*/
						myNew(pairedparameter[i][j].outputPropensities[0], pairedparameter[i][j].numberOfOutputReactions[0]);

						/*	pairedparameter[i][j].numberOfOutputReactions[1]=pairedparameter[i][j].outputReactions[1].size();
						pairedparameter[i][j].outputPropensities[1].resize(pairedparameter[i][j].numberOfOutputReactions[1]);*/
						myNew(pairedparameter[i][j].outputPropensities[1], pairedparameter[i][j].numberOfOutputReactions[1]);
					}
				}
			}
		}
	}
	
	templateList double typeList::getJacobian(int i, int j)//i:species; j:reaction
	{
		if(i==reactantList[j][0])
			return (this->*(Jacobian[j]))(reactantList[j][1], j);
		else if(i==reactantList[j][1])
			return (this->*(Jacobian[j]))(reactantList[j][0], j);
		else
			return 0;
	} 
	
	templateList void typeList::initialize(double startTime, double endTime)//for neg
	{
		std::size_t i;
		SSA::initialize(startTime);
		criticalSpecies.clear();//it seems empty vector doesn't work for insert_element method. So give it one element '-1' at the beginning
		//	noncriticalSpecies.clear();
		revisedSpecies.clear();
		revisedPair.clear();
		noncriticalSpecies=trimed_list;
		/*	for(i=0; i<NumberOfReactants; i++)
		noncriticalSpecies.push_back(trimed_list[i]);*/
		/*	for(i=0; i<NumberOfReactions; i++)
		tagList[i]=1;*/

		accumulator=0;
		for(i=0; i<NumberOfSpecies; i++)
		{
			reevaluatedSpecies[i]=0;
			//	last_stepsize[i]=endTime;
		}

		//	lastIterationSpecies=-1;
		/*	for(i=0; i<NumberOfReactions; i++)
		reevaluatedReaction[i]=0;*/
	}
	
	templateList void typeList::selectTau(double &noncriticalStepsize, double &criticalStepsize)
	{

		std::size_t i, n;
		std::list<std::size_t>::iterator Iterator;

		//for noncritical species
		axpy_prod(currentPropensities, stoichiometry, mu, true);
		axpy_prod(currentPropensities, squaredVj, sigmaSquared, true);

		double numerator, temp1, temp2;
		noncriticalStepsize=std::numeric_limits<double>::max();

		for(Iterator=trimed_list.begin(); Iterator!=trimed_list.end(); Iterator++)
		{
			numerator=epsilon*currentPopulation[*Iterator]/g;//temp--should use G
			numerator=std::max(numerator,1.0);
			temp1=numerator/fabs(mu(*Iterator));
			temp2=numerator*numerator/sigmaSquared(*Iterator);
			stepsize[*Iterator]=std::min(temp1, temp2);
			noncriticalStepsize=std::min(stepsize[*Iterator], noncriticalStepsize);//the noncritical stepsize return the smallest stepsize
		}
		/*	for(Iterator=criticalSpecies.begin(); Iterator!=criticalSpecies.end(); Iterator++)
		{
		numerator=epsilon*currentPopulation[*Iterator]/g;//temp--should use G
		numerator=std::max(numerator,1.0);
		temp1=numerator/fabs(mu(*Iterator));
		temp2=numerator*numerator/sigmaSquared(*Iterator);
		stepsize[*Iterator]=std::min(temp1, temp2);
		}*/

		//	std::sort(stepsize_order.begin(), stepsize_order.end(), boost::bind(&cmp, _1, _2, this));	
		noncriticalSpecies.sort(boost::bind(&cmp, _1, _2, this));
		criticalSpecies.sort(boost::bind(&cmp, _1, _2, this));

		//find revisable species
		bool revised, pair_revised;
		std::list<std::size_t> *p;
		std::list<std::size_t>::iterator position;
		std::vector<int> pair(2);
		std::list<std::vector<int> >::iterator pairIterator;
		std::size_t species_index;
		PairedParameters *pairPointer;
		while((revisedSpecies.size()+2*revisedPair.size())!=NumberOfReactants)
		{
			if(criticalSpecies.empty())
				p=&noncriticalSpecies;
			else if(noncriticalSpecies.empty())
				p=&criticalSpecies;
			else
				p=stepsize[noncriticalSpecies.front()]<stepsize[criticalSpecies.front()]?&noncriticalSpecies:&criticalSpecies;
			species_index=p->front();
			if(revisable[species_index])
			{
				revised=true;
				pair_revised=false;
				/*	for(i=0; i<connectedSpeciesSize[species_index]; i++)
				{
				if(find(revisedSpecies.begin(), revisedSpecies.end(), connectedSpecies[species_index][i])!=revisedSpecies.end())
				{
				revised=false;
				break;
				}
				}*/
				for(Iterator=revisedSpecies.begin(); Iterator!=revisedSpecies.end(); Iterator++)
				{
					if(connectedSpecies[species_index][*Iterator])
					{
						//	connected to a revised species, further examination required
						if(revisable_2nd(species_index, *Iterator)==-2)//impossible pair
						{
							revised=false;
							break;
						}
						else if(pair_revised==true)//already paired up with another revised species
						{
							revised=false;
							break;
						}
						else//record the pair
						{
							position=Iterator;//record the position of the connected revised species so that we can delete it from the list when we move it to the revisedPair list
							pair[0]=species_index;
							pair[1]=*Iterator;
							pair_revised=true;
						}
						/*	{
						revised=false;
						break;
						}//no second order*/
					}
				}

				//	check the paired revised species
				if(revised)
				{
					for(pairIterator=revisedPair.begin(); pairIterator!=revisedPair.end(); pairIterator++)
					{
						if(connectedSpecies[species_index][(*pairIterator)[0]]||connectedSpecies[species_index][(*pairIterator)[1]])
						{
							//	connected to a revised pair, deny immediately
							revised=false;
							break;
						}
					}
				}

				if(revised)
				{
					if(pair_revised)
					{
						revisedPair.push_back(pair);
						revisedSpecies.erase(position);//remove the paird specied from revisedSpecies list
						//for(i=0; i<affectedReactions[species_index]; i++)
						//	tagList[speciesToReaction[species_index][i]]=3;//3 for pair revised
						//for(i=0; i<affectedReactions[pair[1]]; i++) //change the tag of the paired species
						//	tagList[speciesToReaction[pair[1]][i]]=3;//
						pairPointer=&pairedparameter[pair[0]][pair[1]];					
						for(i=0; i<pairPointer->numberOfOutputReactions[0]; i++)
							tagList[pairPointer->outputReactions[0][i]]=3;//3 for pair revised
						for(i=0; i<pairPointer->numberOfOutputReactions[1]; i++) //change the tag of the paired species
							tagList[pairPointer->outputReactions[1][i]]=3;//

						//the reaction connects the pair is labeled 4
					//	n=revisable_2nd_reaction[pair[0]][pair[1]].size();
						n=number_of_revisable_2nd_reaction[pair[0]][pair[1]];
						for(i=0; i<n; i++)
							tagList[revisable_2nd_reaction[pair[0]][pair[1]][i]]=6;
					//	n=revisable_2nd_reaction[pair[1]][pair[0]].size();
						n=number_of_revisable_2nd_reaction[pair[1]][pair[0]];
						for(i=0; i<n; i++)
							tagList[revisable_2nd_reaction[pair[1]][pair[0]][i]]=6;

						if(revisable_2nd(pair[0], pair[1])>-1)
							tagList[revisable_2nd(pair[0], pair[1])]=4;
						if(revisable_2nd(pair[1], pair[0])>-1)
						{
							//if two way, assign 5 to both of them
							if(revisable_2nd(pair[0], pair[1])>-1)
							{
								tagList[revisable_2nd(pair[0], pair[1])]=5;
								tagList[revisable_2nd(pair[1], pair[0])]=5;
							}
							else
								tagList[revisable_2nd(pair[1], pair[0])]=4;
						}

						//the catalyzed reactions are labeled 6
					//	n=revisable_catalyze_reaction[pair[0]].size();
						for(i=0; i<number_of_revisable_catalyze_reaction[pair[0]]; i++)
							tagList[revisable_catalyze_reaction[pair[0]][i]]=6;
					//	n=revisable_catalyze_reaction[pair[1]].size();
						for(i=0; i<number_of_revisable_catalyze_reaction[pair[1]]; i++)
							tagList[revisable_catalyze_reaction[pair[1]][i]]=6;
					}
					else
					{
						revisedSpecies.push_back(species_index);
						for(i=0; i<affectedReactions[species_index]; i++)
							tagList[speciesToReaction[species_index][i]]=2;

						//the catalyzed reactions are labeled 7
					//	n=revisable_catalyze_reaction[species_index].size();
						for(i=0; i<number_of_revisable_catalyze_reaction[species_index]; i++)
							tagList[revisable_catalyze_reaction[species_index][i]]=7;
					}
					p->pop_front();
				}
				else
					break;
			}
			else
				break;
		}
		//	noncriticalStepsize=stepsize(noncriticalSpecies.front());

		//for critical reactions
		//put back the critical tag
		for(Iterator=criticalSpecies.begin(); Iterator!=criticalSpecies.end(); Iterator++)
		{
		//	n=outputReactions[*Iterator].size();
			n=NumberOfOutputReactions[*Iterator];
			for(i=0; i<n; i++)//modify tag for each reaction that take this species as a reactant if neccessary
				tagList[outputReactions[*Iterator][i]]=0;
		}
		criticalPropensitySum=0;
		for (i=0; i<NumberOfReactions; i++)
		{
			if(tagList[i]==0)
				criticalPropensitySum+=currentPropensities[i];
		}
		if(criticalPropensitySum!=0)
			criticalStepsize=randomGenerator.Exponential(1.0/criticalPropensitySum);
		else
			criticalStepsize=-1;

		//merge noncritical and critical species
		mergedSpecies=noncriticalSpecies;
		//	std::list<std::size_t> temp;
		//	temp=criticalSpecies;
		//	mergedSpecies.merge(temp, boost::bind(&cmp, _1, _2, this));
		mergedSpecies.insert(mergedSpecies.end(), criticalSpecies.begin(), criticalSpecies.end());
		//	noncriticalStepsize=stepsize[mergedSpecies.front()];
	}
	
	templateList int typeList::selectReactions(double leapSize, bool runCritical)
	{
		std::size_t i, j, n;
		int previousReactionIndex=-1;
		if(runCritical)//handle both critical and noncritical reactions
		{
			//for critical
			//generate a uniform random number between (0,propensitySum)
			double r=randomGenerator.ContinuousOpen(0,1)*criticalPropensitySum;
			double jsum=0;

			for (i=0; i!=NumberOfReactions; ++i) {
				//for critical
				if(tagList[i]==0)
				{
					previousReactionCounts[i]=0;//not noncritical so set 0
					if(jsum<r)
					{
						previousReactionIndex=i;
						jsum+=currentPropensities[i];
					}
				}
				//for noncritical
				else if(tagList[i]==1)
					previousReactionCounts[i]=randomGenerator.Poisson(leapSize*currentPropensities[i]);
			}
		}
		else//only handle noncritical reactions
		{
			for (i=0; i!=NumberOfReactions; ++i) {
				if(tagList[i]==0)//set 0 since in last round it may be nonzero.
					previousReactionCounts[i]=0;
				else if(tagList[i]==1)
					previousReactionCounts[i]=randomGenerator.Poisson(leapSize*currentPropensities[i]);
			}
		}

		//generate the revised reactions
		std::list<std::size_t>::iterator Iterator;
		//	std::vector<std::size_t> outputlist;
		//	std::vector<double> p;
		double inputNumber, outputNumber, leftNumber, temp1, temp2, temp, EX_T;//, lambda, mu;
		for(Iterator=revisedSpecies.begin(); Iterator!=revisedSpecies.end(); Iterator++)
		{
			//	n=torevisedreaction[*Iterator].size();
			if(c_sum[*Iterator]>0 || number_of_revisable_catalyze_reaction[*Iterator]>0)
			{
				temp1=c_sum[*Iterator]*leapSize;
				temp2=exp(-temp1);

				//	n=outputlist.size();
				if(c_sum[*Iterator]>0)
				{
					inputNumber=0;
				//	n=inputReactions[*Iterator].size();
					n=NumberOfInputReactions[*Iterator];
					//	lambda=0;//lambda should be a_sum[*Iterator]
					for(i=0;i<n;i++)
					{
						inputNumber+=previousReactionCounts[inputReactions[*Iterator][i]];
						//		lambda+=currentPropensities[inputReactions[*Iterator][i]];//not *stoi since only one product produced for one fire
					}
					/*	leftNumber=randomGenerator.Binomial(currentPopulation[*Iterator], exp(-c_sum[*Iterator]*leapSize));

					//	leftNumber+=std::min(inputNumber, randomGenerator.Poisson(a_sum[*Iterator]/c_sum[*Iterator]*(1-exp(-c_sum[*Iterator]*leapSize))));
					leftNumber+=randomGenerator.Binomial(inputNumber, (1-exp(-c_sum[*Iterator]*leapSize))/(c_sum[*Iterator]*leapSize));*/
					//	temp1=c_sum[*Iterator]*leapSize;
					//	temp2=exp(-temp1);
				//	n=outputReactions[*Iterator].size();
					n=NumberOfOutputReactions[*Iterator];
					leftNumber=randomGenerator.Binomial(currentPopulation[*Iterator], temp2);//left number is the population at end time
					leftNumber+=randomGenerator.Binomial(inputNumber, (1-temp2)/temp1);

					outputNumber=inputNumber+currentPopulation[*Iterator]-leftNumber;

					//distribute the output species
					temp=c_sum[*Iterator];//copy c_sum[*Iterator] and restore it after the for loop
					for(i=0; i<n-1; i++)
					{
						//	previousReactionCounts[torevisedreaction[*Iterator][i]]=randomGenerator.Binomial(outputNumber, c[torevisedreaction[*Iterator][i]]/c_sum[*Iterator]);
						previousReactionCounts[outputReactions[*Iterator][i]]=randomGenerator.Binomial(outputNumber, c[outputReactions[*Iterator][i]]/c_sum[*Iterator]);
						//	outputNumber-=previousReactionCounts[torevisedreaction[*Iterator][i]];
						outputNumber-=previousReactionCounts[outputReactions[*Iterator][i]];
						//	c_sum[*Iterator]-=c[torevisedreaction[*Iterator][i]];
						c_sum[*Iterator]-=c[outputReactions[*Iterator][i]];
					}
					//	previousReactionCounts[torevisedreaction[*Iterator][n-1]]=outputNumber;
					previousReactionCounts[outputReactions[*Iterator][n-1]]=outputNumber;
					c_sum[*Iterator]=temp;
				}

				//the catalyzed reactions
			//	n=revisable_catalyze_reaction[*Iterator].size();
				if(number_of_revisable_catalyze_reaction[*Iterator]>0)
				{
					EX_T=currentPopulation[*Iterator]*temp2+(1-temp2)*asum_over_csum[*Iterator];
					if(c_sum[*Iterator]>0)
						temp1=(1-temp2)/c_sum[*Iterator];
					else
						temp1=leapSize;
					temp2=currentPopulation[*Iterator]*temp1+asum_over_csum[*Iterator]*(leapSize-temp1);
					temp2+=leapSize/2*(leftNumber-EX_T);
					for(i=0; i<number_of_revisable_catalyze_reaction[*Iterator]; i++)
						previousReactionCounts[revisable_catalyze_reaction[*Iterator][i]]=randomGenerator.Poisson(temp2*c[revisable_catalyze_reaction[*Iterator][i]]);
				}
			}
		}

		//generate the revised pair reactions
		std::list<std::vector<int> >::iterator pairIterator;
		//	std::vector<std::size_t> outputlist;
		//	std::vector<double> p;
		int forwardReaction, backwardReaction;
		double inputPair[2], leftPair[2], stayPair[2];
		double forwardCounts, backwardCounts;//the number forward and backward reactions from the sampling of the distributions
		double forwardSample, backwardSample;//the number of one directions reactions sampled from the averaged population
		std::vector<double> multinomialSample;
		std::vector<int>::iterator vIterator;
		PairedParameters *parameterPointer;
		for(pairIterator=revisedPair.begin(); pairIterator!=revisedPair.end(); pairIterator++)
		{
			parameterPointer=&(pairedparameter[(*pairIterator)[0]][(*pairIterator)[1]]);
			forwardReaction=revisable_2nd((*pairIterator)[0], (*pairIterator)[1]);
			backwardReaction=revisable_2nd((*pairIterator)[1], (*pairIterator)[0]);

			//the input molecules
			inputPair[0]=0;
			inputPair[1]=0;
			for(i=0; i<parameterPointer->numberOfInputReactions[0]; i++)
				inputPair[0]+=previousReactionCounts[parameterPointer->inputReactions[0][i]];
			for(i=0; i<parameterPointer->numberOfInputReactions[1]; i++)
				inputPair[1]+=previousReactionCounts[parameterPointer->inputReactions[1][i]];

			//update the probabilities
			parameterPointer->update_full_time_dependent(leapSize);

			//	sample
			//multinomial for n0
			forwardCounts=0;
			backwardCounts=0;
			multinomialSample=randomGenerator.Multinomial(currentPopulation[(*pairIterator)[0]], parameterPointer->p[0], 4);
			//if(forwardReaction>-1)
			//	previousReactionCounts[forwardReaction]=multinomialSample[1]+multinomialSample[3];
			if(parameterPointer->c12>0)
				forwardCounts+=multinomialSample[1]+multinomialSample[3];
			leftPair[0]=multinomialSample[2];
			leftPair[1]=multinomialSample[3];
			stayPair[0]=multinomialSample[0];
			stayPair[1]=multinomialSample[1];

			//multinomial for n1
			//	multinomialSample=randomGenerator.Multinomial(currentPopulation[(*pairIterator)[1]], parameterPointer->p[1]);
			randomGenerator.Multinomial(currentPopulation[(*pairIterator)[1]], parameterPointer->p[1], 4, multinomialSample);
			//if(backwardReaction>-1)
			//	previousReactionCounts[backwardReaction]=multinomialSample[0]+multinomialSample[2];
			if(parameterPointer->c21>0)
				backwardCounts+=multinomialSample[0]+multinomialSample[2];
			leftPair[0]+=multinomialSample[2];
			leftPair[1]+=multinomialSample[3];
			stayPair[0]+=multinomialSample[0];
			stayPair[1]+=multinomialSample[1];

			//multinomial for a0
			//	multinomialSample=randomGenerator.Multinomial(inputPair[0], parameterPointer->lambda2p[0]);
			randomGenerator.Multinomial(inputPair[0], parameterPointer->lambda2p[0], 4, multinomialSample);
			//if(forwardReaction>-1)
			//	previousReactionCounts[forwardReaction]+=multinomialSample[1]+multinomialSample[3];
			if(parameterPointer->c12>0)
				forwardCounts+=multinomialSample[1]+multinomialSample[3];
			leftPair[0]+=multinomialSample[2];
			leftPair[1]+=multinomialSample[3];
			stayPair[0]+=multinomialSample[0];
			stayPair[1]+=multinomialSample[1];

			//multinomial for a1
			//	multinomialSample=randomGenerator.Multinomial(inputPair[1], parameterPointer->lambda2p[1]);
			randomGenerator.Multinomial(inputPair[1], parameterPointer->lambda2p[1], 4, multinomialSample);
			//if(backwardReaction>-1)
			//	previousReactionCounts[backwardReaction]+=multinomialSample[0]+multinomialSample[2];
			if(parameterPointer->c21>0)
				backwardCounts+=multinomialSample[0]+multinomialSample[2];
			leftPair[0]+=multinomialSample[2];
			leftPair[1]+=multinomialSample[3];
			stayPair[0]+=multinomialSample[0];
			stayPair[1]+=multinomialSample[1];

			for(j=0; j<2; j++)
			{
#ifdef DEBUG
				multinomialSample=randomGenerator.Multinomial(leftPair[j], parameterPointer->outputPropensities[j]);
#else
				multinomialSample=randomGenerator.Multinomial(leftPair[j], parameterPointer->outputPropensities[j], parameterPointer->numberOfOutputReactions[j]);
#endif
				for(i=0; i<parameterPointer->numberOfOutputReactions[j]; i++)
					previousReactionCounts[parameterPointer->outputReactions[j][i]]=multinomialSample[i];
			}

			//distribute the forward reaction counts
			temp1=std::max(parameterPointer->int_lambda[0]+leapSize/2*(stayPair[0]-parameterPointer->EX_T[0]),0.0);
			temp2=std::max(parameterPointer->int_lambda[1]+leapSize/2*(stayPair[1]-parameterPointer->EX_T[1]),0.0);
			if(forwardCounts<=backwardCounts)//sample forward first, then compute backward
			{
				//	n=parameterPointer->forwardReactions.size();
				if(parameterPointer->numberOfForwardReactions>0)
				{
					forwardSample=randomGenerator.Poisson(temp1*parameterPointer->c12);
#ifdef DEBUG
					multinomialSample=randomGenerator.Multinomial(forwardSample, parameterPointer->forwardPropensities);
#else
					multinomialSample=randomGenerator.Multinomial(forwardSample, parameterPointer->forwardPropensities, parameterPointer->numberOfForwardReactions);
#endif
					for(i=0; i<parameterPointer->numberOfForwardReactions; i++)
						previousReactionCounts[parameterPointer->forwardReactions[i]]=multinomialSample[i];
				}
				else
					forwardSample=0;

				//distribute the backward reaction counts
				//	n=parameterPointer->backwardReactions.size();
				if(parameterPointer->numberOfBackwardReactions>0)
				{
					backwardSample=backwardCounts-forwardCounts+forwardSample;
#ifdef DEBUG
					multinomialSample=randomGenerator.Multinomial(backwardSample, parameterPointer->backwardPropensities);
#else
					multinomialSample=randomGenerator.Multinomial(backwardSample, parameterPointer->backwardPropensities, parameterPointer->numberOfBackwardReactions);
#endif
					for(i=0; i<parameterPointer->numberOfBackwardReactions; i++)
						previousReactionCounts[parameterPointer->backwardReactions[i]]=multinomialSample[i];
				}
			}
			else//sample backward first, then compute forward
			{
				//	n=parameterPointer->backwardReactions.size();
				if(parameterPointer->numberOfBackwardReactions>0)
				{
					backwardSample=randomGenerator.Poisson(temp2*parameterPointer->c21);
#ifdef DEBUG
					multinomialSample=randomGenerator.Multinomial(backwardSample, parameterPointer->backwardPropensities);
#else
					multinomialSample=randomGenerator.Multinomial(backwardSample, parameterPointer->backwardPropensities, parameterPointer->numberOfBackwardReactions);
#endif
					for(i=0; i<parameterPointer->numberOfBackwardReactions; i++)
						previousReactionCounts[parameterPointer->backwardReactions[i]]=multinomialSample[i];
				}
				else
					backwardSample=0;

				//distribute the forwardward reaction counts
				//	n=parameterPointer->forwardReactions.size();
				if(parameterPointer->numberOfForwardReactions>0)
				{
					forwardSample=forwardCounts-backwardCounts+backwardSample;
#ifdef DEBUG
					multinomialSample=randomGenerator.Multinomial(forwardSample, parameterPointer->forwardPropensities);
#else
					multinomialSample=randomGenerator.Multinomial(forwardSample, parameterPointer->forwardPropensities, parameterPointer->numberOfForwardReactions);
#endif
					for(i=0; i<parameterPointer->numberOfForwardReactions; i++)
						previousReactionCounts[parameterPointer->forwardReactions[i]]=multinomialSample[i];
				}
			}

			//the catalyzed reactions
		//	n=revisable_catalyze_reaction[(*pairIterator)[0]].size();
			for(i=0; i<number_of_revisable_catalyze_reaction[(*pairIterator)[0]]; i++)
				previousReactionCounts[revisable_catalyze_reaction[(*pairIterator)[0]][i]]=randomGenerator.Poisson(temp1*parameterPointer->Propensities[revisable_catalyze_reaction[(*pairIterator)[0]][i]]);
		//	n=revisable_catalyze_reaction[(*pairIterator)[1]].size();
			for(i=0; i<number_of_revisable_catalyze_reaction[(*pairIterator)[1]]; i++)
				previousReactionCounts[revisable_catalyze_reaction[(*pairIterator)[1]][i]]=randomGenerator.Poisson(temp2*parameterPointer->Propensities[revisable_catalyze_reaction[(*pairIterator)[1]][i]]);

		}
		return previousReactionIndex;
	}
	
	templateList void typeList::updateTagLists()
	{
		std::size_t i, j;
		//	bool changeTag;
		std::list<std::size_t>::iterator Iterator;
		std::list<std::vector<int> >::iterator pairIterator;

		//from noncritical species to critical species
		for(Iterator=noncriticalSpecies.begin(); Iterator!=noncriticalSpecies.end();)
		{
			if(currentPopulation[*Iterator]<=criticalThreshold)
			{
				criticalSpecies.push_back(*Iterator);//add to critical species list
				Iterator=noncriticalSpecies.erase(Iterator);
			}
			else
				Iterator++;
		}

		//from revised to critical and noncritical
		for(Iterator=revisedSpecies.begin(); Iterator!=revisedSpecies.end(); Iterator++)
		{
			if(currentPopulation[*Iterator]>criticalThreshold)
			{
				noncriticalSpecies.push_back(*Iterator);//add to noncritical species list
				//	for(i=0; i<affectedReactions[*Iterator]; i++)//modify tag for each reaction that take this species as a reactant if neccessary
				//		tagList[speciesToReaction[*Iterator][i]]=1;
			}
			else
				criticalSpecies.push_back(*Iterator);
		}
		revisedSpecies.clear();

		//from revisedpair to critical and noncritical
		for(pairIterator=revisedPair.begin(); pairIterator!=revisedPair.end(); pairIterator++)
		{
			for(j=0; j<2; j++)
			{
				if(currentPopulation[(*pairIterator)[j]]>criticalThreshold)
				{
					noncriticalSpecies.push_back((*pairIterator)[j]);//add to noncritical species list
					//	for(i=0; i<affectedReactions[(*pairIterator)[j]]; i++)//modify tag for each reaction that take this species as a reactant if neccessary
					//		tagList[speciesToReaction[(*pairIterator)[j]][i]]=1;
				}
				else
					criticalSpecies.push_back((*pairIterator)[j]);
			}
		}
		revisedPair.clear();

		//from critical species to noncritical species
		for(Iterator=criticalSpecies.begin(); Iterator!=criticalSpecies.end();)
		{
			if(currentPopulation[*Iterator]>criticalThreshold)
			{
				noncriticalSpecies.push_back(*Iterator);//add to noncritical species list
				//	for(i=0; i<affectedReactions[*Iterator]; i++)//modify tag for each reaction that take this species as a reactant if neccessary
				//		tagList[speciesToReaction[*Iterator][i]]=1;
				Iterator=criticalSpecies.erase(Iterator);//erase from critical list
			}
			else
				Iterator++;
		}

		for(i=0; i<NumberOfReactions; i++)
			tagList[i]=1; //initialize all the taglist with 1 (noncritical reaction)

		//update the taglist, will be done in the selecttau after the revised species are generated
		/*	for(Iterator=criticalSpecies.begin(); Iterator!=criticalSpecies.end(); Iterator++)
		{
		for(i=0; i<affectedReactions[*Iterator]; i++)//modify tag for each reaction that take this species as a reactant if neccessary
		tagList[speciesToReaction[*Iterator][i]]=0;
		}*/
	}
	
	templateList void typeList::getReevaluated()
	{
		std::size_t i, j;
		std::list<std::size_t>::iterator Iterator;
		std::list<std::vector<int> >::iterator pairIterator;

		//get the reevaluatedSpecies from revisedSpecies
		for(Iterator=revisedSpecies.begin(); Iterator!=revisedSpecies.end(); Iterator++)
		{
			for(i=0; i<influencedSpeciesSize[*Iterator]; i++)
				reevaluatedSpecies[influencedSpecies[*Iterator][i]]=accumulator;
			for(i=0; i<affectedReactions[*Iterator]; i++)
			{
				//	reevaluatedReaction[speciesToReaction[*Iterator][i]]=accumulator;
				responsiblespeices[speciesToReaction[*Iterator][i]]=*Iterator;
			}
			reevaluatedSpecies[*Iterator]=-1;//may not necessary
		}

		//get the reevaluatedSpecies from revisedPire
		for(pairIterator=revisedPair.begin(); pairIterator!=revisedPair.end(); pairIterator++)
		{
			for(j=0; j<2; j++)
			{
				for(i=0; i<influencedSpeciesSize[(*pairIterator)[j]]; i++)
					reevaluatedSpecies[influencedSpecies[(*pairIterator)[j]][i]]=accumulator;
				for(i=0; i<affectedReactions[(*pairIterator)[j]]; i++)
				{
					//	reevaluatedReaction[speciesToReaction[*Iterator][i]]=accumulator;
					responsiblespeices[speciesToReaction[(*pairIterator)[j]][i]]=(*pairIterator)[j];
					responsiblepair[speciesToReaction[(*pairIterator)[j]][i]]=(*pairIterator);
				}
			}

			//this may not necessary
			reevaluatedSpecies[(*pairIterator)[0]]=-1;
			reevaluatedSpecies[(*pairIterator)[1]]=-1;
		}
	}

	templateList void typeList::getParameters()
	{
		std::size_t i, j, n;
		std::list<std::size_t>::iterator Iterator;
		std::list<std::vector<int> >::iterator pairIterator;


		//for revisedSpecies
		for(Iterator=revisedSpecies.begin(); Iterator!=revisedSpecies.end(); Iterator++)
		{
			//a_sum
		//	n=inputReactions[*Iterator].size();
			n=NumberOfInputReactions[*Iterator];
			a_sum[*Iterator]=0;
			for(i=0; i<n; i++)
			{
				if(tagList[inputReactions[*Iterator][i]]!=0)
					a_sum[*Iterator]+=currentPropensities[inputReactions[*Iterator][i]];
			}

			//c_sum and c
			//	torevisedreaction[*Iterator].clear();
		//	n=outputReactions[*Iterator].size();
			n=NumberOfOutputReactions[*Iterator];
			c_sum[*Iterator]=0;
			for(i=0; i<n; i++)
			{
				if(tagList[outputReactions[*Iterator][i]]!=0)
				{
					//	torevisedreaction[*Iterator].push_back(outputReactions[*Iterator][i]);
					c[outputReactions[*Iterator][i]]=getJacobian(*Iterator, outputReactions[*Iterator][i]);
					c_sum[*Iterator]+=c[outputReactions[*Iterator][i]];
				}
				else c[outputReactions[*Iterator][i]]=0;
			}

			//catalized reactions
		//	n=revisable_catalyze_reaction[*Iterator].size();
			for(i=0;i<number_of_revisable_catalyze_reaction[*Iterator];i++)
			{
				if(tagList[revisable_catalyze_reaction[*Iterator][i]]==0)
					c[revisable_catalyze_reaction[*Iterator][i]]=0;
				else
					c[revisable_catalyze_reaction[*Iterator][i]]=getJacobian(*Iterator, revisable_catalyze_reaction[*Iterator][i]);
			}
		}

		//for revisedPair
		PairedParameters *parameter;
		std::vector<int>::iterator vIterator;
		for(pairIterator=revisedPair.begin(); pairIterator!=revisedPair.end(); pairIterator++)
		{
			parameter=&(pairedparameter[(*pairIterator)[0]][(*pairIterator)[1]]);
			for(j=0; j<2; j++)
			{
				//a_sum
				a_sum[(*pairIterator)[j]]=0;
				for(i=0; i<parameter->numberOfInputReactions[j]; i++)
				{
					if(tagList[parameter->inputReactions[j][i]]==1)
						a_sum[(*pairIterator)[j]]+=currentPropensities[parameter->inputReactions[j][i]];
				}

				parameter->a[j]=a_sum[(*pairIterator)[j]];

				//c_sum and c
				//	torevisedreaction[(*pairIterator)[j]].clear();
				c_sum[(*pairIterator)[j]]=0;
				for(i=0; i<parameter->numberOfOutputReactions[j]; i++)
				{
					if(tagList[parameter->outputReactions[j][i]]==3)
					{
						c[parameter->outputReactions[j][i]]=getJacobian((*pairIterator)[j], parameter->outputReactions[j][i]);
						parameter->outputPropensities[j][i]=c[parameter->outputReactions[j][i]];
						c_sum[(*pairIterator)[j]]+=parameter->outputPropensities[j][i];
					}
					else
					{
						c[parameter->outputReactions[j][i]]=0;
						parameter->outputPropensities[j][i]=0;
					}
				}
			}
			parameter->c1=c_sum[(*pairIterator)[0]];
			parameter->c2=c_sum[(*pairIterator)[1]];

			//n0 and n1
			parameter->n0=currentPopulation[(*pairIterator)[0]];
			parameter->n1=currentPopulation[(*pairIterator)[1]];

			//c12 and c21
			parameter->c12=0;
			if(revisable_2nd((*pairIterator)[0], (*pairIterator)[1])!=-1) 
			{
			//	n=revisable_2nd_reaction[(*pairIterator)[0]][(*pairIterator)[1]].size();
				n=number_of_revisable_2nd_reaction[(*pairIterator)[0]][(*pairIterator)[1]];
				for(i=0; i<n; i++)
				{
					if(tagList[revisable_2nd_reaction[(*pairIterator)[0]][(*pairIterator)[1]][i]]==0)
					{
						parameter->forwardPropensities[i]=0;
						//	c[revisable_2nd_reaction[(*pairIterator)[0]][(*pairIterator)[1]][i]]=0;
						parameter->Propensities[revisable_2nd_reaction[(*pairIterator)[0]][(*pairIterator)[1]][i]]=0;
					}
					else
					{
						parameter->forwardPropensities[i]=getJacobian((*pairIterator)[0], revisable_2nd_reaction[(*pairIterator)[0]][(*pairIterator)[1]][i]);
						//	c[revisable_2nd_reaction[(*pairIterator)[0]][(*pairIterator)[1]][i]]=parameter->forwardPropensities[i];
						parameter->Propensities[revisable_2nd_reaction[(*pairIterator)[0]][(*pairIterator)[1]][i]]=parameter->forwardPropensities[i];
						parameter->c12+=parameter->forwardPropensities[i];
					}
				}
			}
			parameter->c21=0;
			if(revisable_2nd((*pairIterator)[1], (*pairIterator)[0])!=-1)
			{
			//	n=revisable_2nd_reaction[(*pairIterator)[1]][(*pairIterator)[0]].size();
				n=number_of_revisable_2nd_reaction[(*pairIterator)[1]][(*pairIterator)[0]];
				for(i=0; i<n; i++)
				{
					if(tagList[revisable_2nd_reaction[(*pairIterator)[1]][(*pairIterator)[0]][i]]==0)
					{
						parameter->backwardPropensities[i]=0;
						//	c[revisable_2nd_reaction[(*pairIterator)[1]][(*pairIterator)[0]][i]]=0;
						parameter->Propensities[revisable_2nd_reaction[(*pairIterator)[1]][(*pairIterator)[0]][i]]=0;
					}
					else
					{
						parameter->backwardPropensities[i]=getJacobian((*pairIterator)[1], revisable_2nd_reaction[(*pairIterator)[1]][(*pairIterator)[0]][i]);
						//	c[revisable_2nd_reaction[(*pairIterator)[1]][(*pairIterator)[0]][i]]=parameter->backwardPropensities[i];
						parameter->Propensities[revisable_2nd_reaction[(*pairIterator)[1]][(*pairIterator)[0]][i]]=parameter->backwardPropensities[i];
						parameter->c21+=parameter->backwardPropensities[i];
					}
				}
			}

			//the catalyzed reactions
			for(j=0; j<2; j++)
			{
			//	n=revisable_catalyze_reaction[(*pairIterator)[j]].size();
				for(i=0;i<number_of_revisable_catalyze_reaction[(*pairIterator)[j]];i++)
				{
					if(tagList[revisable_catalyze_reaction[(*pairIterator)[j]][i]]==0)
						parameter->Propensities[revisable_catalyze_reaction[(*pairIterator)[j]][i]]=0;
					else
						parameter->Propensities[revisable_catalyze_reaction[(*pairIterator)[j]][i]]=getJacobian((*pairIterator)[j], revisable_catalyze_reaction[(*pairIterator)[j]][i]);
				}
			}

			//update the time invariant paired parameters
			pairedparameter[(*pairIterator)[0]][(*pairIterator)[1]].update_full_time_invariant();
		}
	}

	templateList double typeList::fmean(std::size_t i, double t)
	{
		/*	double exp_ct=exp(-c_sum[responsiblespeices[i]]*t);
		return c[i]/c_sum[responsiblespeices[i]]*(currentPopulation[responsiblespeices[i]]*(1-exp_ct)+a_sum[responsiblespeices[i]]/c_sum[responsiblespeices[i]]*(c_sum[responsiblespeices[i]]*t-1+exp_ct));*/
		return c_over_csum[i]*(currentPopulation[responsiblespeices[i]]*one_exp_ct[i]+asum_over_csum[responsiblespeices[i]]*(c_sum[responsiblespeices[i]]*t-one_exp_ct[i]));
	}
	
	templateList double typeList::fcatalyst(std::size_t i, double t, double &dmean)
	{
		double e_ct=exp(-c_sum[responsiblespeices[i]]*t);
		double temp=(1-e_ct)/c_sum[responsiblespeices[i]];
		dmean=(currentPopulation[responsiblespeices[i]]*e_ct+a_sum[responsiblespeices[i]]*temp)*c[i];
		return (currentPopulation[responsiblespeices[i]]*temp+asum_over_csum[responsiblespeices[i]]*(t-temp))*c[i];
	}
	
	templateList double typeList::varcatalyst(std::size_t i, double t, double &dvar)
	{
		double e_ct=exp(-c_sum[responsiblespeices[i]]*t);
		double temp=(1-e_ct)/c_sum[responsiblespeices[i]];
		double tempdmean=(currentPopulation[responsiblespeices[i]]*e_ct+a_sum[responsiblespeices[i]]*temp)*c[i];
		double tempmean=(currentPopulation[responsiblespeices[i]]*temp+asum_over_csum[responsiblespeices[i]]*(t-temp))*c[i];
		double varX=currentPopulation[responsiblespeices[i]]*e_ct*(1-e_ct)+a_sum[responsiblespeices[i]]*temp;
		double dvarX=currentPopulation[responsiblespeices[i]]*e_ct*(1-2*e_ct)*(-c_sum[responsiblespeices[i]])+a_sum[responsiblespeices[i]]*e_ct;
		double ctemp=pow(c[i]*t/2, 2);
		dvar=tempdmean+ctemp*dvarX+ctemp*2/t*varX;
		return tempmean+ctemp*varX;
	}
	
	templateList double typeList::dmean(std::size_t i, double t)
	{
		return c_over_csum[i]*(cx[responsiblespeices[i]]*exp_ct[i]+a_sum[responsiblespeices[i]]*one_exp_ct[i]);
	}
	
	templateList double typeList::fvar(std::size_t i, double t)
	{
		return c_over_csum[i]*(currentPopulation[responsiblespeices[i]]*one_exp_ct[i]*(1-c_over_csum[i]+c_over_csum[i]*exp_ct[i])+asum_over_csum[responsiblespeices[i]]*(c_sum[responsiblespeices[i]]*t-one_exp_ct[i]));
	}
	
	templateList double typeList::dvar(std::size_t i, double t)
	{
		return c_over_csum[i]*(currentPopulation[responsiblespeices[i]]*exp_ct[i]*(c_sum[responsiblespeices[i]]-c2[i]*one_exp_ct[i])+a_sum[responsiblespeices[i]]*one_exp_ct[i]);
	}
	
	templateList void typeList::form_formular(std::size_t i)
	{
	//	normalReaction.clear();
		NumberOfNormalReactions=0;
	//	revisedReaction.clear();
		NumberOfRevisedReactions=0;
	//	pairProduct.clear();
		NumberOfPairProducts=0;
	//	pairReactant.clear();
		NumberOfPairReactants=0;
	//	oneDirectionReaction.clear();
		NumberOfOneDirectionReactions=0;
	//	catalyzedReaction.clear();
		NumberOfCatalyzedReactions=0;

		std::size_t j, n;
	//	n=inputReactions[i].size();
		n=NumberOfInputReactions[i];

		//go over the inputReactions
		for(j=0; j<n; j++)
		{
			switch(tagList[inputReactions[i][j]])
			{
			case 1:
			//	normalReaction.push_back(inputReactions[i][j]);
				normalReaction[NumberOfNormalReactions]=inputReactions[i][j];
				NumberOfNormalReactions++;
				break;
			case 2:
			//	revisedReaction.push_back(inputReactions[i][j]);
				revisedReaction[NumberOfRevisedReactions]=inputReactions[i][j];
				NumberOfRevisedReactions++;
				break;
			case 3:
			//	pairProduct.push_back(inputReactions[i][j]);
				pairProduct[NumberOfPairProducts]=inputReactions[i][j];
				NumberOfPairProducts++;
				break;
			case 4:
			//	pairReactant.push_back(inputReactions[i][j]);
				pairReactant[NumberOfPairReactants]=inputReactions[i][j];
				NumberOfPairReactants++;
				break;
			case 5:
				//catalyst won't get contrubution from the catalyzed reaction
				if(reactantList[inputReactions[i][j]][0]!=i&&reactantList[inputReactions[i][j]][1]!=i)
				{
				//	pairReactant.push_back(inputReactions[i][j]);
					pairReactant[NumberOfPairReactants]=inputReactions[i][j];
					NumberOfPairReactants++;
				}
				break;
			case 6:
			//	oneDirectionReaction.push_back(inputReactions[i][j]);
				oneDirectionReaction[NumberOfOneDirectionReactions]=inputReactions[i][j];
				NumberOfOneDirectionReactions++;
				break;
			case 7:
			//	catalyzedReaction.push_back(inputReactions[i][j]);
				catalyzedReaction[NumberOfCatalyzedReactions]=inputReactions[i][j];
				NumberOfCatalyzedReactions++;
			}
		}

		//go over the outputReactions
	//	n=outputReactions[i].size();
		n=NumberOfOutputReactions[i];
		for(j=0; j<n; j++)
		{
			switch(tagList[outputReactions[i][j]])
			{
			case 1:
			//	normalReaction.push_back(outputReactions[i][j]);
				normalReaction[NumberOfNormalReactions]=outputReactions[i][j];
				NumberOfNormalReactions++;
				break;
			case 2:
			//	revisedReaction.push_back(outputReactions[i][j]);
				revisedReaction[NumberOfRevisedReactions]=outputReactions[i][j];
				NumberOfRevisedReactions++;
				break;
			case 3:
			//	pairProduct.push_back(outputReactions[i][j]);
				pairProduct[NumberOfPairProducts]=outputReactions[i][j];
				NumberOfPairProducts++;
				break;
			case 4:
			//	pairReactant.push_back(outputReactions[i][j]);
				pairReactant[NumberOfPairReactants]=outputReactions[i][j];
				NumberOfPairReactants++;
				break;
			case 6:
			//	oneDirectionReaction.push_back(outputReactions[i][j]);
				oneDirectionReaction[NumberOfOneDirectionReactions]=outputReactions[i][j];
				NumberOfOneDirectionReactions++;
				break;
			case 7:
			//	catalyzedReaction.push_back(outputReactions[i][j]);
				catalyzedReaction[NumberOfCatalyzedReactions]=outputReactions[i][j];
				NumberOfCatalyzedReactions++;
			}
		}
	}
	
	templateList double typeList::meanf(std::size_t i, double t, bool &ifpositive)
	{
		std::size_t j;
		double result=0, temp;
		std::vector<int> pair;
		std::vector<int>::iterator Iterator;
	/*	for(Iterator=normalReaction.begin(); Iterator!=normalReaction.end(); Iterator++)
			result+=currentPropensities[*Iterator]*t*stoi(*Iterator, i);*/
		for(j=0; j<NumberOfNormalReactions; j++)
			result+=currentPropensities[normalReaction[j]]*t*stoi(normalReaction[j], i);
	/*	for(Iterator=revisedReaction.begin(); Iterator!=revisedReaction.end(); Iterator++)
			result+=fmean(*Iterator, t)*stoi(*Iterator, i);*/
		for(j=0; j<NumberOfRevisedReactions; j++)
			result+=fmean(revisedReaction[j], t)*stoi(revisedReaction[j], i);
		savedmean=0;
	/*	for(Iterator=pairReactant.begin(); Iterator!=pairReactant.end(); Iterator++)
		{
			pair=responsiblepair[*Iterator];
			pairedparameter[pair[0]][pair[1]].getMean01(t);
			if(revisable_2nd(pair[1], pair[0])>=0)
			{
				result+=pairedparameter[pair[0]][pair[1]].mean*stoi(revisable_2nd(pair[1], pair[0]), i);
				savedmean+=pairedparameter[pair[0]][pair[1]].dmean*stoi(revisable_2nd(pair[1], pair[0]), i);
			}
			else
			{
				result-=pairedparameter[pair[0]][pair[1]].mean*stoi(revisable_2nd(pair[0], pair[1]), i);
				savedmean-=pairedparameter[pair[0]][pair[1]].dmean*stoi(revisable_2nd(pair[0], pair[1]), i);
			}
		}*/
		for(j=0; j<NumberOfPairReactants; j++)
		{
			pair=responsiblepair[pairReactant[j]];
			pairedparameter[pair[0]][pair[1]].getMean01(t);
			if(revisable_2nd(pair[1], pair[0])>=0)
			{
				result+=pairedparameter[pair[0]][pair[1]].mean*stoi(revisable_2nd(pair[1], pair[0]), i);
				savedmean+=pairedparameter[pair[0]][pair[1]].dmean*stoi(revisable_2nd(pair[1], pair[0]), i);
			}
			else
			{
				result-=pairedparameter[pair[0]][pair[1]].mean*stoi(revisable_2nd(pair[0], pair[1]), i);
				savedmean-=pairedparameter[pair[0]][pair[1]].dmean*stoi(revisable_2nd(pair[0], pair[1]), i);
			}
		}
	/*	for(Iterator=pairProduct.begin(); Iterator!=pairProduct.end(); Iterator++)
		{
			pair=responsiblepair[*Iterator];
			if(responsiblespeices[*Iterator]==pair[0])
				pairedparameter[pair[0]][pair[1]].getMean2(t);
			else
				pairedparameter[pair[0]][pair[1]].getMean3(t);
			result+=pairedparameter[pair[0]][pair[1]].mean*c_over_csum[*Iterator]*stoi(*Iterator, i);
			savedmean+=pairedparameter[pair[0]][pair[1]].dmean*c_over_csum[*Iterator]*stoi(*Iterator, i);
		}*/
		for(j=0; j<NumberOfPairProducts; j++)
		{
			pair=responsiblepair[pairProduct[j]];
			if(responsiblespeices[pairProduct[j]]==pair[0])
				pairedparameter[pair[0]][pair[1]].getMean2(t);
			else
				pairedparameter[pair[0]][pair[1]].getMean3(t);
			result+=pairedparameter[pair[0]][pair[1]].mean*c_over_csum[pairProduct[j]]*stoi(pairProduct[j], i);
			savedmean+=pairedparameter[pair[0]][pair[1]].dmean*c_over_csum[pairProduct[j]]*stoi(pairProduct[j], i);
		}
	/*	for(Iterator=oneDirectionReaction.begin(); Iterator!=oneDirectionReaction.end(); Iterator++)
		{
			pair=responsiblepair[*Iterator];
			pairedparameter[pair[0]][pair[1]].getMean(t, *Iterator);
			result+=pairedparameter[pair[0]][pair[1]].mean*stoi(*Iterator, i);
			savedmean+=pairedparameter[pair[0]][pair[1]].dmean*stoi(*Iterator, i);
		}*/
		for(j=0; j<NumberOfOneDirectionReactions; j++)
		{
			pair=responsiblepair[oneDirectionReaction[j]];
			pairedparameter[pair[0]][pair[1]].getMean(t, oneDirectionReaction[j]);
			result+=pairedparameter[pair[0]][pair[1]].mean*stoi(oneDirectionReaction[j], i);
			savedmean+=pairedparameter[pair[0]][pair[1]].dmean*stoi(oneDirectionReaction[j], i);
		}
	/*	for(Iterator=catalyzedReaction.begin(); Iterator!=catalyzedReaction.end(); Iterator++)
		{
			result+=fcatalyst(*Iterator, t, temp)*stoi(*Iterator, i);
			savedmean+=temp*stoi(*Iterator, i);
		}*/
		for(j=0; j<NumberOfCatalyzedReactions; j++)
		{
			result+=fcatalyst(catalyzedReaction[j], t, temp)*stoi(catalyzedReaction[j], i);
			savedmean+=temp*stoi(catalyzedReaction[j], i);
		}
		ifpositive=result>=0?true:false;
		return fabs(result);
	}
	
	templateList double typeList::meandf(std::size_t i, double t, bool ifpositive)
	{
		std::size_t j;
		double result=0;
		std::vector<int>::iterator Iterator;
	/*	for(Iterator=normalReaction.begin(); Iterator!=normalReaction.end(); Iterator++)
			result+=currentPropensities[*Iterator]*stoi(*Iterator, i);*/
		for(j=0; j<NumberOfNormalReactions; j++)
			result+=currentPropensities[normalReaction[j]]*stoi(normalReaction[j], i);
	/*	for(Iterator=revisedReaction.begin(); Iterator!=revisedReaction.end(); Iterator++)
			result+=dmean(*Iterator, t)*stoi(*Iterator, i);*/
		for(j=0; j<NumberOfRevisedReactions; j++)
			result+=dmean(revisedReaction[j], t)*stoi(revisedReaction[j], i);
		result+=savedmean;
		if(ifpositive)
			return result;
		else
			return -result;
	}
	
	templateList double typeList::varf(std::size_t i, double t)
	{
		std::size_t j;
		double result=0, temp;
		std::vector<int> pair;
		std::vector<int>::iterator Iterator;
	/*	for(Iterator=normalReaction.begin(); Iterator!=normalReaction.end(); Iterator++)
			result+=currentPropensities[*Iterator]*t*squarestoi(*Iterator, i);*/
		for(j=0; j<NumberOfNormalReactions; j++)
			result+=currentPropensities[normalReaction[j]]*t*squarestoi(normalReaction[j], i);
	/*	for(Iterator=revisedReaction.begin(); Iterator!=revisedReaction.end(); Iterator++)
			result+=fvar(*Iterator, t)*squarestoi(*Iterator, i);*/
		for(j=0; j<NumberOfRevisedReactions; j++)
			result+=fvar(revisedReaction[j], t)*squarestoi(revisedReaction[j], i);
		savedvar=0;
	/*	for(Iterator=pairReactant.begin(); Iterator!=pairReactant.end(); Iterator++)
		{
			pair=responsiblepair[*Iterator];
			pairedparameter[pair[0]][pair[1]].getVar01(t);
			result+=pairedparameter[pair[0]][pair[1]].var*squarestoi(*Iterator, i);
			savedvar+=pairedparameter[pair[0]][pair[1]].dvar*squarestoi(*Iterator, i);
		}*/
		for(j=0; j<NumberOfPairReactants; j++)
		{
			pair=responsiblepair[pairReactant[j]];
			pairedparameter[pair[0]][pair[1]].getVar01(t);
			result+=pairedparameter[pair[0]][pair[1]].var*squarestoi(pairReactant[j], i);
			savedvar+=pairedparameter[pair[0]][pair[1]].dvar*squarestoi(pairReactant[j], i);
		}
	/*	for(Iterator=pairProduct.begin(); Iterator!=pairProduct.end(); Iterator++)
		{
			pair=responsiblepair[*Iterator];
			if(responsiblespeices[*Iterator]==pair[0])
			{
				pairedparameter[pair[0]][pair[1]].getVar2(t);
				result+=pairedparameter[pair[0]][pair[1]].var*c_over_csum[*Iterator]*squarestoi(*Iterator, i);
				savedvar+=pairedparameter[pair[0]][pair[1]].dvar*c_over_csum[*Iterator]*squarestoi(*Iterator, i);
			}
			else
			{
				pairedparameter[pair[0]][pair[1]].getVar3(t);
				result+=pairedparameter[pair[0]][pair[1]].var*c_over_csum[*Iterator]*squarestoi(*Iterator, i);
				savedvar+=pairedparameter[pair[0]][pair[1]].dvar*c_over_csum[*Iterator]*squarestoi(*Iterator, i);
			}
		}*/
		for(j=0; j<NumberOfPairProducts; j++)
		{
			pair=responsiblepair[pairProduct[j]];
			if(responsiblespeices[pairProduct[j]]==pair[0])
			{
				pairedparameter[pair[0]][pair[1]].getVar2(t);
				result+=pairedparameter[pair[0]][pair[1]].var*c_over_csum[pairProduct[j]]*squarestoi(pairProduct[j], i);
				savedvar+=pairedparameter[pair[0]][pair[1]].dvar*c_over_csum[pairProduct[j]]*squarestoi(pairProduct[j], i);
			}
			else
			{
				pairedparameter[pair[0]][pair[1]].getVar3(t);
				result+=pairedparameter[pair[0]][pair[1]].var*c_over_csum[*Iterator]*squarestoi(*Iterator, i);
				savedvar+=pairedparameter[pair[0]][pair[1]].dvar*c_over_csum[*Iterator]*squarestoi(*Iterator, i);
			}
		}
	/*	for(Iterator=oneDirectionReaction.begin(); Iterator!=oneDirectionReaction.end(); Iterator++)
		{
			pair=responsiblepair[*Iterator];
			pairedparameter[pair[0]][pair[1]].getMean(t, *Iterator);
			result+=pairedparameter[pair[0]][pair[1]].var*squarestoi(*Iterator, i);
			savedvar+=pairedparameter[pair[0]][pair[1]].dvar*squarestoi(*Iterator, i);
		}*/
		for(j=0; j<NumberOfOneDirectionReactions; j++)
		{
			pair=responsiblepair[oneDirectionReaction[j]];
			pairedparameter[pair[0]][pair[1]].getMean(t, oneDirectionReaction[j]);
			result+=pairedparameter[pair[0]][pair[1]].var*squarestoi(oneDirectionReaction[j], i);
			savedvar+=pairedparameter[pair[0]][pair[1]].dvar*squarestoi(oneDirectionReaction[j], i);
		}
	/*	for(Iterator=catalyzedReaction.begin(); Iterator!=catalyzedReaction.end(); Iterator++)
		{
			result+=varcatalyst(*Iterator, t, temp)*squarestoi(*Iterator, i);
			savedvar+=temp*squarestoi(*Iterator, i);
		}*/
		for(j=0; j<NumberOfCatalyzedReactions; j++)
		{
			result+=varcatalyst(catalyzedReaction[j], t, temp)*squarestoi(catalyzedReaction[j], i);
			savedvar+=temp*squarestoi(catalyzedReaction[j], i);
		}
		return result;
	}
	
	templateList double typeList::vardf(std::size_t i, double t)
	{
		std::size_t j;
		double result=0;
		std::vector<int>::iterator Iterator;
	/*	for(Iterator=normalReaction.begin(); Iterator!=normalReaction.end(); Iterator++)
			result+=currentPropensities[*Iterator]*squarestoi(*Iterator, i);*/
		for(j=0; j<NumberOfNormalReactions; j++)
			result+=currentPropensities[normalReaction[j]]*squarestoi(normalReaction[j], i);
	/*	for(Iterator=revisedReaction.begin(); Iterator!=revisedReaction.end(); Iterator++)
			result+=dvar(*Iterator, t)*squarestoi(*Iterator, i);*/
		for(j=0; j<NumberOfRevisedReactions; j++)
			result+=dvar(revisedReaction[j], t)*squarestoi(revisedReaction[j], i);
		result+=savedvar;
		return result;
	}
	
	templateList void typeList::update_time_invariant()
	{
		std::size_t i;
		std::vector<int>::iterator vIterator;
		std::vector<size_t>::iterator sIterator;
		std::list<std::size_t>::iterator lIterator;
		for(lIterator=revisedSpecies.begin(); lIterator!=revisedSpecies.end(); lIterator++)
		{
			asum_over_csum[*lIterator]=a_sum[*lIterator]/c_sum[*lIterator];
			cx[*lIterator]=currentPopulation[*lIterator]*c_sum[*lIterator];
			/*	for(vIterator=torevisedreaction[*lIterator].begin(); vIterator!=torevisedreaction[*lIterator].end(); vIterator++)
			{
			c_over_csum[*vIterator]=c[*vIterator]/c_sum[*lIterator];
			c2[*vIterator]=2*c[*vIterator];
			}*/
		/*	for(sIterator=outputReactions[*lIterator].begin(); sIterator!=outputReactions[*lIterator].end(); sIterator++)
			{
				c_over_csum[*sIterator]=c[*sIterator]/c_sum[*lIterator];
				c2[*sIterator]=2*c[*sIterator];
			}*/
			for(i=0; i<NumberOfOutputReactions[*lIterator]; i++)
			{
				c_over_csum[outputReactions[i]]=c[outputReactions[i]]/c_sum[*lIterator];
				c2[outputReactions[i]]=2*c[outputReactions[i]];
			}
		}

		//for revisedPair
		std::list<std::vector<int> >::iterator pairIterator;
		PairedParameters *parameter;
		std::size_t i, j;
		for(pairIterator=revisedPair.begin(); pairIterator!=revisedPair.end(); pairIterator++)
		{
			parameter=&(pairedparameter[(*pairIterator)[0]][(*pairIterator)[1]]);
			for(j=0; j<2; j++)
			{
				/*	for(vIterator=torevisedreaction[(*pairIterator)[j]].begin(); vIterator!=torevisedreaction[(*pairIterator)[j]].end(); vIterator++)
				c_over_csum[*vIterator]=c[*vIterator]/c_sum[(*pairIterator)[j]];*/
				for(i=0; i<parameter->numberOfOutputReactions[j]; i++)
					c_over_csum[parameter->outputReactions[j][i]]=c[parameter->outputReactions[j][i]]/c_sum[(*pairIterator)[j]];
				/*	for(sIterator=outputReactions[(*pairIterator)[j]].begin(); sIterator!=outputReactions[(*pairIterator)[j]].end(); sIterator++)
				c_over_csum[*sIterator]=c[*sIterator]/c_sum[(*pairIterator)[j]];*/
			}
		}
	}
	
	templateList void typeList::update_time_dependent(double tau)
	{
		std::size_t i;
		std::vector<int>::iterator Iterator;
	/*	for(Iterator=revisedReaction.begin(); Iterator!=revisedReaction.end(); Iterator++)
		{
			exp_ct[*Iterator]=exp(-c_sum[responsiblespeices[*Iterator]]*tau);
			one_exp_ct[*Iterator]=1-exp_ct[*Iterator];
		}*/
		for(i=0; i<NumberOfRevisedReactions; i++)
		{
			exp_ct[revisedReaction[i]]=exp(-c_sum[responsiblespeices[revisedReaction[i]]]*tau);
			one_exp_ct[revisedReaction[i]]=1-exp_ct[revisedReaction[i]];
		}
	}
	
	templateList double typeList::solve_tau(std::size_t i, double tau)
	{
		double numerator, f, tmean, tvar;
		//get t for mean change.
		bool ifpositive;
		numerator=epsilon*currentPopulation[i]/g;//temp--should use G
		numerator=std::max(numerator,1.0);
		double temp;
		int iteration=0, maxiteration=100, stop=100;
		/*	if(stepsize[i]>=tau)
		{*/
		tmean=tau;
		temp=tmean;
		update_time_dependent(tmean);
		f=meanf(i, tmean, ifpositive)-numerator;
		if(f>0)
		{
			/*	if(tau>last_stepsize[i])
			{
			tmean=last_stepsize[i];
			update_time_dependent(tmean);
			f=meanf(i, tmean, ifpositive)-numerator;
			}*/
			while(fabs(f)>0.1 && iteration<=maxiteration && (f>0 || fabs(temp)/tmean>1e-2))
			{
				temp=f/meandf(i, tmean, ifpositive);
				if(tmean<temp)
					tmean/=2;
				else
					tmean-=temp;
				update_time_dependent(tmean);
				f=meanf(i, tmean, ifpositive)-numerator;
				iteration++;
			}
#ifdef PROFILE_SIMULATE_TAU_LEAPING_EXPLICIT_ADAPTIVE
			numberOfEquation++;
			numberOfIteration+=iteration;
#endif
			/*	if(fabs(f)>1e-6)
			{
			iteration=0;
			tmean=100;
			f=meanf(i, tmean, ifpositive)-numerator;
			while(fabs(f)>1e-6 && iteration<=maxiteration)
			{
			temp=f/meandf(i, tmean, ifpositive);
			if(tmean<temp)
			tmean/=2;
			else
			tmean-=temp;
			f=meanf(i, tmean, ifpositive)-numerator;
			iteration++;
			}
			}*/
			if(fabs(f)>0.1)
			{
				tmean=stepsize[i];
				temp=tmean;
				update_time_dependent(tmean);
				iteration=0;//restart
				f=meanf(i, tmean, ifpositive)-numerator;
				while(fabs(f)>1e-6 && iteration<=stop && (f>0 || fabs(temp)/tmean>1e-2))
				{
					temp=f/meandf(i, tmean, ifpositive);
					if(tmean<temp)
						tmean/=2;
					else
						tmean-=temp;
					update_time_dependent(tmean);
					f=meanf(i, tmean, ifpositive)-numerator;
					iteration++;
				}
#ifdef PROFILE_SIMULATE_TAU_LEAPING_EXPLICIT_ADAPTIVE
				numberOfIteration+=iteration;
#endif
				if(iteration>maxiteration)
				{
					std::cout<<"failed to converge, the original tau-leaping stepsize is used \n";
					tmean=stepsize[i];
				}
			}
		}
		//get t for var change.
		iteration=0;
		numerator*=numerator;
		tvar=tmean;
		temp=tmean;
		f=varf(i, tvar)-numerator;
		if(f>0)
		{
			while(fabs(f)>0.1 && iteration<=maxiteration && (f>0 || fabs(temp)/tmean>1e-2))
			{
				temp=f/vardf(i, tvar);
				if(tvar<temp)
					tvar/=2;
				else
					tvar-=temp;
				update_time_dependent(tvar);
				f=varf(i, tvar)-numerator;
				iteration++;
			}
#ifdef PROFILE_SIMULATE_TAU_LEAPING_EXPLICIT_ADAPTIVE
			numberOfEquation++;
			numberOfIteration+=iteration;
#endif
			if(iteration>maxiteration)
			{
				std::cout<<"failed to converge, the original tau-leaping stepsize is used \n";
				tvar=tmean;
			}
		}
		//	}
		/*	else
		{
		tmean=endtime;
		f=meanf(i, tmean, ifpositive)-numerator;
		while(f<0)
		{
		tmean*=2;
		f=meanf(i, tmean, ifpositive)-numerator;
		}

		while(fabs(f)>1e-6 && iteration<=maxiteration)
		{
		temp=f/meandf(i, tmean, ifpositive);
		if(tmean<temp)
		tmean/=2;
		else
		tmean-=temp;
		f=meanf(i, tmean, ifpositive)-numerator;
		iteration++;
		}
		if(fabs(f)>1e-6)
		{
		iteration=0;
		tmean=100;
		f=meanf(i, tmean, ifpositive)-numerator;
		while(fabs(f)>1e-6 && iteration<=maxiteration)
		{
		temp=f/meandf(i, tmean, ifpositive);
		if(tmean<temp)
		tmean/=2;
		else
		tmean-=temp;
		f=meanf(i, tmean, ifpositive)-numerator;
		iteration++;
		}
		}
		if(fabs(f)>1e-6)
		{
		std::cout<<"failed to converge, the original tau-leaping stepsize is used \n";
		tmean=stepsize[i];
		}

		//get t for var change.
		iteration=0;
		numerator*=numerator;
		tvar=tmean;
		f=varf(i, tvar)-numerator;
		if(f>0)
		{
		while(fabs(f)>1e-6 && iteration<=maxiteration)
		{
		temp=f/vardf(i, tvar);
		if(tvar<temp)
		tvar/=2;
		else
		tvar-=temp;
		f=varf(i, tvar)-numerator;
		iteration++;
		}
		if(fabs(f)>1e-6)
		{
		std::cout<<"failed to converge, the original tau-leaping stepsize is used \n";
		tvar=tmean;
		}
		}

		}*/
		return std::min(tmean, tvar);
		/*	stepsize[i]=std::min(tmean, tvar);
		return stepsize[i];*/
	}
	
	templateList double typeList::try_solve_tau(std::size_t i, double tau)
	{
		double numerator, f, tmean, tvar;
		//get t for mean change.
		bool ifpositive;
		numerator=epsilon*currentPopulation[i]/g;//temp--should use G
		numerator=std::max(numerator,1.0);
		double temp;
		int iteration=0, maxiteration=10, stop=100;
		tmean=tau;
		update_time_dependent(tmean);
		f=meanf(i, tmean, ifpositive)-numerator;
		while(fabs(f)>1e-6 && iteration<=maxiteration)
		{
			temp=f/meandf(i, tmean, ifpositive);
			if(tmean<temp)
				tmean/=2;
			else
				tmean-=temp;
			update_time_dependent(tmean);
			f=meanf(i, tmean, ifpositive)-numerator;
			iteration++;
		}
#ifdef PROFILE_SIMULATE_TAU_LEAPING_EXPLICIT_ADAPTIVE
		numberOfEquation++;
		numberOfIteration+=iteration;
#endif
		if(fabs(f)>1e-6)
			return 0;
		//get t for var change.
		iteration=0;
		numerator*=numerator;
		tvar=tmean;
		f=varf(i, tvar)-numerator;
		if(f>0)
		{
			while(fabs(f)>1e-6 && iteration<=maxiteration)
			{
				temp=f/vardf(i, tvar);
				if(tvar<temp)
					tvar/=10;
				else
					tvar-=temp;
				update_time_dependent(tvar);
				f=varf(i, tvar)-numerator;
				iteration++;
			}
#ifdef PROFILE_SIMULATE_TAU_LEAPING_EXPLICIT_ADAPTIVE
			numberOfEquation++;
			numberOfIteration+=iteration;
#endif
			if(fabs(f)>1e-6)
				return 0;
		}
		return std::min(tmean, tvar);
	}
	
	templateList double typeList::reevaluate_tau(double criticalStepsize, double endtime)
	{
		double tau, returned_tau;
		tau=criticalStepsize==-1?endtime:criticalStepsize;//different from iterate from noncriticalStepsize
		std::size_t n=mergedSpecies.size();
		std::list<std::size_t>::iterator Iterator;
		/*	for(Iterator=mergedSpecies.begin(); Iterator!=mergedSpecies.end(); Iterator++)
		{
		if(reevaluatedSpecies[*Iterator]==accumulator)
		tau=std::min(tau, solve_tau(*Iterator, tau));
		else
		tau=std::min(tau, stepsize[*Iterator]);
		}*/
		//if(lastIterationSpecies!=-1 && reevaluatedSpecies[lastIterationSpecies]==accumulator)
		//{
		//	form_formular(lastIterationSpecies);
		//	returned_tau=try_solve_tau(lastIterationSpecies, lastIterationStepsize);
		//	if(returned_tau!=0)//find the solution from the previous step
		//	{
		//		stepsize[lastIterationSpecies]=returned_tau;
		//		if(tau>returned_tau)
		//		{
		//			tau=returned_tau;
		//			lastIterationStepsize=tau;
		//		}
		//		else
		//			lastIterationSpecies=-1;
		//	}
		//	else
		//		lastIterationSpecies=-1;
		//	for(Iterator=mergedSpecies.begin(); Iterator!=mergedSpecies.end(); Iterator++)
		//	{
		//		if(reevaluatedSpecies[*Iterator]!=accumulator)
		//		{
		//			if(tau>stepsize[*Iterator])
		//			{
		//				tau=stepsize[*Iterator];
		//				lastIterationSpecies=-1;
		//			}
		//		}
		//		else
		//		{
		//			if(stepsize[*Iterator]<tau)
		//			{
		//				form_formular(*Iterator);
		//				returned_tau=solve_tau(*Iterator, tau);
		//				if(tau>returned_tau)
		//				{
		//					tau=returned_tau;
		//					lastIterationSpecies=*Iterator;
		//					lastIterationStepsize=tau;
		//				}
		//			}
		//		}
		//	}
		//}
		//else //the stepsize in last step is not from the reevaluatedSpecies species
		//{
		for(Iterator=mergedSpecies.begin(); Iterator!=mergedSpecies.end();)
		{
			if(reevaluatedSpecies[*Iterator]!=accumulator)
			{
				tau=std::min(tau, stepsize[*Iterator]);
				Iterator=mergedSpecies.erase(Iterator);
			}
			else
				Iterator++;
		}
		for(Iterator=mergedSpecies.begin(); Iterator!=mergedSpecies.end(); Iterator++)
		{
			//	if(stepsize[*Iterator]<tau)//only increase stepsize?
			//	{
			form_formular(*Iterator);
			returned_tau=solve_tau(*Iterator, tau);
			if(tau>returned_tau)
			{
				tau=returned_tau;
				//				lastIterationSpecies=*Iterator;
				//				lastIterationStepsize=tau;
			}
			//	}
		}
		//}
		/*	for(Iterator=criticalSpecies.begin(); Iterator!=criticalSpecies.end(); Iterator++)
		{
		if(reevaluatedSpecies[*Iterator]==accumulator)
		tau=std::min(tau, solve_tau(*Iterator));
		else
		tau=std::min(tau, stepsize[*Iterator]);
		}*/
		return tau;
	}

#ifdef PRINT_RELAXED
	templateList void typeList::reevaluate_tau(double criticalStepsize, double endtime, int &confineSpecies)
	{
		double tau, temp;
		tau=criticalStepsize==-1?endtime:criticalStepsize;//different from iterate from noncriticalStepsize
		confineSpecies=-1;
		std::size_t n=mergedSpecies.size();
		std::list<std::size_t>::iterator Iterator;
		for(Iterator=mergedSpecies.begin(); Iterator!=mergedSpecies.end(); Iterator++)
		{
			if(reevaluatedSpecies[*Iterator]==accumulator)
			{
				temp=solve_tau(*Iterator, tau);
				if(tau>temp)
				{
					tau=temp;
					confineSpecies=*Iterator;
				}
			}
			else
			{
				if(tau>stepsize[*Iterator])
				{
					tau=stepsize[*Iterator];
					confineSpecies=*Iterator;
				}
			}
		}
		/*	for(Iterator=criticalSpecies.begin(); Iterator!=criticalSpecies.end(); Iterator++)
		{
		if(reevaluatedSpecies[*Iterator]==accumulator)
		tau=std::min(tau, solve_tau(*Iterator));
		else
		tau=std::min(tau, stepsize[*Iterator]);
		}*/
		return tau;
	}
#endif
	templateList template<typename IntervalOutputType> void typeList::simulate(std::size_t realizations, double startTime, double endTime, IntervalOutputType& output, bool doValidate=true) 
	{

			//validate();
			if (doValidate) {
				if (!SSA::validate(startTime,endTime)) {
					std::cerr << "StochKit ERROR (TauLeapingExplicitAdaptive::simulate): validate() failed, simulation aborted\n";
					exit(1);
				}
			}

			if (!output.initialize(realizations,startTime,endTime,initialPopulation)) {
				std::cerr << "StochKit ERROR (): initialization of output object failed, simulation aborted\n";
				exit(1);
			}

			//delete the product for tau selection and
			//create the reversible pair matrix
			TauLeaping::prepare();

			creatPairMatrix();

			std::vector<double> outputTimes=output.getOutputTimes();
			std::size_t totalIntervals=outputTimes.size();
			std::size_t currentInterval;
			std::size_t failedLeaps;
			double tau;
			double nextTime;
			double reactionsLastLeap;
			std::size_t ssaStepsTaken;
			double criticalStepsize;//for neg
			double noncriticalStepsize;//for neg
			bool runCritical;//for neg
			int criticalIndex;//for neg

#ifdef _MYDEBUG
			double ave_stepsize;
			double old_stepsize;
			double steps;
#endif

#ifdef PRINT_RELAXED
			std::ofstream outRelax("relaxed_species.txt");
			std::list<std::size_t>::iterator printrelax_i;
			int printrelax_index;
#endif

			//for testing
#ifdef PROFILE_SIMULATE_TAU_LEAPING_EXPLICIT_ADAPTIVE
			double timeInSSA=0.0;
			double ssaStartTime=0.0;
			std::size_t totalSSAStepsTaken=0;
			std::size_t totalLeapsTaken=0;
			double totalReactionsDuringLeaps=0;
			std::size_t totalImprovedLeapsTaken=0;
			double showtime;//debug
			//	std::ofstream out("result.txt");//debug out
			//	std::ofstream traj("trajectory.txt");//debug traj
			//	traj<<"0	"<<initialPopulation[0]<<"	"<<initialPopulation[1]<<"	"<<initialPopulation[2]<<"\n";//debug traj
			bool leaping=false;//debug example
			int index[2]={0,0};//debug out
			//	int accumulator=1;//debug traj
			//	int i;//debug traj
			//	std::list<int>::iterator Iterator;//debug critical
			double secondtau=0;//debug out
#endif
#ifdef EACH_TIME
			LARGE_INTEGER begin, middle, end, freq;
			double myiteration=0, myequation=0, myleaps=0, myssa=0, mytimeinSSA=0, middletime, elapsedTime;
			bool go;
#endif
#ifdef OUTPUT
			std::ofstream out("species.txt");
			std::list<std::size_t>::iterator singleIterator;
			std::list<std::vector<int> >::iterator pairIterator;
#endif
			for (std::size_t currentRealization=0; currentRealization!=realizations; ++currentRealization) {
#ifdef EACH_TIME
				randomGenerator.Seed(0);
				QueryPerformanceFrequency(&freq);
				QueryPerformanceCounter(&begin);
				go=true;
#endif
#ifdef _MYDEBUG
				//	std::cout<<currentRealization<<"\n";
				ave_stepsize=0;
				old_stepsize=0;
				steps=0;
#endif

				initialize(startTime, endTime);
#ifdef PROFILE_SIMULATE_TAU_LEAPING_EXPLICIT_ADAPTIVE//debug print
				showtime=startTime;
				//	std::cout<<showtime<<"\n";
#endif
				currentInterval=0;
				failedLeaps=0;
				ssaStepsTaken=0;
				reactionsLastLeap=std::numeric_limits<double>::max();

				while (currentTime<endTime) {
					while ((currentTime>=outputTimes[currentInterval]) && currentInterval<totalIntervals) {
						output.record(currentRealization, currentInterval, currentPopulation);
						++currentInterval;
					}
					/*
					#ifdef EACH_TIME
					if(go)
					{
					if(currentTime>16)
					{
					QueryPerformanceCounter(&middle);
					middletime=(double)(middle.QuadPart-begin.QuadPart)/freq.QuadPart;
					go=false;
					}
					}
					#endif
					*/
					//	std::cout<<currentTime<<"\n";
#ifdef PRINT_RELAXED
					//	std::cout<<currentTime<<"\n";
#endif

					if (reactionsLastLeap<threshold) {
						//do SSA
						failedLeaps=0;//not taking leaps
						if (ssaStepsTaken<SSASteps) {
							if (ssaStepsTaken==0) {
								currentTime+=SSA::selectStepSize();
								++ssaStepsTaken;
#ifdef PROFILE_SIMULATE_TAU_LEAPING_EXPLICIT_ADAPTIVE
								ssaStartTime=currentTime;
								if(currentTime-showtime>1)//debug print
								{
									showtime=currentTime;
									//	std::cout<<showtime<<"	SSA Step\n";
								}
								/* 	if(currentTime>accumulator)//debug traj
								{
								accumulator+=1;
								traj<<currentTime<<"	"<<currentPopulation[0]<<"	"<<currentPopulation[1]<<"	"<<currentPopulation[2]<<"\n";//debug traj
								}*/
#endif
							}
							else {
								SSA::fireReaction(SSA::selectReaction());
								currentTime+=SSA::selectStepSize();
								++ssaStepsTaken;
#ifdef PROFILE_SIMULATE_TAU_LEAPING_EXPLICIT_ADAPTIVE
								++totalSSAStepsTaken;
								if(currentTime-showtime>1)//debug print
								{
									showtime=currentTime;
									//	std::cout<<showtime<<"	SSA Step\n";
								}
								/* 	if(currentTime>accumulator)//debug traj
								{
								accumulator+=1;
								traj<<currentTime<<"	"<<currentPopulation[0]<<"	"<<currentPopulation[1]<<"	"<<currentPopulation[2]<<"\n";//debug traj
								}*/
#endif
							}
						}
						else { // we've taken enough SSA steps, try to go back to tau-leaping
							SSA::fireReaction(SSA::selectReaction());
							reactionsLastLeap=std::numeric_limits<double>::max();
							ssaStepsTaken=0;
#ifdef PROFILE_SIMULATE_TAU_LEAPING_EXPLICIT_ADAPTIVE
							timeInSSA+=(currentTime-ssaStartTime);
#endif
						}
					}
					else
					{ 
						//do tau-leaping
						accumulator++;//for revise
						//	std::cout<<accumulator<<"\n";
						updateTagLists();//for neg
						selectTau(noncriticalStepsize, criticalStepsize);//for neg
#ifdef USE_OLDSTEP
						tau=endTime;
						runCritical=false;
						if(!noncriticalSpecies.empty())
							tau=std::min(stepsize[*(noncriticalSpecies.begin())], tau);
						if(criticalStepsize!=-1)
						{
							if(tau<criticalStepsize)
								runCritical=false;
							else
							{
								runCritical=true;
								tau=criticalStepsize;
							}
						}
#else
						getParameters();
						update_time_invariant();
						if(noncriticalStepsize<criticalStepsize||criticalStepsize==-1)//for neg
						{
#ifdef PRINT_RELAXED
							outRelax<<"current_time: "<<currentTime<<"	relaxed_species(index&stepsize):	";
							for(printrelax_i=revisedSpecies.begin(); printrelax_i!=revisedSpecies.end(); printrelax_i++)
								outRelax<<*printrelax_i<<" & "<<stepsize[*printrelax_i]<<"	";
#endif
							getReevaluated();

#ifdef _MYDEBUG
							old_stepsize+=noncriticalStepsize;
#endif

#ifndef PRINT_RELAXED
							tau=reevaluate_tau(criticalStepsize, endTime);
#else
							tau=reevaluate_tau(criticalStepsize, endTime, printrelax_index);
#endif

#ifdef _MYDEBUG
							//	std::cout<<noncriticalStepsize<<"\n";
							ave_stepsize+=tau;
							steps++;
#endif
#ifdef PRINT_RELAXED
							outRelax<<"relaxed_stepsize(index&stepsize):	"<<printrelax_index<<" & "<<tau<<"\n";
#endif
							if(tau<criticalStepsize||criticalStepsize==-1)
								runCritical=false;
							else
								runCritical=true;
						}
						else
						{
							tau=criticalStepsize;
							runCritical=true;
						}
#endif
						if (currentTime+tau>outputTimes[currentInterval]) {
							//don't leap past next output time
							tau=outputTimes[currentInterval]-currentTime;
							nextTime=outputTimes[currentInterval];
							runCritical=false;//for neg
						}	
						else {
							nextTime=currentTime+tau;
						}
						criticalIndex=selectReactions(tau, runCritical);
						reactionsLastLeap=norm_1(previousReactionCounts);
						if(runCritical)//for neg
							reactionsLastLeap++;
#ifdef PROFILE_SIMULATE_TAU_LEAPING_EXPLICIT_ADAPTIVE
						totalReactionsDuringLeaps+=reactionsLastLeap;
#endif
						if (TauLeaping::fireReactions(criticalIndex)) {
							currentTime=nextTime;
							failedLeaps=0;
							//	last_stepsize=stepsize;
#ifdef PROFILE_SIMULATE_TAU_LEAPING_EXPLICIT_ADAPTIVE
							++totalLeapsTaken;
#endif
						}
						else {
							++failedLeaps;
							if (failedLeaps==3) {
								std::cout << "StochKit WARNING (TauLeapingExplicitAdaptive::simulate): rejected three or more consecutive leaps, consider reducing epsilon.\n";
							}
							if (failedLeaps==10) {
								std::cerr << "StochKit ERROR (TauLeapingExplicitAdaptive::simulate): rejected ten consecutive leaps, terminating simulation.\n";
								exit(1);
							}
						}
					}
#ifdef OUTPUT
					out<<currentTime<<" ";
					for(singleIterator=revisedSpecies.begin(); singleIterator!=revisedSpecies.end(); singleIterator++)
						out<<*singleIterator<<" ";
					for(pairIterator=revisedPair.begin(); pairIterator!=revisedPair.end(); pairIterator++)
						out<<(*pairIterator)[0]<<"-"<<(*pairIterator)[1]<<" ";
					out<<std::endl;
#endif
				}

				//modified for visual studio
				while(currentTime>=outputTimes[currentInterval])
				{
					output.record(currentRealization, currentInterval, currentPopulation);
					++currentInterval;
					if(currentInterval>=totalIntervals)
						break;
				}
				/*	while ((currentTime>=outputTimes[currentInterval]) && currentInterval<totalIntervals) {
				output.record(currentRealization, currentInterval, currentPopulation);
				++currentInterval;
				}*/
				// Add the last ssa session to the timeInSSA
#ifdef EACH_TIME
				QueryPerformanceCounter(&end);
				elapsedTime=(double)(end.QuadPart-begin.QuadPart)/freq.QuadPart;
				std::cout<<"the simulation takes "<<elapsedTime<<" seconds. Iteration number: "<<numberOfIteration-myiteration<<" average per equation: "<<(numberOfIteration-myiteration)/(numberOfEquation-myequation)<<" leaps: "<<totalLeapsTaken-myleaps<<" ssa steps: "<<totalSSAStepsTaken-myssa<<" time in SSA: "<<timeInSSA-mytimeinSSA<<std::endl;
				myiteration=numberOfIteration;
				myequation=numberOfEquation;
				myleaps=totalLeapsTaken;
				myssa=totalSSAStepsTaken;
				mytimeinSSA=timeInSSA;
#endif
#ifdef PROFILE_SIMULATE_TAU_LEAPING_EXPLICIT_ADAPTIVE 
				if (ssaStepsTaken!=0)
					timeInSSA+=(endTime-ssaStartTime);
#endif

#ifdef _MYDEBUG
				std::cout<<"	"<<ave_stepsize/steps<<"	"<<old_stepsize/steps<<"	"<<steps<<"\n";
#endif

			}
#ifdef OUTPUT
			out.close();
#endif
#ifdef PROFILE_SIMULATE_TAU_LEAPING_EXPLICIT_ADAPTIVE
			std::cout << "StochKit MESSAGE (TauLeapingExplicitAdaptive::simulate) simulate() spent "<<timeInSSA<< " out of " << (endTime-startTime)*(double)realizations << " in SSA\n";
			std::cout << "and took "<<totalImprovedLeapsTaken<<" improved leaps and "<<totalLeapsTaken<<" leaps totaling "<<totalReactionsDuringLeaps<<" reactions and "<<totalSSAStepsTaken<<" SSA steps\n";
			std::cout << "Total number of iterations: "<<numberOfIteration<< ". Average iterations for each equation: "<<numberOfIteration/numberOfEquation<<"\n";
			/* 	for(i=0;i<NumberOfReactions;i++)//debug
			std::cout<<tagList[i]<<"	";
			std::cout<<"\n";*/
#endif
#ifndef DEBUG
			freeMemory();
#endif
	}//end simulate
#ifndef DEBUG
	templateList void typeList::freeMemory()
	{
		int i;
		myDelete[] (inputReactions, NumberOfSpecies);
		delete[] NumberOfInputReactions
		myDelete[] (OutputReactions, NumberOfSpecies);
		delete[] NumberOfOutputReactions;

		delete[] influencedSpeciesSize;
		delete[] revisable;
		delete[] stepsize;
		delete[] responsiblespeices;
		myDelete(responsiblepair, NumberOfReactions);
		delete[] c;
		delete[] c_sum;
		delete[] a_sum;
		delete[] c_in;
		delete[] c_out;
		for(i=0; i<NumberOfSpecies; i++)
			myDelete(revisable_2nd_reaction[i], NumberOfSpecies);
		delete[] revisable_2nd_reaction;
		myDelete(revisable_catalyze_reaction, NumberOfSpecies);
		delete[] normalReaction;
		delete[] revisedReaction;
		delete[] pairProduct;
		delete[] pairReactant;
		delete[] oneDirectionReaction;
		delete[] catalyzedReaction;
		myDelete(reactantList, NumberOfReactions);
		delete[] rateList;
		delete[] c_over_csum;
		delete[] asum_over_csum;
		delete[] cx;
		delete[] exp_ct;
		delete[] one_exp_ct;
		delete[] c2;
		myDelete(pairedparameter, NumberOfSpecies);
	}
#endif
}