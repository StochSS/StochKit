/******************************************************************************
 *  FILE:    Random.h                                                         *
 ******************************************************************************/

#ifndef STOCHKIT_RANDOM_H
#define STOCHKIT_RANDOM_H

#include <iostream>
#include <vector>
#include <limits>
#include <ctime>
#include <cmath>
#include <numeric>
#include <boost/random.hpp>
#include <boost/numeric/ublas/vector.hpp>

namespace STOCHKIT{
	class RandomGenerator
	{
		private:

			boost::mt19937 generatorUniform;

		public:
			// constructor
			RandomGenerator(){}; 

			//! allow compiler-generated copy constructor
			//! allow compiler-generated assignment operator

			// without specified seed, use current time
			void Seed()
			{
				generatorUniform.seed((boost::uint32_t) time(NULL));
			}

			// with specified seed RanSeed
			void Seed(boost::uint32_t RanSeed)
			{
				generatorUniform.seed(RanSeed);
			}

			// discrete uniform random number generator
			double DiscreteUniform()
			{
				return generatorUniform();
			}

			// continuous uniform random number in (a,b)
			double ContinuousOpen(double a, double b)
			{
				boost::uniform_01<> ContinuousZeroOneOpen;
				boost::variate_generator<boost::mt19937&, boost::uniform_01<> > generatorContinuousZeroOneOpen(generatorUniform, ContinuousZeroOneOpen);
				return (generatorContinuousZeroOneOpen()*(b-a) + a);
			}

			// exponential random number
			double Exponential(double mean)
			{
				double lambda;
				if(mean==0)
				{
					lambda=std::numeric_limits<double>::max();
					std::cout<<"0 mean is given to the exponential random number generator"<<"\n";
				}
				else if(mean>=std::numeric_limits<double>::max())
				{
					return mean;
				}
				else
					lambda=1/mean;
				boost::exponential_distribution<> ExponentialDistribution(lambda);
				boost::variate_generator<boost::mt19937&, boost::exponential_distribution<> > generatorExponential(generatorUniform, ExponentialDistribution);
				return generatorExponential();
			}

			// normal random number
			double Normal(double mean, double var)
			{
				boost::normal_distribution<> NormalDistribution(mean, sqrt(var));
				boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > generatorNormal(generatorUniform, NormalDistribution);
				return generatorNormal();
			}

			// poisson random number
			double Poisson(double mean)
			{
				double result;
				if(mean==0)
					return 0;
				else if(mean<=700)
				{
					boost::poisson_distribution<> PoissonDistribution(mean);
					boost::variate_generator<boost::mt19937&, boost::poisson_distribution<> > generatorPoisson(generatorUniform, PoissonDistribution);
					return generatorPoisson();
				}
				else
				{
					result=Normal(mean, mean);
					return floor(result+0.5);
				}
			}

			// binomial random number
			double Binomial(double n, double p)
			{
				double result, mean;

				p=std::min(1.0, p);
				p=std::max(0.0, p);

				if(n<=20)
				{
					boost::binomial_distribution<> BinomialDistribution((int)n, p);
					boost::variate_generator<boost::mt19937&, boost::binomial_distribution<> > generatorBinomial(generatorUniform, BinomialDistribution);
					return generatorBinomial();
				}
				else if(n<100)
				{
					mean=n*p;
					if(p<=0.05)
						result=std::min(n, Poisson(mean));
					else if(p>=0.95)
						result=n-std::min(n, Poisson(n-mean));
					else
						result=std::min(n, Normal(mean, mean*(1-p)));
					return floor(result+0.5);
				}
				else
				{
					mean=n*p;
					if(mean<=10)
						result=std::min(n, Poisson(mean));
					else if(n-mean<=10)
						result=n-std::min(n, Poisson(n-mean));
					else
						result=std::min(n, Normal(mean, mean*(1-p)));
					return floor(result+0.5);
				}
			}

			// multinomial random number
			template<typename vectorType>
			std::vector<double> Multinomial(double n, vectorType &p)
			{
				int i, m=p.size();
				std::vector<double> result(m, 0.0);
				if(n==0)
					return result;
#ifdef _DEBUG
				for(i=0; i<m; i++)
				{
					if(p[i]<0)
					{
						if(p[i]<-1e-8)
						{
							std::cerr<<"Multinomial error: negative probability \n";
							exit(1);
						}
						else
							p[i]=0;
					}
				}
				/*	if(p_sum<=0)
				{
				std::cerr<<"Multinomial error: zero probability sum \n";
				exit(1);
				}*/
#endif
				double p_sum=accumulate(p.begin(), p.end(), 0.0);
				i=0;
				while(n>0)
				{
					result[i]=Binomial(n, p[i]/p_sum);
					n-=result[i];
					p_sum-=p[i];
					i++;
				}
				return result;
			}

			
			template<typename vectorType>
			std::vector<double> Multinomial(double n, vectorType &p, std::size_t sizeOfp)
			{
				int i;
				std::vector<double> result(sizeOfp, 0.0);
				if(n==0)
					return result;
#ifdef _DEBUG
				for(i=0; i<sizeOfp; i++)
				{
					if(p[i]<0)
					{
						if(p[i]<-1e-8)
						{
							std::cerr<<"Multinomial error: negative probability \n";
							exit(1);
						}
						else
							p[i]=0;
					}
				}
				/*	if(p_sum<=0)
				{
				std::cerr<<"Multinomial error: zero probability sum \n";
				exit(1);
				}*/
#endif
				double p_sum=std::accumulate(p, p+sizeOfp, 0.0);
				i=0;
				while(n>0)
				{
					result[i]=Binomial(n, p[i]/p_sum);
					n-=result[i];
					p_sum-=p[i];
					i++;
				}
				return result;
			}
			
			template<typename vectorType>
			void Multinomial(double n, vectorType &p, std::vector<double> &result)
			{
				int i, m=p.size();
				if(n==0)
				{
					for(i=0; i<m; i++)
						result[i]=0;
				}
				else
				{
#ifdef _DEBUG
					for(i=0; i<m; i++)
					{
						if(p[i]<0)
						{
							if(p[i]<-1e-8)
							{
								std::cerr<<"Multinomial error: negative probability \n";
								exit(1);
							}
							else
								p[i]=0;
						}
					}
					/*	if(p_sum<=0)
					{
					std::cerr<<"Multinomial error: zero probability sum \n";
					exit(1);
					}*/
#endif
					double p_sum=accumulate(p.begin(), p.end(), 0.0);
					i=0;
					while(n>0)
					{
						result[i]=Binomial(n, p[i]/p_sum);
						n-=result[i];
						p_sum-=p[i];
						i++;
					}
					for(;i<m;i++)
						result[i]=0;
				}
			}
			
			template<typename vectorType>
			void Multinomial(double n, vectorType &p, std::size_t sizeOfp, std::vector<double> &result)
			{
				int i;
				if(n==0)
				{
					for(i=0; i<sizeOfp; i++)
						result[i]=0;
				}
				else
				{
#ifdef _DEBUG
					for(i=0; i<sizeOfp; i++)
					{
						if(p[i]<0)
						{
							if(p[i]<-1e-8)
							{
								std::cerr<<"Multinomial error: negative probability \n";
								exit(1);
							}
							else
								p[i]=0;
						}
					}
					/*	if(p_sum<=0)
					{
					std::cerr<<"Multinomial error: zero probability sum \n";
					exit(1);
					}*/
#endif
					double p_sum=std::accumulate(p, p+sizeOfp, 0.0);
					i=0;
					while(n>0)
					{
						result[i]=Binomial(n, p[i]/p_sum);
						n-=result[i];
						p_sum-=p[i];
						i++;
					}
					for(;i<sizeOfp;i++)
						result[i]=0;
				}
			}

			//ordered Multinomial
			template<typename vectorType, typename matrixType>
			boost::numeric::ublas::vector<double> Multinomial(std::vector<std::vector<std::size_t> > &order, vectorType population, matrixType &p, matrixType &sample)
			{
				std::size_t i, j, k, m=population.size();
				double p_sum, n;
				boost::numeric::ublas::vector<double> result(m);
				result.clear();
				for(i=0; i<m; i++)
				{
					p_sum=sum(row(p, i));
					n=population[i];
					j=0;
					while(n>0)
					{
						sample(i, order[i][j])=Binomial(n, p(i, order[i][j])/p_sum);
						result[order[i][j]]+=sample(i, order[i][j]);
						n-=sample(i, order[i][j]);
						if(p_sum>1e-8)
							p_sum-=p(i, order[i][j]);
						else
						{
							p_sum=0;
							for(k=j+1; k<m; k++)
								p_sum+=p(i, order[i][k]);
						}
						j++;
					}
				}
				return result;
			}	

			//ordered Multinomial, don't record exact samples
				//create partial sum vector
			template<typename vectorType, typename matrixType>
			boost::numeric::ublas::vector<double> Multinomial(/*std::vector<std::vector<std::size_t> > &order*/ std::size_t **order, vectorType population, matrixType &p)
			{
				std::size_t i, j, k, m=population.size();
				double p_sum, n, sample, *partialSum;
				partialSum=new double[m];
				boost::numeric::ublas::vector<double> result(m);
				result.clear();
				for(i=0; i<m; i++)
				{
					partialSum[m-1]=p(i, order[i][m-1]);
					for(j=m-2;j>=0;j--)
						partialSum[j]=partialSum[j+1]+p(i, order[i][j]);
				//	p_sum=sum(row(p, i));
					n=population[i];
					j=0;
					while(n>0)
					{
					/*	if(j==m)
						{
							result[order[i][j-1]]+=n;
							n=0;
						}
						else
						{*/
						//	sample=Binomial(n, p(i, order[i][j])/p_sum);
							sample=Binomial(n, p(i, order[i][j])/partialSum[j]);
							result[order[i][j]]+=sample;
							n-=sample;
						/*	if(p_sum>1e-8)
								p_sum-=p(i, order[i][j]);
							else
							{
								p_sum=0;
								for(k=j+1; k<m; k++)
									p_sum+=p(i, order[i][k]);
							}*/
							j++;
					//	}
					}
				}
				delete[] partialSum;
				return result;
			}
			
				//pass partial sum, single cell sample
			template<typename vectorType>
#ifdef DEBUG
			boost::numeric::ublas::vector<double> Multinomial(std::vector<std::size_t> &order, double population, vectorType &p, std::vector<double> &partialSum, boost::numeric::ublas::vector<double> &result)
#else
			boost::numeric::ublas::vector<double> Multinomial(std::size_t order[], double population, vectorType &p, double partialSum[], boost::numeric::ublas::vector<double> &result)			
#endif
			{
				std::size_t j, m=result.size();
				partialSum[m-1]=p(order[m-1]);
				for(j=m-1;j>0;j--)
					partialSum[j-1]=partialSum[j]+p(order[j-1]);
				//now j=0;
				result.clear();
				while(population>0)
				{
					result[order[j]]=Binomial(population, p(order[j])/partialSum[j]);
					population-=result[order[j]];
					j++;
				}
				return result;
			}

				//pass partial sum, all cells sample
			template<typename vectorType, typename matrixType>
#ifdef DEBUG
			boost::numeric::ublas::vector<double> Multinomial(std::vector<std::vector<std::size_t> > &order, vectorType &population, matrixType &p, std::vector<double> &partialSum, boost::numeric::ublas::vector<double> &result)
#else
			boost::numeric::ublas::vector<double> Multinomial(std::size_t **order, vectorType &population, matrixType &p, double partialSum[], boost::numeric::ublas::vector<double> &result)
#endif
			{
				std::size_t i, j, m=population.size();
				double n, sample;
				result.clear();
				for(i=0; i<m; i++)
				{
					if(population[i]>0)
					{
						partialSum[m-1]=p(i, order[i][m-1]);
						for(j=m-1;j>0;j--)
							partialSum[j-1]=partialSum[j]+p(i, order[i][j-1]);
						n=population[i];
						//now j=0;
						while(n>0)
						{
							sample=Binomial(n, p(i, order[i][j])/partialSum[j]);
							result[order[i][j]]+=sample;
							n-=sample;
							j++;
						}
					}
				}
				return result;
			}
			
				//pass partial sum, all cells sample, record number of binomial samplings in each multinomial sample
			template<typename vectorType, typename matrixType>
#ifdef DEBUG
			boost::numeric::ublas::vector<double> Multinomial(std::vector<std::vector<std::size_t> > &order, vectorType &population, matrixType &p, std::vector<double> &partialSum, boost::numeric::ublas::vector<double> &result, double &numberOfBinomialSamplings, double &numberOfMultinomialSamplings)
#else
			boost::numeric::ublas::vector<double> Multinomial(std::size_t **order, vectorType &population, matrixType &p, double partialSum[], boost::numeric::ublas::vector<double> &result, double &numberOfBinomialSamplings, double &numberOfMultinomialSamplings)
#endif
			{
				std::size_t i, j, m=population.size();
				double n, sample;
				result.clear();
				for(i=0; i<m; i++)
				{
					if(population[i]>0)
					{
						partialSum[m-1]=p(i, order[i][m-1]);
						for(j=m-1;j>0;j--)
							partialSum[j-1]=partialSum[j]+p(i, order[i][j-1]);
						n=population[i];
						//now j=0;
						while(n>0)
						{
							sample=Binomial(n, p(i, order[i][j])/partialSum[j]);
							result[order[i][j]]+=sample;
							n-=sample;
							j++;
						}
						numberOfBinomialSamplings+=j;
						numberOfMultinomialSamplings++;
					}
				}
				return result;
			}

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
#ifndef BOOST_RANDOM_NO_STREAM_OPERATORS
			// output random generator state to an ostream
			template<class CharT, class Traits>
			friend std::basic_ostream<CharT,Traits>&
				operator<<(std::basic_ostream<CharT,Traits>& os, const RandomGenerator& randomGenerator);

			// input random generator state from an istream
			// WARNING: no check to make sure this works, make sure the istream only contains relevant information, use at your own risk
			template<class CharT, class Traits>
			friend std::basic_istream<CharT,Traits>&
				operator>>(std::basic_istream<CharT,Traits>& is, RandomGenerator& randomGenerator);
#endif
#endif
	};//end class
}//end namespace

#endif
