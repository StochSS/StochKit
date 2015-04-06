#ifndef MYFUN_H
#define MYFUN_H
#include<vector>

#ifdef DEBUG
		template<typename valType>
		inline void myNew(std::vector<valType> &v, std::size_t n)
		{
			try{
				v.resize(n);
			}
			catch (std::bad_alloc)
			{
				assert(0);
				std::cerr << "bad_alloc caught: size=" <<n<< '\n';;
				exit(1);
			}
		}
		
		template<typename valType>
		inline void myNew(std::vector<valType> &v, std::size_t n, const valType &val)
		{
			try{
				v.resize(n, val);
			}
			catch (std::bad_alloc)
			{
				assert(0);
				std::cerr << "bad_alloc caught: size=" <<n<< '\n';;
				exit(1);
			}
		}

		template<typename valType, typename vectorType>
		inline void myNew(std::vector<valType> &destiny, vectorType &source)
		{
			try{
				destiny.insert(destiny.end(), source.begin(), source.end());
			}
			catch (std::bad_alloc)
			{
				assert(0);
				std::cerr << "bad_alloc caught: size=" <<source.size()<< '\n';;
				exit(1);
			}
		}

		template<typename valType>
		inline void myNew(std::vector<std::vector<valType> >&v, std::size_t m, std::size_t n)
		{
			try{
				v.resize(m, std::vector<valType>(n) );
			}
			catch (std::bad_alloc)
			{
				assert(0);
				std::cerr << "bad_alloc caught: size=" <<m<<"*"<<n<<"="<<m*n<< '\n';
				exit(1);
			}
		}
		
		template<typename valType>
		inline void myNew(std::vector<std::vector<valType> >&v, std::size_t m, std::size_t n, const valType &val)
		{
			try{
				v.resize(m, std::vector<valType>(n, val) );
			}
			catch (std::bad_alloc)
			{
				assert(0);
				std::cerr << "bad_alloc caught: size=" <<m<<"*"<<n<<"="<<m*n<< '\n';
				exit(1);
			}
		}
		
		template<typename valType, typename MatrixType>
		inline void myNew(std::vector<std::vector<valType> >&destiny, MatrixType &source)
		{
			int n=source.size();
			try{
				destiny.resize(n);
				for(std::size_t i=0; i<n; i++)
					destiny[i].insert(destiny[i].end(), source[i].begin(), source[i].end());
			}
			catch (std::bad_alloc)
			{
				assert(0);
				std::cerr << "bad_alloc caught"<< '\n';
				exit(1);
			}
		}

		template<typename valType, typename vectorType>
		inline void myCopy(vectorType &source, std::vector<valType> &destiny)
		{
			std::copy(source.begin(), source.end(), destiny.begin());
		}
#else
		template<typename valType>
		inline void myNew(valType* &v, std::size_t n)
		{
			if(n>0)
			{
				try{
					v=new valType[n];
				}
				catch (std::bad_alloc)
				{
					assert(0);
					std::cerr << "bad_alloc caught: size=" <<n<< '\n';
					exit(1);
				}
			}
		}

		template<typename valType>
		inline void myNew(valType* &v, std::size_t n, const valType &val)
		{
			if(n>0)
			{
				try{
					v=new valType[n];
					std::fill(v, v+n, val);
				}
				catch (std::bad_alloc)
				{
					assert(0);
					std::cerr << "bad_alloc caught: size=" <<n<< '\n';
					exit(1);
				}
			}
		}
		
		template<typename valType, typename vectorType>
		inline void myNew(valType* &destiny, vectorType &source)
		{
			try{
				v=new valType[source.size()];
				std::copy(source.begin(), source.end(), destiny);
			}
			catch (std::bad_alloc)
			{
				assert(0);
				std::cerr << "bad_alloc caught: size=" <<source.size()<< '\n';
				exit(1);
			}
		}

		template<typename valType>
		inline void myNew(valType** &v, std::size_t m, std::size_t n)
		{
			if(n>0)
			{
				try{
					v=new valType*[m];
					for(std::size_t i=0; i<m; i++)
						v[i]=new valType[n];
				}
				catch (std::bad_alloc)
				{
					assert(0);
					std::cerr << "bad_alloc caught: size=" <<m<<"*"<<n<<"="<<m*n<< '\n';
					exit(1);
				}
			}
		}
		
		template<typename valType>
		inline void myNew(valType** &v, std::size_t m, std::size_t n, const valType &val)
		{
			if(n>0)
			{
				try{
					v=new valType*[m];
					for(std::size_t i=0; i<m; i++)
					{
						v[i]=new valType[n];
						std::fill(v[i], v[i]+n, val);
					}
				}
				catch (std::bad_alloc)
				{
					assert(0);
					std::cerr << "bad_alloc caught: size=" <<m<<"*"<<n<<"="<<m*n<< '\n';
					exit(1);
				}
			}
		}
		
		template<typename valType, typename matrixType>
		inline void myNew(valType** &destiny, matrixType &source)
		{
			int m=source.size();
			try{
				destiny=new valType*[m];
				for(std::size_t i=0; i<m; i++)
				{
					destiny[i]=new valType[source[i].size()];
					std::copy(source[i].begin(), source[i].end(), destiny[i]);
				}
			}
			catch (std::bad_alloc)
			{
				assert(0);
				std::cerr << "bad_alloc caught"<< '\n';
				exit(1);
			}
		}

		template<typename valType, typename vectorType>
		inline void myCopy(vectorType &source, valType* &destiny)
		{
			std::copy(source.begin(), source.end(), destiny);
		}

		template<typename valType>
		inline void myDelete(valType** &v, std::size_t m)
		{
			int i;
			for(i=0; i<m; i++)
				delete[] v[i];
			delete[] v;
		}
#endif
#endif