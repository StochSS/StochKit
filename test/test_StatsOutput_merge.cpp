/*
 *  test_SSA_Direct.cpp
 *  
 */
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/thread.hpp>
#include <boost/shared_ptr.hpp>
#include <stdio.h>
#include <iomanip>
#include "SSA_Direct.h"
#include "TestModel.h"
#include "StatsAndTrajectoriesOutput.h"
#include "Random.h"

using namespace std;

extern bool veryClose(double,double,double);

int main()
{

  typedef boost::numeric::ublas::vector<double> denseVector;
  boost::shared_ptr<SimpleMessageHandler> messageHandler(new SimpleMessageHandler());
  
  typedef StatsAndTrajectoriesOutput<denseVector,SimpleMessageHandler> outputType;

  denseVector x(1);

  outputType x1(messageHandler);
  x1.initialize(4,0,1,x);
  outputType x2(messageHandler);
  x2.initialize(3,0,1,x);
  outputType x3(messageHandler);//combined obs of x1 and x2
  x3.initialize(7,0,1,x);

  x[0]=2.0;
  x1.record(0,0,x);
  x3.record(0,0,x);
  x[0]=4.0;
  x1.record(1,0,x);
  x3.record(1,0,x);
  x[0]=7.0;
  x1.record(2,0,x);
  x3.record(2,0,x);
  x[0]=16.0;
  x1.record(3,0,x);
  x3.record(3,0,x);
  x[0]=3.0;
  x2.record(0,0,x);
  x3.record(4,0,x);
  x[0]=9.0;
  x2.record(1,0,x);
  x3.record(5,0,x);
  x[0]=27.0;
  x2.record(2,0,x);
  x3.record(6,0,x);

  //verify
  std::size_t passedTests=0;
  if (x1.stats.data[1][0][0]==7.25) {
    passedTests++;
  }
  else {
    std::cout << "x1 mean FAILED" << std::endl;
  }
  if (x1.stats.data[0][0][0]==7.25) {
    passedTests++;
  }
  else {
    std::cout << "x1 'old mean' FAILED" << std::endl;
  }
  if (x2.stats.data[1][0][0]==13.0) {
    passedTests++;
  }
  else {
    std::cout << "x2 mean FAILED" << std::endl;
  }
  if (x2.stats.data[0][0][0]==13.0) {
    passedTests++;
  }
  else {
    std::cout << "x2 'old mean' FAILED" << std::endl;
  }

  if (veryClose(x1.stats.data[3][0][0],114.75,1E-13)) {
    passedTests++;
  }
  else {
    std::cout << "x1 'newS' FAILED" << std::endl;
  }

  if (veryClose(x1.stats.data[2][0][0],114.75,1E-13)) {
    passedTests++;
  }
  else {
    std::cout << "x1 'oldS' FAILED" << std::endl;
  }
  if (x2.stats.data[3][0][0]==312.0) {
    passedTests++;
  }
  else {
    std::cout << "x2 'newS' FAILED" << std::endl;
  }
  if (x2.stats.data[2][0][0]==312.0) {
    passedTests++;
  }
  else {
    std::cout << "x2 'oldS' FAILED" << std::endl;
  }

  x1.merge(x2);

  //if (x1.stats.data[1][0][0]==x3.stats.data[1][0][0]) {//will fail due to roundoff error
  if (veryClose(x1.stats.data[1][0][0],x3.stats.data[1][0][0],1E-13)) {
    passedTests++;
  }
  else {
    std::cout << "merge mean="<<std::setprecision(18)<<x1.stats.data[1][0][0]<<" should be "<<x3.stats.data[1][0][0]<<" FAILED" << std::endl;
  }
  if (veryClose(x1.stats.data[3][0][0]/7.0,x3.stats.data[3][0][0]/7.0,1E-13)) {
    passedTests++;
  }
  else {
    std::cout << "merge variance="<<std::setprecision(18)<<x1.stats.data[3][0][0]/7.0<<" should be "<<x3.stats.data[3][0][0]/7.0<<" FAILED" << std::endl;
  }

  std::cout << "passed " << passedTests << " out of 10 tests" << std::endl;
  if (passedTests==10) {
    std::cout << "PERFECT" << std::endl;
  }
  else {
    std::cout << "VALIDATION FAILED" << std::endl;
  }

  //do a controlled test
  //draw many values from round(Normal(10000,5))
  
  STOCHKIT::RandomGenerator gen;
  std::size_t totalObs=1000000;

  outputType out1(messageHandler);
  out1.initialize(totalObs/2,0,1,x);
  outputType out2(messageHandler);
  out2.initialize(totalObs/2,0,1,x);
  outputType master(messageHandler);//combined obs of x1 and x2
  master.initialize(totalObs,0,1,x);

  for (std::size_t i=0; i!=totalObs/2; ++i) {
    x[0]=gen.Normal(10000,10000);
    out1.record(i,0,x);
    master.record(i,0,x);
  }
  for (std::size_t i=0; i!=totalObs/2; ++i) {
    x[0]=gen.Normal(10000,10000);
    out2.record(i,0,x);
    master.record(i+totalObs/2,0,x);    
  }

  out1.merge(out2);
  std::cout << "compare means: true="<<std::setprecision(18)<<master.stats.data[1][0][0]<<" vs. merged="<<std::setprecision(18)<<out1.stats.data[1][0][0]<<std::endl;
  std::cout << "compare variances: true="<<master.stats.data[3][0][0]/(double)totalObs<<" vs. merged="<<out1.stats.data[3][0][0]/(double)totalObs<<std::endl;
  //merging two output objects seems to lead to accurate variance estimate--should look at more decimal places and try merging multiple output objects
  

  return 0;
}

bool veryClose(double x, double y,double eps=1E-13) {
  return fabs(x-y)<=eps;
}
