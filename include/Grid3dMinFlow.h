#ifndef GRID3DMINFLOW_H
#define GRID3DMINFLOW_H

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>

#include <map>
#include <vector>
#include <iostream>

#include "LPSolver.h"


template <typename TPrecision>
class Grid3dMinFlow {




  public:
    typedef typename Eigen::Tensor<TPrecision, 3> TensorXp;
    typedef typename Eigen::Tensor<int, 3> TensorXi;
    typedef typename Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> MatrixXi;


  private:
    TPrecision cost;
    MatrixXi id2spatial;

  public:

    Grid3dMinFlow(){
      cost=-1;
    };

    virtual ~Grid3dMinFlow(){
    };





    std::map< std::pair<long, long>, double > solve(const TensorXp &X1,
        const TensorXp &X2, LPSolver &solver, double lambda=0, bool scaleMass = true){


      auto dims = X1.dimensions();
      long nNodes = dims[0]*dims[1]*dims[2];
      TensorXi nId(dims[0], dims[1], dims[2]);


      std::vector<TPrecision> sources( nNodes );
      std::vector<TPrecision> targets( nNodes );

      //These statements crash
      //TensorXp X1n = X1 / X1.sum();
      //TensorXp X2n = X2 / X2.sum();

      Eigen::Tensor<TPrecision, 0> s1 = X1.sum();
      Eigen::Tensor<TPrecision, 0> s2 = X2.sum();
      TPrecision sum1 = s1(0);
      TPrecision sum2 = s2(0);

      long nSources = 0;
      long nTargets = 0;

      MatrixXi sid2spatial = MatrixXi( nNodes, 3 );
      MatrixXi tid2spatial = MatrixXi( nNodes, 3 );
      for(int i=0; i<dims[0]; i++){
        for(int j=0; j<dims[1]; j++){
          for(int k=0; k<dims[2]; k++){
            TPrecision v1 = X1(i, j, k);
            TPrecision v2 = X2(i, j, k);
            nId(i, j, k) = -1;

            if(v1 == 0 && v2 == 0){
              continue;
            }


            //diff[idCounter] = X1n(i, j, k) - X2n(i, j, k);
            TPrecision diff = 0;
            if(scaleMass){
              diff = v1 / sum1 - v2 / sum2;
            }
            else{
              diff = v1 - v2;
            }
            if( diff < 0 ){
              targets[ nTargets ] = diff;
              tid2spatial( nTargets, 0 ) = i;
              tid2spatial( nTargets, 1 ) = j;
              tid2spatial( nTargets, 2 ) = k;
              nTargets++;
            }
            else{
              sources[ nSources ] = diff;
              sid2spatial( nSources, 0 ) = i;
              sid2spatial( nSources, 1 ) = j;
              sid2spatial( nSources, 2 ) = k;
              nSources++;
            }
          }
        }
      }

      sources.resize( nSources );
      targets.resize( nTargets );


      id2spatial = MatrixXi( nSources + nTargets, 3 );
      for(long n=0; n < nSources; n++){
        int i = sid2spatial( n, 0 );
        int j = sid2spatial( n, 1 );
        int k = sid2spatial( n, 2 );
        id2spatial( n, 0 ) = i;
        id2spatial( n, 1 ) = j;
        id2spatial( n, 2 ) = k;
        nId(i, j, k) = n;
      }
      for(long n=0; n < nTargets; n++){
        int i = tid2spatial( n, 0 );
        int j = tid2spatial( n, 1 );
        int k = tid2spatial( n, 2 );
        id2spatial( n + nSources, 0 ) = i;
        id2spatial( n + nSources, 1 ) = j;
        id2spatial( n + nSources, 2 ) = k;
        nId( i, j, k ) = n + nSources;
      }

      std::cout << "nSources: " << nSources << std::endl;
      std::cout << "nSources: " << sources.size() << std::endl;
      std::cout << "nTargets: " << nTargets << std::endl;
      std::cout << "nTargets: " << targets.size() << std::endl;


      std::vector< long > asInd( nNodes * 27 );
      std::vector< long > atInd( nNodes * 27 );
      std::vector< double > coeff( nNodes * 27 );

      long index = 0;
      for(int i=0; i<dims[0]; i++){
        for(int j=0; j<dims[1]; j++){
          for(int k=0; k<dims[2]; k++){

            for(int ii=-1; ii < 2; ii+=1){
              int iii = i + ii;
              if(iii < 0 || iii == dims[0]){
                continue;
              }

              for(int jj=-1; jj < 2; jj+=1){
                int jjj = j + jj;
                if(jjj < 0 || jjj == dims[1] ){
                  continue;
                }

                for(int kk=-1; kk < 2; kk+=1){
                  int kkk = k + kk;
                  if(kkk < 0 || kkk == dims[2]){
                    continue;
                  }


                  if( iii == i && jjj == j && kkk == k){
                    continue;
                  }
                  long id1 = nId(i, j, k);
                  long id2 = nId(iii, jjj, kkk);
                  if( id1 == -1 || id2 == -1){
                    continue;
                  };
                  asInd[index] = id1;
                  atInd[index] = id2;

                  coeff[index] = sqrt( ii*ii + jj*jj + kk*kk);
                  if(lambda > 0){
                    coeff[index] += lambda * fabs( X1(i, j, k) - X2(iii, jjj, kkk) );
                  }

                  ++index;

                }
              }
            }

          }
        }
      }


      std::cout << index << std::endl;

      asInd.resize(index);
      atInd.resize(index);
      coeff.resize(index);

      std::cout << coeff.size() << std::endl;

      solver.createLP( nSources, nTargets );
      solver.addRows( nSources+nTargets );
      for(long i=0; i<nSources; i++){
         solver.setRowBounds(i, sources[i]);
      }
      for(long i=0; i<nTargets; i++){
         solver.setRowBounds( i +nSources, targets[i]);
      }

      solver.addColumns( coeff.size() );
      for(long i=0; i<coeff.size(); i++){
        solver.setColumnBoundsLower(i, 0);
        solver.setColumnObjective(i, coeff[i] );
        solver.setColumnCoefficients( i, asInd[i], atInd[i]);
      }

      solver.solveLP();

      std::cout << "Solution status: " << solver.isOptimal() << std::endl;

      std::map<std::pair<long, long>, double> plan;
      index = 0;
      for(long i = 0; i < asInd.size(); i++){
        double w = solver.getColumnPrimal(i);
        if(w > 0){
          plan[ std::pair<long, long>(asInd[i], atInd[i]) ] = w;
        }
      }

      return plan;
    };


    TPrecision getCost(){
      return cost;
    };

    MatrixXi &getId2Spatial(){
      return id2spatial;
    };



};


#endif


