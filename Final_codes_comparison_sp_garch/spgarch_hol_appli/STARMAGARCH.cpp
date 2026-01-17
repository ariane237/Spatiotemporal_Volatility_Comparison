#include <TMB.hpp>                                // Links in the TMB libraries
#include <math.h>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // DATA ARRAY/MATRIX
  //DATA_MATRIX(x);
  DATA_MATRIX(y);
  DATA_ARRAY(W);
  DATA_VECTOR(init);


  // DATA_MATRIX(w); //neighbour matrix size m x m
  //DATA_SCALAR(h0);

  // PARAMETERS
PARAMETER(mu)
  PARAMETER_MATRIX(phi);
  PARAMETER_MATRIX(theta);

  //PARAMETER(sigma);
  PARAMETER(omega);
  PARAMETER_MATRIX(alpha);
  PARAMETER_MATRIX(beta);

  //Type mu = 0.0;
  // dimension of x
    array<Type>x(y.rows(),y.cols());
  //array<Type>x2(y_dim(0),y_dim(1));
  array<Type>sigma(y.rows(),y.cols());
  array<Type>ell(y.rows(),y.cols()-1);

   // for(int k=0;k<y_dim(0);k++){
  //  x(k,0)=0.0;
  //  x2(k,0)=0.0;
  //  sigma(k,0)=init(k);
  //  }
  // Declare the "objective function" (neg. log. likelihood)
  Type f;
  f=0;
  //for (int t = 0; t < phi.cols();t++) {
  for(int k=0; k < y.rows(); k++) {
    x(k,0)=Type(0.0);//y(k,0)-mu;
    sigma(k,0)=init(k);
  }
  //}
  for (int t = 0; t < y.cols()-1; t++) {
    for(int k=0; k<y.rows(); k++){
      x(k,t+1) = y(k,t+1)-mu;
      for(int i = 0; i < phi.cols(); i++){
        for(int j = 0; j < phi.rows(); j++){
          for(int u =0; u < y.rows(); u++) {
            if(t-i>=0){
            x(k,t+1) -=  phi(j,i)*W(k,u,j) * (y(u,t-i)-mu);
            }
          }
        }
      }
      //x(k,t+1)= y(k,t+1)-yhat(k,t+1);
      for(int i = 0; i < theta.cols(); i++){
        for(int j = 0; j < theta.rows(); j++){
          for(int u =0; u < y.rows(); u++) {
            if(t-i>=0){
            x(k,t+1) -=  theta(j,i) * W(k,u,j) * x(u,t-i);
            }
          }
        }
      }
      }
  }
  for (int t = 0; t < y.cols()-1; t++) {
    for(int k=0; k<y.rows(); k++){
      sigma(k,t+1) = omega;
      for(int u =0; u < y.rows(); u++) {
      for(int i = 0; i < alpha.cols(); i++){
        for(int j = 0; j < alpha.rows(); j++){
          if(t-i>=0){
            sigma(k,t+1) +=  alpha(j,i) * W(k,u,j) * x(u,t-i) * x(u,t-i);
          }
          }
        }
      for(int i = 0; i < beta.cols(); i++){
        for(int j = 0; j < beta.rows(); j++){
          if(t-i>=0){
          sigma(k,t+1) +=  beta(j,i) * W(k,u,j) * sigma(u,t-i);
          }
        }
      }

      }
      f -= dnorm(x(k,t+1),Type(0.0),sqrt(sigma(k,t+1)),true);
      ell(k,t) = -dnorm(x(k,t+1),Type(0.0),sqrt(sigma(k,t+1)),true);
    }
  }


  // SIMULATE {
  //     for(int k =0; k < y.rows(); k++){
  //       x(k,0) = rnorm(Type(0.0),Type(1.0));
  //       y(k,0) = mu;
  //       sigma(k,0)=init(k);
  //     }
  // 
  //   for (int t = 0; t < y.cols()-1; t++) {
  //     for(int k=0; k < y.rows(); k++){
  //       sigma(k,t+1) = omega;
  //       for(int u =0; u < y.rows(); u++) {
  //         for(int i = 0; i < alpha.cols(); i++){
  //           for(int j = 0; j < alpha.rows(); j++){
  //             if(t-i >=0){
  //               sigma(k,t+1) +=  alpha(j,i) * W(k,u,j) * x(u,t-i)*x(u,t-i);
  //             }
  //            }
  //         }
  //         for(int i = 0; i < beta.cols(); i++){
  //           for(int j = 0; j < beta.rows(); j++){
  //             if(t-i >=0){
  //               sigma(k,t+1) +=  beta(j,i) * W(k,u,j) * sigma(u,t-i);
  //             }
  //           }
  //         }
  //       }
  //       x(k,t+1) = rnorm(Type(0.0), sqrt(sigma(k,t+1)));
  // 
  //     }
  //   }
  //   for (int t = 0; t < y.cols()-1; t++) {
  //     for(int k=0; k<y.rows(); k++){
  //       y(k,t+1) = mu + x(k,t+1);
  //       for(int u =0; u < y.rows(); u++) {
  //       for(int i = 0; i < phi.cols(); i++){
  //         for(int j = 0; j < phi.rows(); j++){
  //           if(t-i>=0){
  //               y(k,t+1) +=  phi(j,i)*W(k,u,j) * (y(u,t-i)-mu);
  //           }
  //         }
  //       }
  //       for(int i = 0; i < theta.cols(); i++){
  //         for(int j = 0; j < theta.rows(); j++){
  //           if(t-i>=0){
  //             y(k,t+1) +=  theta(j,i)*W(k,u,j) * x(u,t-i);
  //           }
  //         }
  //       }
  //       }
  //     }
  //   }
  // }
  REPORT(y);
  REPORT(x);
  REPORT(sigma);
  ADREPORT(ell);
  //REPORT(yhat);
  return f;
}


