// $Id: pca2.hpp,v 1.1 2011/01/31 16:06:58 samn Exp $ 
#ifndef __pca_hpp__
#define __pca_hpp__

#include "diag.hpp"
#include "A2D.h"
#include "WCMath.h"
#include <vector>

typedef float PCA_T;

/** implementation of pca transformation. to use this class do it in
    this order:

    1. create a PCA object with PCA(dim) where dim is the
       dimensionality of the data to be transformed.

    2. then you can use putData to provide data for estimation of the
       covariance matrix.

    3. once you have put all data using the putData-method you have to
       tell the pca that you are ready: dataEnd(). This terminates the
       calculation of the covariance matrix.

    4. now the next step is to call calcPCA which uses methods from
       lapack to perform singular value decomposition

    5. when this is finished you can use the methods transform and
       backtransform to transform vectors. using transform you can create
       vectors of arbitrary dimension (<dim) which best represent the
       given input under the given assumptions (e.g. train data from
       putData). And you can use backtransform to transform arbitrary
       small vectors back into the original vector space.

    NOTE: you can load and save your PCA models. but this makes sense
          only after calcPCA.

*/


class PCA {
public:

  /** constructor which also makes the object know the dimensionality which is to be used. */
  PCA(int dim=0);

  void init(int dim);

  /** load the transformation from a file.
   */
  void load(const ::std::string &filename);

  /** save the transformation to a file.
   */
  void save(const ::std::string &filename) const;

  template< class T >
  void putData(const ::std::vector<T> &in)
  {	++counter_;
	T tmp;	
	int i,j;
	
	for(i=0;i<dim_;++i)	mean_[i]+=in[i];

	for(i=0;i<dim_;++i) {
		tmp=in[i];
		for(j=0;j<dim_;++j)
			covariance_[i][j]+=tmp*in[j];
	}
  }

  template< class T >
  void putData(const T* in)
  {
	++counter_;
	T tmp;
	int i,j;

	for(i=0;i<dim_;++i) mean_[i]+=in[i];

	for(i=0;i<dim_;++i) {
		tmp=in[i];
		for(j=0;j<dim_;++j) 
			covariance_[i][j]+=tmp*in[j];
	}
  }

  template< class T > 
  void putData(const ::std::vector< std::vector< T > >& vData)
  {
	counter_ = vData.size();
	mean_ = vector<double>(dim_,0.0);
	CovarMat<T>(vData,counter_,dim_,covariance_,mean_);
  }

  template< class T >
  void putData( A2D< T >& vData)
  {
	counter_ = vData.Rows();
	mean_ = vector<double>(dim_,0.0);
	CovarMat(vData,counter_,dim_,covariance_,mean_);
  }

  void putData( A2D< float >& vData)
  {
	counter_ = vData.Rows();
	mean_ = vector<double>(dim_,0.0);
	CovarMat<float,double>(vData,counter_,dim_,covariance_,mean_);
  }

  /** make the object know, that there will be no further data. After
      this it is not a good idea to call putData again, and also it is
      probably not a good idea to call calcPCA, mean, covariance
      before this. */
  void dataEnd();

  /** start the svd decomposition.
   *
   * Attention: Does not take care whether the data is in a valid
   * state for this. Will apply the svd to anything contained in the
   * covariance. Before this is started the data has to be put into
   * with putData, and the process of doing this has to be finished
   * with dataEnd();
   */
  void calcPCA();

  /** return the dimensionality of the data to be put into */
  const int dim() const;

  /** return the mean
   * Attention: does not take care whether the process of estimation was
   * already ended. Will return whatever is contained in the mean matrix.
   */
  const ::std::vector<double> & mean() const;

  /** return the covariance matrix
   *
   * Attention: does not take care whether the process of estimation was
   * already ended. Will return whatever is contained in the covariance matrix.
   */
  const A2D<double>& covariance() const ;

  /** return the eigenvalues */
  const ::std::vector<double>& eigenvalues() const {return eigenvalues_;}

  /** return the eigenvector i*/
  const double* eigenvector(const int i) const {return covariance_[i];}

  /** how many data elements have been put into this pca. */
  const int counter() const;

  /**
   * transform a vector due to the calculated pca.
   * @param in vector to be transformed
   * @param dim the dimensionality of the vector to be returned. If 0 no dimensionality reduction will be applied.
   *
   * Attention: This method does not take care if the pca is in a
   * valid state for transformation. If the svd has not yet been
   * applied the vector will be multiplied by whatever is contained in
   * the covariance matrix.
   */
  template< class T >
  const ::std::vector<T> transform(const ::std::vector<T> &in, int dim=0) const
  {
	if(dim==0) dim=in.size();
	vector<PCA_T> result(dim,0.0);
	DBG(20) << "Transforming to dim " << dim << endl;

	unsigned int i,j;

	for(i=0;i<dim;++i) {
		for(j=0;j<in.size();++j) 
			result[i]+=covariance_[i][j]*(in[j]-mean_[j]);
	}
	return result;
  }

  template< class T >
  void transform(T* in, int insz, T* out, int dim=0) const
  {
	if(dim==0) dim=insz;
	DBG(20) << "Transforming to dim " << dim << endl;

	unsigned int i,j;

	for(i=0;i<dim;++i) {
		out[i]=0.0;
		for(j=0;j<insz;++j)
			out[i]+=covariance_[i][j]*(in[j]-mean_[j]);
	}
  }


  /**
   * transform a vector from pca-transformed space back into the
   * untransformed space. If the input vector is shorter than the
   * dimensionality of the untransformed space, zero padding is
   * automatically done.
   *
   * @param in vector to be transformed
   *
   * Attention: This method does not take care if the pca is in a
   * valid state for transformation. If the svd has not yet been
   * applied the vector will be multiplied by whatever is contained in
   * the covariance matrix.
   */
	/** as PCA-matrices are rotation matrices, for this calculation it is
	*  not necessary to explicitly invert them, but instead it holds:
	*
	*   M^{-1}=transpose(M), and thus M×transpose(M)=Id
	*/
  template< class T >
  const ::std::vector<T> backTransform(::std::vector<T> in) const
  {
	if(int(in.size()) < dim_) { in.resize(dim_,0);} // zeropadding
	vector<PCA_T> result(dim_,0.0);
	DBG(20) << "Backtransforming" << endl;

	int i,j;

	for(i=0;i<dim_;++i) {
		for(j=0;j<dim_;++j)
			result[i]+=covariance_[j][i]*in[j];
		result[i]+=mean_[i];
	}
	return result;
  }


private:
  int dim_;
  long int counter_;

  ::std::vector<double> mean_;
  A2D<double> covariance_;
  ::std::vector<double> eigenvalues_;
};

#endif
