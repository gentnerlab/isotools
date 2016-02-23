// $Id: pca2.cpp,v 1.1 2011/01/31 16:06:45 samn Exp $ 
#include <vector>
#include <fstream>
#include <algorithm>
#include "pca2.hpp"
#include "WCMath.h"
using namespace std;

void svd(int m, int n, double** u, double w[], double** v, int* ierr)
/*
 *   This subroutine is a translation of the Algol procedure svd,
 *   Num. Math. 14, 403-420(1970) by Golub and Reinsch.
 *   Handbook for Auto. Comp., Vol II-Linear Algebra, 134-151(1971).
 *
 *   This subroutine determines the singular value decomposition
 *        t
 *   A=usv  of a real m by n rectangular matrix, where m is greater
 *   then or equal to n.  Householder bidiagonalization and a variant
 *   of the QR algorithm are used.
 *  
 *
 *   On input.
 *
 *      m is the number of rows of A (and u).
 *
 *      n is the number of columns of A (and u) and the order of v.
 *
 *      u contains the rectangular input matrix A to be decomposed.
 *
 *   On output.
 *
 *      w contains the n (non-negative) singular values of a (the
 *        diagonal elements of s).  they are unordered.  if an
 *        error exit is made, the singular values should be correct
 *        for indices ierr+1,ierr+2,...,n.
 *
 *      u contains the matrix u (orthogonal column vectors) of the
 *        decomposition.
 *        if an error exit is made, the columns of u corresponding
 *        to indices of correct singular values should be correct.
 *
 *      v contains the matrix v (orthogonal) of the decomposition.
 *        if an error exit is made, the columns of v corresponding
 *        to indices of correct singular values should be correct.
 *
 *      ierr is set to
 *        zero       for normal return,
 *        k          if the k-th singular value has not been
 *                   determined after 30 iterations,
 *        -1         if memory allocation fails.
 *
 *   Questions and comments should be directed to B. S. Garbow,
 *   Applied Mathematics division, Argonne National Laboratory
 *
 *   Modified to eliminate machep
 *
 *   Translated to C by Michiel de Hoon, Human Genome Center,
 *   University of Tokyo, for inclusion in the C Clustering Library.
 *   This routine is less general than the original svd routine, as
 *   it focuses on the singular value decomposition as needed for
 *   clustering. In particular,
 *     - We require m >= n
 *     - We calculate both u and v in all cases
 *     - We pass the input array A via u; this array is subsequently
 *       overwritten.
 *     - We allocate for the array rv1, used as a working space,
 *       internally in this routine, instead of passing it as an
 *       argument. If the allocation fails, svd sets *ierr to -1
 *       and returns.
 *   2003.06.05
 */
{ int i, j, k, i1, k1, l1, its;
  double c,f,h,s,x,y,z;
  int l = 0;
  double g = 0.0;
  double scale = 0.0;
  double anorm = 0.0;
  double* rv1 = (double*)malloc(n*sizeof(double));
  if (!rv1)
  { *ierr = -1;
    return;
  }
  *ierr = 0;
  /* Householder reduction to bidiagonal form */
  for (i = 0; i < n; i++)
  { l = i + 1;
    rv1[i] = scale * g;
    g = 0.0;
    s = 0.0;
    scale = 0.0;
    for (k = i; k < m; k++) scale += fabs(u[k][i]);
    if (scale != 0.0)
    { for (k = i; k < m; k++)
      { u[k][i] /= scale;
        s += u[k][i]*u[k][i];
      }
      f = u[i][i];
      g = (f >= 0) ? -sqrt(s) : sqrt(s);
      h = f * g - s;
      u[i][i] = f - g;
      if (i < n-1)
      { for (j = l; j < n; j++)
        { s = 0.0;
          for (k = i; k < m; k++) s += u[k][i] * u[k][j];
          f = s / h;
          for (k = i; k < m; k++) u[k][j] += f * u[k][i];
        }
      }
      for (k = i; k < m; k++) u[k][i] *= scale;
    }
    w[i] = scale * g;
    g = 0.0;
    s = 0.0;
    scale = 0.0;
    if (i<n-1)
    { for (k = l; k < n; k++) scale += fabs(u[i][k]);
      if (scale != 0.0)
      { for (k = l; k < n; k++)
        { u[i][k] /= scale;
          s += u[i][k] * u[i][k];
        }
        f = u[i][l];
        g = (f >= 0) ? -sqrt(s) : sqrt(s);
        h = f * g - s;
        u[i][l] = f - g;
        for (k = l; k < n; k++) rv1[k] = u[i][k] / h;
        for (j = l; j < m; j++)
        { s = 0.0;
          for (k = l; k < n; k++) s += u[j][k] * u[i][k];
          for (k = l; k < n; k++) u[j][k] += s * rv1[k];
        }
        for (k = l; k < n; k++)  u[i][k] *= scale;
      }
    }
    anorm = max(anorm,fabs(w[i])+fabs(rv1[i]));
  }
  /* accumulation of right-hand transformations */
  for (i = n-1; i>=0; i--)
  { if (i < n-1)
    { if (g != 0.0)
      { for (j = l; j < n; j++) v[j][i] = (u[i][j] / u[i][l]) / g;
        /* double division avoids possible underflow */
        for (j = l; j < n; j++)
        { s = 0.0;
          for (k = l; k < n; k++) s += u[i][k] * v[k][j];
          for (k = l; k < n; k++) v[k][j] += s * v[k][i];
        }
      }
    }
    for (j = l; j < n; j++)
    { v[i][j] = 0.0;
      v[j][i] = 0.0;
    }
    v[i][i] = 1.0;
    g = rv1[i];
    l = i;
  }
  /* accumulation of left-hand transformations */
  for (i = n-1; i >= 0; i--)
  { l = i + 1;
    g = w[i];
    if (i!=n-1)
      for (j = l; j < n; j++) u[i][j] = 0.0;
    if (g!=0.0)
    { if (i!=n-1)
      { for (j = l; j < n; j++)
        { s = 0.0;
          for (k = l; k < m; k++) s += u[k][i] * u[k][j];
          /* double division avoids possible underflow */
          f = (s / u[i][i]) / g;
          for (k = i; k < m; k++) u[k][j] += f * u[k][i];
        }
      }
      for (j = i; j < m; j++) u[j][i] /= g;
    }
    else
      for (j = i; j < m; j++) u[j][i] = 0.0;
    u[i][i] += 1.0;
  }
  /* diagonalization of the bidiagonal form */
  for (k = n-1; k >= 0; k--)
  { k1 = k-1;
    its = 0;
    while(1)
    /* test for splitting */
    { for (l = k; l >= 0; l--)
      { l1 = l-1;
        if (fabs(rv1[l]) + anorm == anorm) break;
        /* rv1[0] is always zero, so there is no exit
         * through the bottom of the loop */
        if (fabs(w[l1]) + anorm == anorm)
        /* cancellation of rv1[l] if l greater than 0 */
        { c = 0.0;
          s = 1.0;
          for (i = l; i <= k; i++)
          { f = s * rv1[i];
            rv1[i] *= c;
            if (fabs(f) + anorm == anorm) break;
            g = w[i];
            h = sqrt(f*f+g*g);
            w[i] = h;
            c = g / h;
            s = -f / h;
            for (j = 0; j < m; j++)
            { y = u[j][l1];
              z = u[j][i];
              u[j][l1] = y * c + z * s;
              u[j][i] = -y * s + z * c;
            }
          }
          break;
        }
      }
      /* test for convergence */
      z = w[k];
      if (l==k) /* convergence */
      { if (z < 0.0)
        /*  w[k] is made non-negative */
        { w[k] = -z;
          for (j = 0; j < n; j++) v[j][k] = -v[j][k];
        }
        break;
      }
      else if (its==30)
      { *ierr = k;
        break;
      }
      else
      /* shift from bottom 2 by 2 minor */
      { its++;
        x = w[l];
        y = w[k1];
        g = rv1[k1];
        h = rv1[k];
        f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
        g = sqrt(f*f+1.0);
        f = ((x - z) * (x + z) + h * (y / (f + (f >= 0 ? g : -g)) - h)) / x;
        /* next qr transformation */
        c = 1.0;
        s = 1.0;
        for (i1 = l; i1 <= k1; i1++)
        { i = i1 + 1;
          g = rv1[i];
          y = w[i];
          h = s * g;
          g = c * g;
          z = sqrt(f*f+h*h);
          rv1[i1] = z;
          c = f / z;
          s = h / z;
          f = x * c + g * s;
          g = -x * s + g * c;
          h = y * s;
          y = y * c;
          for (j = 0; j < n; j++)
          { x = v[j][i1];
            z = v[j][i];
            v[j][i1] = x * c + z * s;
            v[j][i] = -x * s + z * c;
          }
          z = sqrt(f*f+h*h);
          w[i1] = z;
          /* rotation can be arbitrary if z is zero */
          if (z!=0.0)
          { c = f / z;
            s = h / z;
          }
          f = c * g + s * y;
          x = -s * g + c * y;
          for (j = 0; j < m; j++)
          { y = u[j][i1];
            z = u[j][i];
            u[j][i1] = y * c + z * s;
            u[j][i] = -y * s + z * c;
          }
        }
        rv1[l] = 0.0;
        rv1[k] = f;
        w[k] = x;
      }
    }
  }
  free(rv1);
  return;
}

PCA::PCA(int dim) {
  init(dim);
}

void PCA::init(int dim)
{
  dim_=dim;
  mean_=vector<double>(dim_,0.0);
  //covariance_.Init(dim_,dim_);
  covariance_=A2D<double>(dim_,dim_);
  covariance_.Fill(0.0);
  counter_=0;
}

void PCA::save(const string &filename) const{
  /*ogzstream ofs; ofs.open(filename.c_str());
  if(!ofs.good()) {
    ERR << "[PCA::save] Cannot write PCA to '" << filename << "'."<< endl;
  } else {
    ofs << dim_ << " " <<dim_ << endl;
    ofs << "mean "  ;
    for(int i=0;i<dim_;++i) {
      ofs << mean_[i] << " "; }
    ofs << endl;

    for(int y=0;y<dim_;++y) {
      ofs << y;
      for(int x=0;x<dim_;++x) {
        ofs << " " << covariance_[y][x];
      }
      ofs << endl;
    }
    ofs << -1 << endl;
  }*/
}

void PCA::load(const string &filename) {
  /*igzstream ifs; ifs.open(filename.c_str());
  int dimy, dimx, posy;

  if(!ifs.good()) {
    ERR << "[PCA::load] Cannot open '"<<filename<<"' to read PCA." << endl;
  } else {
    ifs >> dimx >> dimy;
    if(dimx==dimy) dim_=dimx;
    else ERR << "[PCA::load] Strange: PCA not square: dimx="<< dimx << " dimy=" << dimy << endl;
    ::std::string tmpstr;
    ifs >> tmpstr;
    mean_=vector<PCA_T>(dimy);
    for(int i=0;i<dimy;++i) {
      ifs >> mean_[i];
    }
    covariance_=vector< vector<PCA_T> >(dimy,vector<PCA_T>(dimx,0.0));
    for(int y=0;y<dimy;++y) {
      ifs >> posy;
      if(posy!=y) {ERR << "[PCA::load] Strangeness in reading PCA." << endl;}
      for(int x=0;x<dimx;++x) {
        ifs >> covariance_[y][x];
      }
    }
    ifs >> posy;
    if(posy!=-1) {ERR << "[PCA::load] Reading PCA was strange. Did not end with -1." << endl;}
  }*/
}

void PCA::dataEnd() {
  int i,j;
  
  for(i=0;i<dim_;++i) mean_[i]/=counter_;

  for(i=0;i<dim_;++i)
    for(j=0;j<dim_;++j) 
      covariance_[i][j]/=counter_;

  for(i=0;i<dim_;++i) 
    for(j=0;j<dim_;++j) 
      covariance_[i][j]-=mean_[i]*mean_[j];
}

void PCA::calcPCA() {

  vector<double> w(dim_);
  A2D<double> v(dim_,dim_);
  v.Fill(0.0);
/*
  double** cov = covariance_.GetP();
  Write2LogPlain("covariance:\n");
  for(int i=0;i<dim_;i++){ for(int j=0;j<dim_;j++){
	  Write2LogPlain("%.4f ",cov[i][j]);
  } Write2LogPlain("\n");}
*/
  int ierr = 0;
  svd(w.size(),w.size(),covariance_.GetP(),&w[0],v.GetP(),&ierr);

  //negate PCs so will be consistent with Matlab
  for(int i=0;i<dim_;i++) for(int j=0;j<dim_;j++){
	  covariance_[i][j]=-covariance_[i][j];
  }

  //NR::svdcmp(covariance_, w, v);
/*
  Write2LogPlain("Eigenvalues: ");
  for(unsigned int i=0;i<w.size();++i) {
    Write2LogPlain("%g ",w[i]);
  } Write2LogPlain("\n");
  Write2LogPlain("sum of eigenvalues = %g",Sum(w));

  Write2LogPlain("transform before sorting:\n");
  for(int i=0;i<dim_;++i) {
    for(int j=0;j<dim_;++j) {
	  Write2LogPlain(" %g",covariance_[i][j]);
  }Write2LogPlain("\n");
  }
*/
  vector< pair<double, vector<double> > > toSort;
  for(unsigned int i=0;i<w.size();++i) {
    vector<double> tmp(dim_);
    for(int j=0;j<dim_;++j) tmp[j]=covariance_[j][i];
    toSort.push_back(pair<double,vector<double> >(w[i],tmp));
  }

  sort(toSort.rbegin(),toSort.rend());

  for(unsigned int i=0;i<toSort.size();++i) {
    copy(toSort[i].second.begin(),toSort[i].second.end(),covariance_[i]);
    w[i]=toSort[i].first;
  }
/*
  Write2LogPlain("Eigenvalues: ");
  for(unsigned int i=0;i<w.size();++i) {
    Write2LogPlain("%g ",w[i]);
  } Write2LogPlain("\n");
  Write2LogPlain("sum of eigenvalues = %g",Sum(w));

  Write2LogPlain("transform after sorting:\n");
  for(int i=0;i<dim_;++i) {
    for(int j=0;j<dim_;++j) {
	  Write2LogPlain(" %g",covariance_[i][j]);
  }Write2LogPlain("\n");
  }
*/
  const bool bSave = false;
  if(bSave)
  {	FILE* fp = fopen("__pcs.txt","w");
	for(int r=0;r<covariance_.Rows();r++){
	  for(int c=0;c<covariance_.Cols();c++)
		  fprintf(fp,"%g ",covariance_[r][c]);
	  fprintf(fp,"\n");
	}
	fclose(fp);
  }

  eigenvalues_=w;

}

const int PCA::dim() const {
  return dim_;
}

const vector<double> & PCA::mean() const {
  return mean_;
}

const A2D<double>& PCA::covariance() const {
  return covariance_;
}

const int PCA::counter() const {
  return counter_;
}
