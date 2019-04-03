package org.pankratzlab.ngspca;

import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
// https://raw.githubusercontent.com/yunjhongwu/matrix-routines/master/randomized_svd.java
// https://stackoverflow.com/questions/8677946/handle-large-sized-matrix-in-java/8678180
// org.jbl
import org.jblas.Decompose;
import org.jblas.DoubleMatrix;
import org.jblas.Singular;

public class RandomizedSVD {

  // Compute a (truncated) randomized SVD of a JBLAS DoubleMatrix
  private int numComponents = 2;
  private int niters = 5;

  private int numOversamples = 10;
  private boolean transpose = false;
  public RealMatrix[] rsvd = new RealMatrix[3];

  public RandomizedSVD(int numComponents, int niters) {
    this.numComponents = numComponents;
    this.niters = niters;
  }

  public void fit(RealMatrix A) {
    transpose = A.getRowDimension() > A.getColumnDimension();
    rsvd[0] = MatrixUtils.createRealMatrix(A.getRowDimension(), numComponents);
    rsvd[1] = MatrixUtils.createRealMatrix(1, numComponents);
    rsvd[2] = MatrixUtils.createRealMatrix(A.getColumnDimension(), numComponents);
    if (transpose) {
      A = A.transpose();
    }

    RealMatrix C = A.multiply(A.transpose());
    RealMatrix Q = DoubleMatrix.randn(A.rows, Math.min(A.rows, numComponents + numOversamples));
    for (int i = 0; i < niters; i++) {
      Q=C.multiply(Q);
    
      Q =  new LUDecomposition(Q).getL();
    }
    Q=C.multiply(Q);

//    C.mmuli(Q, Q);
    Q = Decompose.qr(Q).q;
    DoubleMatrix[] svd = Singular.fullSVD(Q.transpose().mmul(A));
    DoubleMatrix W = Q.mmul(svd[0]);

    if (transpose) {
      for (int i = 0; i < numComponents; i++) {
        rsvd[0].putColumn(i, svd[2].getColumn(i));
        rsvd[1].put(i, svd[1].get(i));
        rsvd[2].putColumn(i, W.getColumn(i));
      }
    } else {
      for (int i = 0; i < numComponents; i++) {
        rsvd[0].putColumn(i, W.getColumn(i));
        rsvd[1].put(i, svd[1].get(i));
        rsvd[2].putColumn(i, svd[2].getColumn(i));
      }
    }
  }
}
