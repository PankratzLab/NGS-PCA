package org.pankratzlab.ngspca;

import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.QRDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;
import org.apache.commons.math3.random.MersenneTwister;
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
    RealMatrix Q = randn(A.getRowDimension(),
                         Math.min(A.getRowDimension(), numComponents + numOversamples));
    for (int i = 0; i < niters; i++) {
      Q = C.multiply(Q);
      Q = new LUDecomposition(Q).getL();
    }
    Q = C.multiply(Q);

    //    C.mmuli(Q, Q);

    Q = new QRDecomposition(Q).getQ();
    //    DoubleMatrix[] svd = Singular.fullSVD(Q.transpose().mmul(A));
    //    Returns:
    //      A DoubleMatrix[3] array of U, S, V such that A = U * diag(S) * V'
    SingularValueDecomposition svd = new SingularValueDecomposition(Q.transpose().multiply(A));
    RealMatrix W = Q.multiply(svd.getU());

    if (transpose) {
      for (int i = 0; i < numComponents; i++) {
        rsvd[0].putColumn(i, svd.getV().getColumn(i));
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

  private static RealMatrix randn(int rows, int columns) {
    RealMatrix m = MatrixUtils.createRealMatrix(rows, columns);
    MersenneTwister twister = new MersenneTwister(42);

    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < columns; j++) {
        m.setEntry(i, j, twister.nextDouble());
      }
      //      m.data[i] =twister.nextDouble();
    }

    return m;
  }
}
