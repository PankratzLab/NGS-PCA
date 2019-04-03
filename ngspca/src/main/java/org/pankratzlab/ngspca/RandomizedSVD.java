package org.pankratzlab.ngspca;

import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.QRDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;
import org.apache.commons.math3.random.MersenneTwister;

// https://raw.githubusercontent.com/yunjhongwu/matrix-routines/master/randomized_svd.java
// https://stackoverflow.com/questions/8677946/handle-large-sized-matrix-in-java/8678180 ... (Don't
// do it :) )
//
public class RandomizedSVD {

  //  https://arxiv.org/pdf/1608.02148.pdf
  //  https://arxiv.org/pdf/0909.4061.pdf
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
    rsvd[1] = MatrixUtils.createRealMatrix(numComponents, 1);
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
    Q = new QRDecomposition(Q).getQ();
    SingularValueDecomposition svd = new SingularValueDecomposition(Q.transpose().multiply(A));
    RealMatrix W = Q.multiply(svd.getU());

    if (transpose) {
      for (int i = 0; i < numComponents; i++) {
        rsvd[0].setColumn(i, svd.getV().getColumn(i));
        rsvd[1].setEntry(i, 0, svd.getSingularValues()[i]);
        rsvd[2].setColumn(i, W.getColumn(i));
      }
    } else {
      for (int i = 0; i < numComponents; i++) {
        rsvd[0].setColumn(i, W.getColumn(i));
        rsvd[1].setEntry(i, 0, svd.getSingularValues()[i]);
        rsvd[2].setColumn(i, svd.getV().getColumn(i));
      }
    }
  }

  /**
   * @param rows
   * @param columns
   * @return a {@link RealMatrix} populated with random values (deterministic random using
   *         {@link MersenneTwister})
   */
  private static RealMatrix randn(int rows, int columns) {
    RealMatrix m = MatrixUtils.createRealMatrix(rows, columns);
    MersenneTwister twister = new MersenneTwister(42);

    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < columns; j++) {
        m.setEntry(i, j, twister.nextDouble());
      }
    }
    return m;
  }
}
