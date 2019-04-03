package org.pankratzlab.ngspca;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.charset.Charset;
import java.util.StringJoiner;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.io.FileUtils;
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

  private static final int RANDOM_SEED = 42;
  //  https://arxiv.org/pdf/1608.02148.pdf
  //  https://arxiv.org/pdf/0909.4061.pdf
  // Compute a (truncated) randomized SVD of a JBLAS DoubleMatrix
  private int numComponents = 2;
  private int niters = 5;

  private int numOversamples = 10;
  private boolean transpose = false;
  private RealMatrix[] rsvd = new RealMatrix[3];
  private final Logger log;

  /**
   * Column names of the original input data
   */
  private final String[] originalColNames;
  /**
   * Row names of the original input data
   */
  private final String[] originalRowNames;

  public RandomizedSVD(String[] originalColNames, String[] originalRowNames, int numComponents,
                       int niters, Logger log) {
    this.numComponents = numComponents;
    this.niters = niters;
    this.originalColNames = originalColNames;
    this.originalRowNames = originalRowNames;
    this.log = log;
  }

  public void fit(RealMatrix A) {
    log.info("Initializing matrices");
    transpose = A.getRowDimension() > A.getColumnDimension();
    rsvd[0] = MatrixUtils.createRealMatrix(A.getRowDimension(), numComponents);
    rsvd[1] = MatrixUtils.createRealMatrix(numComponents, 1);
    rsvd[2] = MatrixUtils.createRealMatrix(A.getColumnDimension(), numComponents);
    if (transpose) {
      log.info("Transposing, since row N > column N");
      A = A.transpose();
    }

    RealMatrix C = A.multiply(A.transpose());
    log.info("Selecting randomized Q");

    RealMatrix Q = randn(A.getRowDimension(),
                         Math.min(A.getRowDimension(), numComponents + numOversamples));

    log.info("Beginning LU decomp iterations");
    for (int i = 0; i < niters; i++) {
      Q = C.multiply(Q);
      Q = new LUDecomposition(Q).getL();
    }
    Q = C.multiply(Q);
    Q = new QRDecomposition(Q).getQ();
    log.info("SVD of reduced matrix");

    SingularValueDecomposition svd = new SingularValueDecomposition(Q.transpose().multiply(A));

    RealMatrix W = Q.multiply(svd.getU());

    log.info("Setting SVD V/W/U results");

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
      log.info("Finished SVD");
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
    MersenneTwister twister = new MersenneTwister(RANDOM_SEED);

    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < columns; j++) {
        m.setEntry(i, j, twister.nextDouble());
      }
    }
    return m;
  }

  /**
   * @param file dump the PCs to this text file
   * @param log
   */
  void dumpPCsToText(String file, Logger log) {
    //
    RealMatrix v = transpose ? rsvd[0] : rsvd[2];

    String[] pcNames = SVD.getNumberedColumnHeader("PC", v.getRowDimension());

    dumpMatrix(file, v, "SAMPLE", pcNames, originalColNames, true, log);
  }

  private static void dumpMatrix(String file, RealMatrix m, String rowTitle,
                                 String[] outputColumnNames, String[] outputRowNames,
                                 boolean transposed, Logger log) {
    try (PrintWriter writer = new PrintWriter(new File(file))) {

      StringJoiner joiner = new StringJoiner("\t");
      joiner.add(rowTitle);
      for (String colName : outputColumnNames) {
        joiner.add(colName);
      }
      writer.println(joiner.toString());
      log.info(outputRowNames.length + " output rows by " + outputColumnNames.length
               + " output columns");

      for (int outputRow = 0; outputRow < outputRowNames.length; outputRow++) {
        StringJoiner sample = new StringJoiner("\t");
        sample.add(outputRowNames[outputRow]);
        for (int outputColumn = 0; outputColumn < outputColumnNames.length; outputColumn++) {
          if (transposed) {
            sample.add(Double.toString(m.getEntry(outputColumn, outputRow)));
          } else {
            sample.add(Double.toString(m.getEntry(outputRow, outputColumn)));
          }
        }
        writer.println(sample.toString());

      }
    } catch (FileNotFoundException e) {

      log.log(Level.SEVERE, "unable to write to file " + file, e);

    }
  }

  /**
   * @param file loadings will be computed and dumped to this file
   * @param log
   */
  void computeAndDumpLoadings(String file, Logger log) {
    RealMatrix loadingData = transpose ? rsvd[2] : rsvd[0];
    String[] loadingNames = SVD.getNumberedColumnHeader("Loading",
                                                        loadingData.getColumnDimension());
    dumpMatrix(file, loadingData, "MARKER", loadingNames, originalRowNames, false, log);

  }

  /**
   * @param file singular values will be dumped to this file
   * @param log
   */
  void dumpSingularValuesToText(String file, Logger log) {
    StringJoiner joiner = new StringJoiner("\n");
    joiner.add("SINGULAR_VALUES\tPC");

    for (int component = 0; component < numComponents; component++) {
      joiner.add(component + "\t" + Double.toString(rsvd[1].getEntry(component, 0)));
    }

    try {
      FileUtils.writeStringToFile(new File(file), joiner.toString(), Charset.defaultCharset(),
                                  false);
    } catch (IOException e) {
      log.log(Level.SEVERE, "unable to write to file " + file, e);
    }

  }
}
