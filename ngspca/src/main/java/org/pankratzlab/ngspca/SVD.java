package org.pankratzlab.ngspca;

import java.io.Serializable;
import java.util.logging.Logger;
import org.ejml.alg.dense.decomposition.svd.SvdImplicitQrDecompose_D64;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.SingularOps;

/**
 * class to perform singular value decomposition on {@link DenseMatrix64F}
 */
class SVD implements Serializable {

  /**
   * 
   */
  private static final long serialVersionUID = 5517545544997813383L;
  private final String[] colNames;
  private final String[] rowNames;

  private SvdImplicitQrDecompose_D64 svd;

  /**
   * @param colNames
   * @param rowNames
   */
  SVD(String[] colNames, String[] rowNames) {
    super();
    this.colNames = colNames;
    this.rowNames = rowNames;
  }

  /**
   * @param dm {@link DenseMatrix64F} that has been pre-normalized, if desired
   * @param log
   */
  void computeSVD(DenseMatrix64F dm, Logger log) {
    if (dm.getNumCols() != colNames.length) {
      throw new IllegalArgumentException("Mismatched column lengths");
    }
    if (dm.getNumRows() != rowNames.length) {
      throw new IllegalArgumentException("Mismatched row lengths");
    }
    log.info("Computing EJML PCs");
    this.svd = new SvdImplicitQrDecompose_D64(false, true, true, false);
    svd.decompose(dm);
    log.info("Finished Computing EJML PCs");

    log.info("Sorting matrices to descending order");

    SingularOps.descendingOrder(svd.getU(null, false), false, svd.getW(null), svd.getV(null, true),
                                true);

  }

  /**
   * @param numComponents number of singular values to return from the W matrix
   * @return
   */
  double[] getSingularValues(int numComponents) {
    double[] singularValues = new double[numComponents];
    for (int i = 0; i < numComponents; i++) {
      singularValues[i] = svd.getW(null).get(i, i);
    }
    return singularValues;
  }

  private void dumpPCsToText(String file) {
    //
    // PrintWriter writer = Files.getAppropriateWriter(file);
    // StringJoiner joiner = new StringJoiner("\t");
    // joiner.add("SAMPLE");
    // for (int i = 0; i < v.getColumnDimension(); i++) {
    // joiner.add("PC" + (i + 1));
    // }
    // writer.println(joiner.toString());
    //
    // for (int i = 0; i < v.getRowDimension(); i++) {
    // StringJoiner sample = new StringJoiner("\t");
    // sample.add(colNames[i]);
    // for (int j = 0; j < v.getColumnDimension(); j++) {
    // sample.add(Double.toString(v.getEntry(j, i)));
    // }
    // writer.println(sample.toString());
    //
    // }
    // writer.close();
  }
}
