package org.pankratzlab.ngspca;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.Serializable;
import java.util.HashMap;
import java.util.Map;
import java.util.StringJoiner;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.ejml.alg.dense.decomposition.svd.SvdImplicitQrDecompose_D64;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;
import org.ejml.ops.SingularOps;

/**
 * class to perform singular value decomposition on {@link DenseMatrix64F}
 */
class SVD implements Serializable {

  /**
   * 
   */
  private static final long serialVersionUID = 5517545544997813383L;
  /**
   * Column names of the original input data
   */
  private final String[] originalColNames;
  /**
   * Row names of the original input data
   */
  private final String[] originalRowNames;

  private SvdImplicitQrDecompose_D64 svd;

  /**
   * Number of components that are stored
   */
  private int numComponents;

  /**
   * @param colNames
   * @param rowNames
   */
  SVD(String[] colNames, String[] rowNames) {
    super();
    this.originalColNames = colNames;
    this.originalRowNames = rowNames;
  }

  /**
   * @param dm {@link DenseMatrix64F} that has been pre-normalized, if desired
   * @param numberOfComponentsToStore Number of components (PCs) that are stored
   * @param log
   */

  void computeSVD(DenseMatrix64F dm, int numberOfComponentsToStore, Logger log) {
    if (dm.getNumCols() != originalColNames.length) {
      throw new IllegalArgumentException("Mismatched column lengths");
    }
    if (dm.getNumRows() != originalRowNames.length) {
      throw new IllegalArgumentException("Mismatched row lengths");
    }
    log.info("Computing EJML PCs");
    this.svd = new SvdImplicitQrDecompose_D64(false, true, true, true);
    svd.decompose(dm);
    log.info("Finished Computing EJML PCs");

    //    svd.getU(null, false)

    log.info("Sorting matrices to descending order");

    SingularOps.descendingOrder(null, false, svd.getW(null), svd.getV(null, true), true);

    this.numComponents = Math.min(numberOfComponentsToStore,
                                  Math.min(svd.getW(null).numRows, svd.getV(null, true).numCols));
    if (numComponents < numberOfComponentsToStore) {
      log.info(numberOfComponentsToStore + " PCs requested, but will only be able to compute"
               + numComponents);
    }
    log.info("Subsetting matrices to " + numComponents + " components");

    svd.getV(null, true).reshape(numComponents, dm.getNumCols(), true);

  }

  /**
   * @param numComponents number of singular values to return from the W matrix
   * @return
   */
  double[] getSingularValues() {
    double[] singularValues = new double[numComponents];
    for (int i = 0; i < numComponents; i++) {
      singularValues[i] = svd.getW(null).get(i, i);
    }
    return singularValues;
  }

  private void dumpPCsToText(String file, Logger log) {
    //
    try (PrintWriter writer = new PrintWriter(new File(file))) {

      StringJoiner joiner = new StringJoiner("\t");
      joiner.add("SAMPLE");
      DenseMatrix64F v = svd.getV(null, true);
      for (int i = 0; i < v.getNumCols(); i++) {
        joiner.add("PC" + (i + 1));
      }
      writer.println(joiner.toString());

      for (int i = 0; i < v.getNumRows(); i++) {
        StringJoiner sample = new StringJoiner("\t");
        sample.add(originalColNames[i]);
        for (int j = 0; j < v.getNumCols(); j++) {
          sample.add(Double.toString(v.get(j, i)));
        }
        writer.println(sample.toString());

      }
    } catch (FileNotFoundException e) {

      log.log(Level.SEVERE, "unable to write to file " + file, e);

    }

  }

  //  private void computeLoadings(DenseMatrix64F dm) {
  //    //    Will have all markers, but not all "PCs" all the time
  //    DenseMatrix64F loadingData = new DenseMatrix64F(dm.numRows, numComponents);
  //
  //    Map<String, Integer> rowMap = new HashMap<>();
  //    for (int row = 0; row < m.getDenseMatrix().numRows; row++) {
  //      rowMap.put(m.getNameForRowIndex(row), row);
  //    }
  //
  //    Map<String, Integer> loadingMap = new HashMap<>();
  //    for (int component = 0; component < numComponents; component++) {
  //      loadingMap.put("LOADING" + component, component);
  //    }
  //    for (int row = 0; row < m.getDenseMatrix().numRows; row++) {
  //
  //      DenseMatrix64F rowData = new DenseMatrix64F(1, m.getDenseMatrix().numCols);
  //      CommonOps.extractRow(m.getDenseMatrix(), row, rowData);
  //
  //      for (int component = 0; component < numComponents; component++) {
  //        DenseMatrix64F componentData = new DenseMatrix64F(1, v.getDenseMatrix().numCols);
  //
  //        double loading = getLoading(w.getEntry(component, component), rowData.data,
  //                                    CommonOps.extractRow(v.getDenseMatrix(), component,
  //                                                         componentData).data);
  //        loadingData.add(row, component, loading);
  //      }
  //    }
  //    this.loadings = new NamedRealMatrix(rowMap, loadingMap, loadingData);
  //  }

  //  private static double getLoading(double singularValue, double[] data, double[] basis) {
  //    double loading = 0;
  //    for (int i = 0; i < basis.length; i++) {
  //      loading += data[i] * basis[i];
  //    }
  //    return loading / singularValue;
  //  }

}
