package org.pankratzlab.ngspca;

import java.io.Serializable;
import java.util.logging.Logger;

import org.apache.commons.math3.linear.DiagonalMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.ejml.alg.dense.decomposition.svd.SvdImplicitQrDecompose_D64;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.SingularOps;

/**
 * class to perform singular value decomposition on {@link DenseMatrix64F}
 *
 */
class SVD implements Serializable {

	/**
	 * 
	 */
	private static final long serialVersionUID = 5517545544997813383L;
	private final DenseMatrix64F dm;
	private final String[] colNames;
	private final String[] rowNames;

	// With V,W, and original data M we can always compute U
	private DenseMatrix64F v;
	private double[] w;

	/**
	 * @param dm
	 *            {@link DenseMatrix64F} that has been pre-normalized, if desired
	 * @param colNames
	 * @param rowNames
	 */
	SVD(DenseMatrix64F dm, String[] colNames, String[] rowNames) {
		super();
		this.dm = dm;
		this.colNames = colNames;
		this.rowNames = rowNames;
		if (dm.getNumCols() != colNames.length) {
			throw new IllegalArgumentException("Mismatched column lengths");
		}
		if (dm.getNumRows() != rowNames.length) {
			throw new IllegalArgumentException("Mismatched row lengths");
		}
	}

	void computeSVD(Logger log) {

		log.info("Computing EJML PCs");
		SvdImplicitQrDecompose_D64 svd = new SvdImplicitQrDecompose_D64(false, false, true, false);
		svd.decompose(dm);
		log.info("Finished Computing EJML PCs");

		this.v = svd.getV(null, true);

		DenseMatrix64F tmpW = svd.getW(null);
		SingularOps.descendingOrder(null, false, tmpW, v, true);
		int numSingular = Math.min(tmpW.numRows, tmpW.numCols);
		this.w = new double[numSingular];
		for (int i = 0; i < numSingular; i++) {
			w[i] = tmpW.get(i, i);
		}
	}
}
