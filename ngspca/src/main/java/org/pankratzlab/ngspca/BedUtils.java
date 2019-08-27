package org.pankratzlab.ngspca;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;
import java.util.stream.Collectors;
import org.apache.commons.lang3.StringUtils;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.OverlapDetector;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.FeatureReader;
import htsjdk.tribble.TribbleException;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;

public class BedUtils {

  private BedUtils() {

  }

  /**
   * @param file load autosomal {@link BEDFeature}s from this file and convert to UCSC format
   * @param excluder autosomal regions that return true for
   *          {@link BEDOverlapDetector#overlapsAny(Locatable)} will not be included
   * @return {@link List} of ucsc regions
   */
  static List<String> loadAutosomalUCSC(String file, BEDOverlapDetector excluder) {
    return BedUtils.loadAll(file).stream().filter(BedUtils::autosomal)
                   .filter(excluder::overlapsNone).map(BedUtils::getBedUCSC)
                   .collect(Collectors.toList());
  }

  /**
   * @param file load all {@link BEDFeature}s in this file
   * @return {@link List} of {@link BEDFeature}s
   */
  static List<BEDFeature> loadAll(String file) {
    BEDFileReader reader = new BEDFileReader(file, false);
    List<BEDFeature> result = reader.iterator().stream().collect(Collectors.toList());
    reader.close();
    return result;

  }

  /**
   * TODO, currently not an index-based query. Reads all records and filters
   * 
   * @param file load {@link BEDFeature}s from this file
   * @param ucscRegions {@link Set} of ucsc formatted regions to load
   * @return
   */
  static BedRegionResult loadSpecificRegions(String file, Set<String> ucscRegions) {
    BEDFileReader reader = new BEDFileReader(file, false);
    CloseableIterator<BEDFeature> iter = reader.iterator();
    List<BEDFeature> result = iter.stream().filter(bf -> ucscRegions.contains(getBedUCSC(bf)))
                                  .collect(Collectors.toList());
    iter.close();
    reader.close();

    return new BedRegionResult(file, result);

  }

  /**
   * Stores the results of loading a bed file to memory
   */
  static class BedRegionResult {

    final String file;
    final List<BEDFeature> features;

    /**
     * @param file
     * @param features
     */
    private BedRegionResult(String file, List<BEDFeature> features) {
      super();
      this.file = file;
      this.features = features;
    }

  }

  /**
   * @param bedFeature
   * @return UCSC representation of this {@link BEDFeature}
   */
  static String getBedUCSC(BEDFeature bedFeature) {
    return bedFeature.getContig() + ":" + bedFeature.getStart() + "-" + bedFeature.getEnd();

  }

  private static boolean autosomal(BEDFeature bedFeature) {
    String contig = bedFeature.getContig().replaceAll("chr", "");
    return StringUtils.isNumeric(contig) && Integer.parseInt(contig) < 23;
  }

  /**
   * @param file bed file
   * @param requireIndex , if the .tbi index is required for query
   * @return {@link BEDFileReader}
   */
  static BEDFileReader getReader(String file, boolean requireIndex) {
    return new BEDFileReader(file, requireIndex);
  }

  static class BEDFileReader implements Iterable<BEDFeature> {

    private final FeatureReader<BEDFeature> reader;

    /**
     * TODO : currently only handles .tbi indices
     */
    private BEDFileReader(final String file, final boolean requireIndex) {
      reader = AbstractFeatureReader.getFeatureReader(file, new BEDCodec(), requireIndex);
    }

    // /** Queries for records within the region specified. */
    // public CloseableIterator<BEDFeature> query(final String chrom, final int
    // start, final int end) {
    // try {
    // return reader.query(chrom, start, end);
    // } catch (final IOException ioe) {
    // throw new TribbleException("Could not create an iterator from a feature
    // reader.", ioe);
    // }
    // }

    public void close() {
      try {
        reader.close();
      } catch (final IOException ioe) {
        throw new TribbleException("Could not close a bed context feature reader.", ioe);
      }
    }

    /** Returns an iterator over all records in this bed file */
    public CloseableIterator<BEDFeature> iterator() {
      try {
        return reader.iterator();
      } catch (final IOException ioe) {
        throw new TribbleException("Could not create an iterator from a feature reader.", ioe);
      }
    }

  }

  /**
   * Helper to test for overlapping regions sourced from a bed file
   */
  static class BEDOverlapDetector {

    private OverlapDetector<BEDFeature> detector;
    private int numOverlapped;

    /**
     * @param bedFile
     */
    BEDOverlapDetector(String bedFile, Logger log) {
      super();
      numOverlapped = 0;
      if (bedFile != null) {
        if (!FileOps.fileExists(bedFile)) {
          String err = "BED file to exclude did not exist: " + bedFile;
          log.severe(err);
          throw new IllegalArgumentException(err);
        }
        List<BEDFeature> regions = loadAll(bedFile);
        log.info("Loaded " + regions.size() + " regions to detect overlaps from " + bedFile);
        detector = OverlapDetector.create(loadAll(bedFile));
      } else {
        detector = OverlapDetector.create(new ArrayList<>());
      }
    }

    boolean overlapsAny(Locatable query) {
      boolean olap = detector.overlapsAny(query);
      if (olap) {
        numOverlapped++;
      }
      return olap;
    }

    boolean overlapsNone(Locatable query) {
      return !overlapsAny(query);
    }

    /**
     * @return the numExcluded
     */
    int getNumExcluded() {
      return numOverlapped;
    }

  }
}
