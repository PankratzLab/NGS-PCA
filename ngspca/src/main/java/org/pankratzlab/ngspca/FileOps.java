package org.pankratzlab.ngspca;

import java.io.File;
import java.util.List;
import java.util.stream.Collectors;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;

public class FileOps {

  static String stripDirectoryAndExtension(String path, String extension) {

    return FilenameUtils.getName(path).replaceAll(extension, "");
  }

  static boolean fileExists(String path) {
    File f = new File(path);
    return f.exists() && !f.isDirectory();
  }

  static List<String> listFilesWithExtension(String inputDir, String[] extensions) {
    return FileUtils.listFiles(new File(inputDir), extensions, true).stream()
                    .map(File::getAbsolutePath).collect(Collectors.toList());
  }
}
