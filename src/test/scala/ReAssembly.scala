import java.io.File

import nl.biopet.utils.biowdl.Pipeline

trait ReAssembly extends Pipeline {
  def inputAssembly: File

  def read1: File

  def read2: Option[File]

  def startFile: File = new File("ReAssembly.wdl")

  override def inputs: Map[String, Any] =
    super.inputs ++
      Map(
        "ReAssembly.read1" -> read1.getAbsolutePath,
        "ReAssembly.inputAssembly" -> inputAssembly.getAbsolutePath,
        "ReAssembly.outputDir" -> outputDir.getAbsolutePath
      ) ++
      read2.map("ReAssembly.read2" -> _.getAbsolutePath)

}
