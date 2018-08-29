import java.io.File
import nl.biopet.utils.biowdl.fixtureFile

class ReAssemblyTest extends ReAssemblySuccess {
  def inputAssembly: File = fixtureFile("reference", "reference.fasta")

  def read1: File = fixtureFile("samples", "wgs1", "R1.fq.gz")

  def read2: Option[File] = Some(fixtureFile("samples", "wgs1", "R2.fq.gz"))
}
