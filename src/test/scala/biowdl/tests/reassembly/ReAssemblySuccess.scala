/*
 * Copyright (c) 2018 Biowdl
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of
 * this software and associated documentation files (the "Software"), to deal in
 * the Software without restriction, including without limitation the rights to
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
 * the Software, and to permit persons to whom the Software is furnished to do so,
 * subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
 * FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
 * COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
 * IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package biowdl.tests.reassembly

import java.io.File

import nl.biopet.utils.biowdl.PipelineSuccess

trait ReAssemblySuccess extends ReAssembly with PipelineSuccess {
  val unmapppedReadsFilteredFile: File = createFile(
    "unmapppedReadsFiltered.bam")
  val readsMappedToInputAssemblyFile: File = createFile(
    "ReadsMappedToInputAssembly.bam")
  val filteredR1File: File = createFile("filtered_reads1.fq.gz")
  val filteredR2File: Option[File] =
    createOptionalFile(read2.isDefined, "filtered_reads2.fq.gz")

  val scaffoldsFastaFile: File = createFile("spades", "scaffolds.fasta")
  val contigsFastFile: File = createFile("spades", "contigs.fasta")
}
