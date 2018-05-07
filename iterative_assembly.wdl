
import "tasks/bwa.wdl" as bwa
import "tasks/samtools.wdl" as samtools
import "tasks/spades.wdl" as spades
import "tasks/biopet.wdl" as biopet


workflow iterativeAssembly {
    File inputAssembly
    File read1
    File? read2
    String outputDir

    call bwa.index as bwaIndex {
        input:
            fasta = inputAssembly
    }
    call bwa.BwaMem as bwaMem {
        input:
            inputR1 = read1,
            inputR2 = read2,
            referenceFasta = bwaIndex.indexedFasta,
            indexFiles = bwaIndex.indexFiles,
            outputPath= outputDir + "/" + "ReadsMappedToInputAssembly.bam"
    }

    call samtools.fastq as selectMappedReads {
        input:
            inputBam = bwaMem.bamFile,
            outputRead1 = outputDir + "/alignedReads/reads1.fastq.gz",
            outputRead2 = if (defined(read2)) then outputDir + "/alignedReads/read2.fastq.gz" else read2,
            excludeFilter = 516, # UNMAP,QCFAIL will be filtered out
            includeFilter = if (defined(read2)) then 2 else read2 # If paired, use reads that have PROPER_PAIR flag
    }

    # If paired. Reads need to be synced after selection.
    if (defined(read2)) {
        call biopet.FastqSync as syncSelectedReads {
            input:
                ref1 = read1,
                ref2 = select_first([read2]),
                in1 = selectMappedReads.read1,
                in2 = select_first([selectMappedReads.read2]),
                out1path = outputDir + "/mappedReads/reads1.fastq.gz",
                out2path = outputDir + "/mappedReads/reads2.fastq.gz"
        }
    }

    File selectedReads1 = select_first([syncSelectedReads.out1, selectMappedReads.read1 ])
    File selectedReads2 = select_first([syncSelectedReads.out2, selectMappedReads.read2 ])

    call spades.spades {
        input:
            read1 = selectedReads1,
            read2 = selectedReads2,
            outputDir = outputDir + "/spades"
    }

    output {
        File scaffolds = spades.scaffolds
        File contigs = spades.contigs
    }
}