
import "tasks/bwa.wdl" as bwa
import "tasks/samtools.wdl" as samtools
import "tasks/spades.wdl" as spades


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
            includeFilter = if (defined(read2)) then 2 else read2 # If paired, only properly paired reads should be used.
    }

    call spades.spades {
        input:
            read1 = selectMappedReads.read1,
            read2 = selectMappedReads.read2,
            outputDir = outputDir + "/spades"
    }

    output {
        File scaffolds = spades.scaffolds
        File contigs = spades.contigs
    }
}