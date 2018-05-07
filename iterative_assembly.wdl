
import "tasks/bwa.wdl" as bwa
import "tasks/samtools.wdl" as samtools
import "tasks/spades.wdl" as spades
import "tasks/biopet.wdl" as biopet

# This workflow takes an existing assembly and reads that have already passed QC.
# It maps the reads back to the assembly and uses the mapped reads to create a new assembly.
workflow iterativeAssembly {
    File inputAssembly
    File read1
    File? read2
    String outputDir

    # First index the assembly
    call bwa.index as bwaIndex {
        input:
            fasta = inputAssembly
    }

    # Map the reads back to the assembly.
    call bwa.BwaMem as bwaMem {
        input:
            inputR1 = read1,
            inputR2 = read2,
            referenceFasta = bwaIndex.indexedFasta,
            indexFiles = bwaIndex.indexFiles,
            outputPath= outputDir + "/" + "ReadsMappedToInputAssembly.bam"
    }

    # Get the reads that mapped to the assembly. This means filtering out the UNMAP flag.
    # QCFAIL is also filtered out. For paired reads, the PROPER_PAIR flag is used.
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

    # This will default to the single mapped reads if there is no read2
    File selectedReads1 = select_first([syncSelectedReads.out1, selectMappedReads.read1 ])
    # This will default to none when there is no read2
    File selectedReads2 = select_first([syncSelectedReads.out2, read2 ])

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