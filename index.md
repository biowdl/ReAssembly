---
layout: default
title: Home
version: develop
latest: true
---

A pipeline that tries to improve your assembly.
1. Takes your current assembly
2. Indexes it
3. Maps reads back to it
4. Extracts the reads from the resulting bam file
5. Constructs a new assembly using the extracted reads.


## Usage
In order to run the complete multisample pipeline, you can
run `ReAssembly.wdl` using
[Cromwell](http://cromwell.readthedocs.io/en/stable/):
```bash
java -jar cromwell-<version>.jar run -i inputs.json ReAssembly.wdl
```

The inputs JSON can be generated using WOMtools as described in the [WOMtools
documentation](http://cromwell.readthedocs.io/en/stable/WOMtool/). Note that
not some inputs should not be used! See [this page](inputs.md) for more
information.

The primary inputs are described below, additional inputs (such as precommands
and JAR paths) are available. Please use the above mentioned WOMtools command
to see all available inputs.

| field | type | default | |
|-|-|-|-|
| inputAssembly | `File` | | The assembly which is to be reassembled. |
| outputDir | `String` | | The output directory. |
| read1 | `File` | | The first end FASTQ files. |
| read2 | `File?` | | The second-end FASTQ files. |

>All inputs have to be preceded by `ReAssembly.`.
Type is indicated according to the WDL data types: `File` should be indicators
of file location (a string in JSON). Types ending in `?` indicate the input is
optional, types ending in `+` indicate they require at least one element.

## Output
This pipeline will produce a new assembly, consisting of a scaffolds file and
contigs file.

## Contact
<p>
  <!-- Obscure e-mail address for spammers -->
For any question related to this pipeline, please use the
<a href='https://github.com/biowdl/virus-assembly/issues'>github issue tracker</a>
or contact
 <a href='http://sasc.lumc.nl/'>the SASC team</a> directly at: <a href='&#109;&#97;&#105;&#108;&#116;&#111;&#58;&#115;&#97;&#115;&#99;&#64;&#108;&#117;&#109;&#99;&#46;&#110;&#108;'>
&#115;&#97;&#115;&#99;&#64;&#108;&#117;&#109;&#99;&#46;&#110;&#108;</a>.
</p>
