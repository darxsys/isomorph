# isomorph

Isomorph is a software for transcript abundance estimation using
RNA-Seq reads and transcripts reconstructed de novo. It is an implementation
of a model proposed [here](http://bioinformatics.oxfordjournals.org/content/26/4/493.full) and 
[here](http://www.biomedcentral.com/1471-2105/12/323). The model has been adapted and modified to
fit the needs of our own lab and will in the future be changed and rewritten.

## Dependencies

Isomorph currently relies on [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) and 
[seqan](https://github.com/seqan/seqan) open-source library.
bowtie2 has to be included in your PATH or you can provide bowtie path directly to isomorph.
Seqan should be placed inside the isomorph directory and named only seqan (without any version number).
To install it, just do `git clone` on seqan repository into the isomorph folder. 

## Installation

To install, check out the repository and run `make`. This will generate a binary which
you can then include in the `PATH`, if wanted.

To create the documentation, type `make docs`. Make sure to have [doxygen](http://www.stack.nl/~dimitri/doxygen/) 
and [graphviz](http://www.graphviz.org/) installed.

## Usage

Isomorph should be very easy to use. If for example, you have a set of single end reads along with reference
transcripts, just run the following:

		./isomorph -R reads.fq -T transcripts.fa

Detailed description of parameters is available by running 
		isomorph -h

## History

Current version: 0.2

## License

This software is licensed according to GNU General Public License.
The author cannot be considered responsible for any damages resulting from the use
of this software.