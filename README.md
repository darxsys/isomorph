# isomorph

Isomorph is a software for transcript abundance estimation using
RNA-Seq reads and transcripts reconstructed de novo. It is an implementation
of a model proposed [here](http://bioinformatics.oxfordjournals.org/content/26/4/493.full) and 
[here](http://www.biomedcentral.com/1471-2105/12/323). The model has been adapted and modified to
fit the needs of our own lab and will in the future be changed and rewritten.

## Installation

To install, check out the repository and run make. This will generate a binary which
you can then be included in the PATH, if wanted.

To create the documentation, type make docs. Make sure to have doxygen and graphviz installed.

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