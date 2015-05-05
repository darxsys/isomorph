#include <iostream>

#include <seqan/stream.h>
#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/sequence.h>

using namespace seqan;

int main() {
	std::cout << CharString("Hello SeqAn!") << std::endl;
    return 0;
}

