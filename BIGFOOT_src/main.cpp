// ==========================================================================
//                               	RAIDER
// ==========================================================================
// Author: Nathan Figueroa <figuernd@miamioh.edu>
// ==========================================================================

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>

#include <iostream>
#include <sstream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <unordered_map>

#include "raider.h"

using namespace std;

typedef pair<uint, seqan::CharString> idThreshold;

const string OUTPUT_SEEDS_FILENAME = "seeds";
const string OUTPUT_SUMMARY_FILENAME = "summary_info";


// --------------------------------------------------------------------------
// Class AppOptions
// --------------------------------------------------------------------------

// This struct stores the options from the command line.
struct AppOptions {
	// Verbosity level.  0 -- quiet, 1 -- normal, 2 -- verbose
	int verbosity;
	// Minimum repeat length
	uint L;
	// Minimum number of repeats to be significant
	uint count;
	// Minimum percentage required for a subsequent base to be significant
	float T;
	// Minimum percent identity required for two Lmer counts to be considered
	// close enough to possibly be in the same larger fragment
	float I;


	seqan::CharString sequence_file;
	seqan::CharString output_directory;

	AppOptions() :
			verbosity(1) {
	}
};

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------
seqan::ArgumentParser::ParseResult parseCommandLine(AppOptions & options, int argc, char const ** argv) {
	// Setup ArgumentParser.
	seqan::ArgumentParser parser("RAIDER");
	// Set short description, version, and date.
	setShortDescription(parser, "RAIDER - Rapid Ab Initio Detection of Elementary Repeats");
	setVersion(parser, "2.0");
	setDate(parser, "August 2014");

	// Define usage line and long description.
	addUsageLine(parser, "[\\fIOPTIONS\\fP] \"\\fISEQUENCE_FILE\\fP\"  \"\\fIOUTPUT_DIRECTORY\\fP\"");
	addDescription(
			parser,
			"RIADER parses the given sequence file to identify de novo repeats. Minimum repeat size and other options can be configured as described below.");

	// We require two arguments.
	addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::STRING, "SEQUENCE_FILE"));
	addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::STRING, "OUTPUT_DIRECTORY"));

	addOption(
			parser,
			seqan::ArgParseOption("l", "L", "Minimum repeat length. Defaults to 24.",
					seqan::ArgParseOption::INTEGER));
	addOption(
			parser,
			seqan::ArgParseOption("c", "count", "Minimum number of repeats in a family. Defaults to 2.",
					seqan::ArgParseOption::INTEGER));
	addOption(
				parser,
				seqan::ArgParseOption("i", "I", "Minimum percent required before an Lmer count is too low to be compatible with an Lmer of greater count. Default is 0.75",
						seqan::ArgParseOption::DOUBLE));
	addOption(
				parser,
				seqan::ArgParseOption("T", "T", "Minimum percent required for a subsequent base to be considered part of the same family. Default is 0.75.",
						seqan::ArgParseOption::DOUBLE));
	addOption(parser, seqan::ArgParseOption("q", "quiet", "Set verbosity to a minimum."));
	addOption(parser, seqan::ArgParseOption("v", "verbose", "Enable verbose output."));

	// Add Examples Section.
	addTextSection(parser, "Examples");
	addListItem(parser, "\\fBraider\\fP \\fB-v\\fP \\fIchr23.fa\\fP \"\\fIchr23_out\\fP\"",
			"Call with chr23.fa and verbose output.");

	// Parse command line.
	seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

	// Only extract  options if the program will continue after parseCommandLine()
	if (res != seqan::ArgumentParser::PARSE_OK)
		return res;

	// Extract option values.
	options.verbosity = 1;
	if (isSet(parser, "verbose"))
		options.verbosity = 2;
	if (isSet(parser, "quiet"))
		options.verbosity = 0;

	seqan::getArgumentValue(options.sequence_file, parser, 0);
	seqan::getArgumentValue(options.output_directory, parser, 1);

	if (isSet(parser, "L"))
		seqan::getOptionValue(options.L, parser, "L");
	else
		options.L = 24;

	if (isSet(parser, "count"))
		seqan::getOptionValue(options.count, parser, "count");
	else
		options.count = 2;

	if (isSet(parser, "I"))
		seqan::getOptionValue(options.I, parser, "I");
	else
		options.I = 0.75;

	if (isSet(parser, "T"))
		seqan::getOptionValue(options.I, parser, "T");
	else
		options.T = 0.75;

	// Ensure a trailing /
	if (options.output_directory[seqan::length(options.output_directory) - 1] != '/') {
		seqan::append(options.output_directory, "/");
	}

	return seqan::ArgumentParser::PARSE_OK;
}

/**
 * Given a sequence stream, read in all sequences and concatenate into a master sequence
 */
void concatenateSequences(seqan::SequenceStream &seqStream, vector<idThreshold> &thresholds,
		seqan::Dna5String &outSequence, int verbosity) {
	seqan::CharString id;
	if (seqan::readRecord(id, outSequence, seqStream) != 0) {
		throw std::runtime_error("Unable to read record");
	}
	if (verbosity > 0) {
		cout << "Preparing " << id << endl;
	}
	uint currentThreshold = seqan::length(outSequence);
	thresholds.push_back(make_pair(currentThreshold, id));
	seqan::Dna5String other;
	while (seqan::readRecord(id, other, seqStream) == 0) {
		if (verbosity > 0) {
			cout << "Preparing " << id << endl;
		}
		currentThreshold += seqan::length(other);
		thresholds.push_back(make_pair(currentThreshold, id));
		seqan::append(outSequence, other);
	}
}

/**
 * Given a FASTA name, read all sequences from file and compile into master sequence. Also return
 * list of thresholds: indices into the master sequence paired with sequence IDs to mark where each
 * individual sequence begins.
 */
void compileFasta(seqan::CharString file, seqan::Dna5String &sequence, vector<idThreshold> &thresholds, int verbosity) {
	if (verbosity > 0) {
		cout << "Loading sequence..." << endl;
	}

	seqan::SequenceStream seqStream(seqan::toCString(file));
	if (!seqan::isGood(seqStream)) {
		throw std::invalid_argument("Unknown problem with sequence stream");
	}

	concatenateSequences(seqStream, thresholds, sequence, verbosity);
}


/**
 * Print the command line arguments back to the user.
 */
void printArgs(AppOptions &options) {
	if (options.verbosity > 0) {
		cout << "__ARGUMENTS____________________________________________________________________" << endl
				<< "VERBOSITY\t" << options.verbosity << endl
				<< "MIN_LENGTH\t" << options.L << endl
				<< "MIN_COUNT\t" << options.count << endl
				<< "SEQUENCE_FILE\t" << options.sequence_file << endl
				<< "OUTPUT_DIRECTORY\t" << options.output_directory << endl;
	}
}

void writeResults(vector<pair<string, int> > &results, AppOptions &options) {
	string resultsFilePath(seqan::toCString(options.output_directory));
	resultsFilePath.append(OUTPUT_SEEDS_FILENAME);

	int maxLen = 0;
	int maxCount = 0;
	uint repCount = 0;

	if (options.verbosity > 0) {
		cout << "Writing seeds to " << resultsFilePath << endl;
	}

	ofstream seedFile;
	seedFile.open(resultsFilePath.c_str());
	for (uint i = 0; i < results.size(); i++) {
		pair<string, int> p = results[i];
		int length = p.first.length();
		seedFile <<">ID:" <<i <<"-LEN:" <<length <<"-CNT:" <<p.second <<"\n";
		seedFile <<p.first <<"\n";

		if (length > maxLen)
			maxLen = length;
		if (p.second > maxCount)
			maxCount = p.second;
		repCount += p.second;
	}
	seedFile.close();


	string summaryFilePath(seqan::toCString(options.output_directory));
	summaryFilePath.append(OUTPUT_SUMMARY_FILENAME);
	ofstream summaryFile;
	summaryFile.open(summaryFilePath.c_str());

	summaryFile << "#SEEDS\t" << results.size() << endl;
	summaryFile << "#REPEATS\t" << repCount << endl;
	summaryFile << "LONGEST_SEED\t" << maxLen << endl;
	summaryFile << "HIGHEST_COUNT\t" << maxCount << endl;
	summaryFile.close();
}


int main(int argc, char const ** argv) {
	AppOptions options;
	parseCommandLine(options, argc, argv);
	printArgs(options);

	// Load sequences into master sequence, mark where each sequence ID begins within the master
	vector<idThreshold> thresholds;
	seqan::Dna5String sequence;
	compileFasta(options.sequence_file, sequence, thresholds, options.verbosity);

	if (options.verbosity > 0) {
		cout << "BASE PAIRS\t" << seqan::length(sequence) << endl;
	}

	vector<pair<string, int> > results;
	raider(sequence, options.L, options.count, options.I, options.T, results);

	if (options.verbosity > 0) {
		cout << "Writing results..." << endl;
	}
	writeResults(results, options);

	return 0;
}
