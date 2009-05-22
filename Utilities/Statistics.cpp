#include "Statistics.h"

/*!
 * \brief
 * Write brief comment for Statistics here.
 * 
 * \throws <exception class>
 * Description of criteria for throwing this exception.
 * 
 * Write detailed description for Statistics here.
 * 
 * \remarks
 * Write remarks for Statistics here.
 * 
 * \see
 * Separate items with the '|' character.
 */
Statistics::Statistics(void)
{
	recording = false;
}

/*!
 * \brief
 * Write brief comment for ~Statistics here.
 * 
 * \throws <exception class>
 * Description of criteria for throwing this exception.
 * 
 * Write detailed description for ~Statistics here.
 * 
 * \remarks
 * Write remarks for ~Statistics here.
 * 
 * \see
 * Separate items with the '|' character.
 */
Statistics::~Statistics(void)
{
	if (recording) fout.close();
}


void Statistics::startRecording(const char* filename, bool append)
{
	if (recording)
	{
		fout.close();
	}
	// open the file, append if necessary
	if (append)	fout.open(filename, ios_base::app);
	else fout.open(filename, ios_base::out);
	recording = true;
}

void Statistics::stopRecording()
{
	// open the file, append if necessary
	fout.close();
	recording = false;
}

void Statistics::addEntry(int id, vector<double> &x, double fitness)
{
	// we can store everything in a big array, create the separate id
	// for local search/global search or different usage
	// or even calculate any statistics that we want
	if (!recording) return;

	// a simple thing we can do
	fout << id << " " << fitness << endl;
}
