/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2005 by <your name/ organization here>
 */
#pragma once

#ifndef _Statistics_H_
#define _Statistics_H_

#include <fstream>
#include <vector>
using namespace std;

class Statistics
{
public:
	Statistics(void);
	~Statistics(void);

	/*!
	 * \brief
	 * Write brief comment for addEntry here.
	 * 
	 * \param x
	 * Description of parameter x.
	 * 
	 * \param fitness
	 * Description of parameter fitness.
	 * 
	 * \throws <exception class>
	 * Description of criteria for throwing this exception.
	 * 
	 * Write detailed description for addEntry here.
	 * 
	 * \remarks
	 * Write remarks for addEntry here.
	 * 
	 * \see
	 * Separate items with the '|' character.
	 */
	virtual void addEntry(int id, vector<double>& x, double fitness);

	/*!
	 * \brief
	 * Write brief comment for startRecording here.
	 * 
	 * \throws <exception class>
	 * Description of criteria for throwing this exception.
	 * 
	 * Write detailed description for startRecording here.
	 * 
	 * \remarks
	 * Write remarks for startRecording here.
	 * 
	 * \see
	 * Separate items with the '|' character.
	 */
	virtual void startRecording(const char* filename, bool append = false);
	
	/*!
	 * \brief
	 * Write brief comment for stopRecording here.
	 * 
	 * \throws <exception class>
	 * Description of criteria for throwing this exception.
	 * 
	 * Write detailed description for stopRecording here.
	 * 
	 * \remarks
	 * Write remarks for stopRecording here.
	 * 
	 * \see
	 * Separate items with the '|' character.
	 */
	virtual void stopRecording();

private:
	bool recording;	
	ofstream fout;

};

#endif