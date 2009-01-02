/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2008 by <Quang Huy / NTU>
 */

#ifndef _Scaling_Linear_H_
#define _Scaling_Linear_H_

#include <vector>
#include "Scaling.h"

using namespace std;

/*!
 * \brief
 * Implementation of the linear scaling strategy. 
 *
 * \see
 * Scaling
 */

class Scaling_Linear : public Scaling {
public:	
	virtual bool operator() (vector<double> & v);
	
private:	
};

#endif
