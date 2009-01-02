/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2008 by <Quang Huy / NTU>
 */
#ifndef _Scaling_H_
#define _Scaling_H_

#include <vector>
using namespace std;

/*!
 * \brief
 * The super class for the scaling operators.
 *  
 * \remark
 *
 * To be used in Selection operators.
 *
 * \see
 * Scaling_Linear | Selection
 */

class Scaling {

public:	

	/*!
	 * \brief
	 * Wrapper operator for the scale function.
	 *	 
	 */
	virtual bool operator()(vector< double > & v)
	{
		return scale(v);
	}

	/*!
	 * \brief
	 * Implementation of the scaling strategy.
	 * 	
	 */
	virtual bool scale(vector< double > & v)
	{
		return false;
	}
	
private:	
};

#endif
