#ifndef _CONFIG_H_
#define _CONFIG_H_

#include "parameter.h"
#define COMMENT_CHAR '#'//comment sign

class config
{
private:
	ifstream *infile;
public:
	config(void);
	config(const string & filename);
	string getValue(const string & name); //config name, return value
};

#endif