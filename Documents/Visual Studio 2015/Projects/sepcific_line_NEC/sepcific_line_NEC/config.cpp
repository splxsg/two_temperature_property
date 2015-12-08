#include "config.h"

bool IsSpace(char c)//determine if it is space 
{
	if (' ' == c || '\t' == c)
		return true;
	return false;
}

bool IsCommentChar(char c)//determine if it is comment
{
	switch (c) {
	case COMMENT_CHAR:
		return true;
	default:
		return false;
	}
}

void Trim(string & str)// trim space 
{
	if (str.empty()) {
		return;
	}
	int i, start_pos, end_pos;
	for (i = 0; i < str.size(); ++i) {
		if (!IsSpace(str[i])) {
			break;
		}
	}
	if (i == str.size()) { // all space
		str = "";
		return;
	}

	start_pos = i;

	for (i = (int)str.size() - 1; i >= 0; --i) {
		if (!IsSpace(str[i])) {
			break;
		}
	}
	end_pos = i;

	str = str.substr(start_pos, end_pos - start_pos + 1);
}


string config::getValue(const string & name)
{
	string line;
	string new_line;
	infile->clear();
	infile->seekg(0);
	while (getline(*infile, line))
	{
		if (line.empty())
		{
			return "";
		}
		int start_pos = 0, end_pos = (int)line.size() - 1, pos;
		if ((pos = (int)line.find(COMMENT_CHAR)) != -1)
		{
			if (0 == pos)
			{
				return "";
			}
			end_pos = pos - 1;
		}
		new_line = line.substr(start_pos, start_pos + 1 - end_pos);  // 预处理，删除注释部分   
		if ((pos = (int)new_line.find('=')) == -1)
		{
			return "";
		}
		string na = new_line.substr(0, pos);
		Trim(na);
		if (na == name)
		{
			string value = new_line.substr(pos + 1, end_pos + 1 - (pos + 1));
			Trim(value);
			return  value;
		}
	}
	return "";
	//infile->close;
}

config::config(const string & filename)
{
	infile = new ifstream(filename.c_str());
	if (!infile)
	{
		cout << "Can't open config file" << endl;
	}

}

config::config(void)
{
}
