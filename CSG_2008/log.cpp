
#include "prelude.h"

#include <cstdlib>
using std::atexit;

#include <ctime>
#include <fstream>
using std::ofstream;
using std::endl;

namespace {

	ofstream error_log_stream;

	void on_exit()
	{
		error_log_stream << "Ending error logging at " << endl;
		std::time_t time_var = std::time(NULL);
		error_log_stream << std::ctime(&time_var) << endl;
		error_log_stream.close();
	}

}

void logInit()
{
	error_log_stream.open("error_log.txt", std::ios_base::app | std::ios_base::out);
	error_log_stream << "Begining error logging at " << endl;
	std::time_t time_var = std::time(NULL);
	error_log_stream << std::ctime(&time_var);// << endl;

	atexit(on_exit);
}

// exposes the error log for writing
std::ostream& err()
{
	static bool first = true;
	if (first) {
		logInit();
		first = false;
	}
	return error_log_stream;
}
