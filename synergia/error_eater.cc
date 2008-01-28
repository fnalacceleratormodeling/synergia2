#include <iostream>
#include <sstream>
#include <streambuf>
#undef _POSIX_C_SOURCE
#include <boost/python.hpp>

using namespace boost::python;

class Error_eater {
private:
  std::stringstream* error_eater;
  std::streambuf* real_error_buffer;
  bool eating;
public:
  Error_eater()
  {
    eating = false;
  }
  void start()
  {
    real_error_buffer = std::cerr.rdbuf();
    error_eater = new std::stringstream;
    std::cerr.rdbuf(error_eater->rdbuf());
    eating = true;
  }
  void stop()
  {
    if(eating){
      std::cerr.rdbuf(real_error_buffer);
      delete error_eater;
    }
    eating = false;
  }
  void dump()
  {
    if(eating){
      std::cerr.rdbuf(real_error_buffer);
      std::cerr << error_eater->str() << std::endl;
      stop();
      start();
    }
  }
  bool is_eating()
  {
    return eating;
  }
  ~Error_eater()
  {
    if (eating) {
      delete error_eater;
    }
  }
};

void test()
{
  std::cout << "message sent to stdout\n";
  std::cerr << "message sent to stderr\n";
}

BOOST_PYTHON_MODULE(error_eater)
{
  def("test",&test);
  class_<Error_eater>("Error_eater",init<>() )
    .def("start",&Error_eater::start)
    .def("stop",&Error_eater::stop)
    .def("dump",&Error_eater::dump)
    .def("is_eating",&Error_eater::is_eating);
}

