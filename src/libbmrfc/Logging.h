/*
 * Logging.h
 *
 *  Created on: May 31, 2012
 *      Author: barnesrobe
 */

#ifndef LOGGING_H_
#define LOGGING_H_

// Logging is available
#include <log4cxx/logger.h>
#include <log4cxx/stream.h>
// provide access to configure common appenders
#include <log4cxx/defaultconfigurator.h>
#include <log4cxx/basicconfigurator.h>
#include <log4cxx/propertyconfigurator.h>

#define LOG4CXX_DECLARE(logger)      static log4cxx::LoggerPtr logger;
#define LOG4CXX_INIT(logger, lname)  log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger(lname));
#define LOG4CXX_CONF()  log4cxx::BasicConfigurator::configure();

#include <sstream>
#define LOG4CXX_TRACE_STREAM(logger, args) { \
   std::stringstream sst; \
   sst << args; \
   LOG4CXX_TRACE(logger, sst.str()); \
}
#define LOG4CXX_DEBUG_STREAM(logger, args) { \
   std::stringstream sst; \
   sst << args; \
   LOG4CXX_DEBUG(logger, sst.str()); \
}
#define LOG4CXX_INFO_STREAM(logger, args) { \
   std::stringstream sst; \
   sst << args; \
   LOG4CXX_INFO(logger, sst.str()); \
}
#define LOG4CXX_WARN_STREAM(logger, args) { \
   std::stringstream sst; \
   sst << args; \
   LOG4CXX_WARN(logger, sst.str()); \
}
#define LOG4CXX_ERROR_STREAM(logger, args) { \
   std::stringstream sst; \
   sst << args; \
   LOG4CXX_ERROR(logger, sst.str()); \
}
#define LOG4CXX_FATAL_STREAM(logger, args) { \
   std::stringstream sst; \
   sst << args; \
   LOG4CXX_FATAL(logger, sst.str()); \
}

#endif /* LOGGING_H_ */
