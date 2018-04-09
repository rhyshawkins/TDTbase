//
//    Simply Logging Library  
//    
//    Copyright (C) 2014 - 2018 Rhys Hawkins
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <errno.h>
#include <limits.h>

#include <time.h>

#include "slog.h"

static char *log_file = NULL;

static char *timestamp(void);

static const char *LEVEL_MESSAGE[] = {
  "error",
  "warning",
  "info",
  "debug",
  ""
};
  
int slog_set_output_file(const char *filename,
			 int flags)
{
  FILE *fp;

  if (log_file != NULL) {
    free(log_file);
  }

  log_file = malloc(sizeof(char) * PATH_MAX);
  if (log_file == NULL) {
    fprintf(stderr, "slog_set_output_file: failed to allocate log file buffer\n");
    return -1;
  }

  if (flags && SLOG_FLAGS_CLEAR) {
    fp = fopen(filename, "w");
  } else {
    fp = fopen(filename, "a");
  }

  if (fp == NULL) {
    free(log_file);
    log_file = NULL;
    fprintf(stderr, "log: failed to create/open log file\n");
    return -1;
  }
  
  fprintf(fp, "%s: begin\n", timestamp());
  fclose(fp);

  if (realpath(filename, log_file) == NULL) {
    fprintf(stderr, "slog_set_output_file: failed to get absolution file name\n");
    free(log_file);
    log_file = NULL;
    return -1;
  }
  
  return 0;
}

int slog(int level,
	 const char *source_file,
	 const char *function,
	 int lineno,
	 const char *msg,
	 ...)
{
  va_list ap;
  int r;
  
  va_start(ap, msg);
  r = vslog(level, source_file, function, lineno, msg, ap);
  va_end(ap);

  return r;
}

int vslog(int level,
	  const char *source_file,
	  const char *function,
	  int lineno,
	  const char *msg,
	  va_list ap)
{
  FILE *fp;
  
  if (log_file != NULL) {
    fp = fopen(log_file, "a");
    if (fp == NULL) {
      fprintf(stderr, "slog: failed to open log file in append mode %s, %s\n", log_file, strerror(errno));
      return -1;
    }
  } else {
    fp = stderr;
  }

  fprintf(fp, "%s:%s:%s:%s:%4d:",
	  timestamp(),
	  LEVEL_MESSAGE[level],
	  source_file,
	  function,
	  lineno);

  vfprintf(fp, msg, ap);
  fprintf(fp, "\n");
  
  if (log_file != NULL) {
    fclose(fp);
  }

  return 0;
}

static FILE *slogopen = NULL;

FILE *slogstart(const char *source_file,
		 const char *function,
		 int lineno)
{
  if (slogopen == NULL) {
    if (log_file != NULL) {
      slogopen = fopen(log_file, "a");
      if (slogopen == NULL) {
	fprintf(stderr, "slog: failed to open log file in append mode %s, %s\n", log_file, strerror(errno));
	return NULL;
      }
    } else {
      slogopen = stderr;
    }
  }

  fprintf(slogopen,
	  "%s:%s:%s:%s:%4d:\n",
	  timestamp(),
	  "LONG",
	  source_file,
	  function,
	  lineno);

  return slogopen;
}

void slogend()
{
  if (slogopen != NULL) {

    fprintf(slogopen, "\n");
    
    if (slogopen != stderr) {
      fclose(slogopen);
    }
    slogopen = NULL;
  }
}
  
		   


static FILE *large_message_fp = NULL;

int slog_large_message_start(int level,
			     const char *source_file,
			     const char *function,
			     int lineno,
			     const char *msg,
			     ...)
{
  va_list ap;
    
  if (log_file != NULL) {
    large_message_fp = fopen(log_file, "a");
    if (large_message_fp == NULL) {
      fprintf(stderr, "log: failed to open log file in append mode\n");
      return -1;
    }
  } else {
    large_message_fp = stderr;
  }

  fprintf(large_message_fp, "%s:%s:%s:%s:%4d:",
	  timestamp(),
	  LEVEL_MESSAGE[level],
	  source_file,
	  function,
	  lineno);

  va_start(ap, msg);
  vfprintf(large_message_fp, msg, ap);
  va_end(ap);

  fprintf(large_message_fp, "\n  ");

  return 0;
}

int slog_large_message_write(const char *msg, ...)
{
  va_list ap;
  
  if (large_message_fp == NULL) {
    fprintf(stderr, "slog_large_message_write: unassigned output\n");
    return -1;
  }
  
  va_start(ap, msg);
  vfprintf(large_message_fp, msg, ap);
  va_end(ap);

  return 0;
}

int slog_large_message_newline(void)
{
  if (large_message_fp == NULL) {
    fprintf(stderr, "slog_large_message_write: unassigned output\n");
    return -1;
  }

  fprintf(large_message_fp, "\n  ");
  return 0;
}

int slog_large_message_end(void)
{
  if (large_message_fp == NULL) {
    fprintf(stderr, "slog_large_message_write: unassigned output\n");
    return -1;
  }

  fprintf(large_message_fp, "\n");

  if (large_message_fp != stderr) {
    fclose(large_message_fp);
  }
  large_message_fp = NULL;
  return 0;
}

static char *timestamp(void)
{
  const char *TIME_FORMAT = "%Y-%m-%d %H:%M:%S";
  static char buffer[256];
  time_t tmp;
  struct tm *t;

  tmp = time(NULL);
  t = localtime(&tmp);
  
  if (strftime(buffer, sizeof(buffer), TIME_FORMAT, t) == 0) {
    fprintf(stderr, "log::timestamp: failed to format time\n");
    return NULL;
  }

  return buffer;
}


