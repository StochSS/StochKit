/**
 * @file    SBMLExtensionException.h
 * @brief   SBMLExtensionExceptions class, the exception class for package extension
 * @author  Akiya Jouraku
 *
 * $Id: $
 * $HeadURL: $
 *
 *<!---------------------------------------------------------------------------
 * This file is part of libSBML.  Please visit http://sbml.org for more
 * information about SBML, and the latest version of libSBML.
 *
 * Copyright 2009 California Institute of Technology.
 * 
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation.  A copy of the license agreement is provided
 * in the file named "LICENSE.txt" included with this software distribution
 * and also available online as http://sbml.org/software/libsbml/license.html
 *------------------------------------------------------------------------- -->
 *
 */

#ifndef SBMLExtensionException_h
#define SBMLExtensionException_h

#include <sbml/common/common.h>

#ifdef __cplusplus

#include <string>
#include <stdexcept>

LIBSBML_CPP_NAMESPACE_BEGIN

/**
 * An exception class for SBMLExtensionException.
 */
class LIBSBML_EXTERN SBMLExtensionException : public std::invalid_argument
{
public:

  /** 
   * constructor 
   */
  SBMLExtensionException (const std::string& errmsg) throw();

#ifndef SWIG
  virtual ~SBMLExtensionException () throw();
#endif
};

LIBSBML_CPP_NAMESPACE_END

#endif  /* __cplusplus */
#endif  /* SBMLExtensionException_h */
