/**
 * @file    SBMLExtensionRegister.h
 * @brief   Definition of SBMLExtensionRegister, the template class for registering
 *          an extension package to SBMLExtensionRegistry class.
 * @author  Akiya Jouraku
 *
 * $Id: SBMLExtensionRegister.h 10667 2010-01-16 10:20:44Z ajouraku $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/branches/libsbml-5/src/extension/SBMLExtensionRegister.h $
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
 */

#ifndef SBMLExtensionRegister_h
#define SBMLExtensionRegister_h

#include <sbml/extension/SBMLExtension.h>
#include <sbml/extension/SBMLExtensionRegistry.h>

#ifdef __cplusplus

LIBSBML_CPP_NAMESPACE_BEGIN

template<class SBMLExtensionType>
class LIBSBML_EXTERN SBMLExtensionRegister
{
public:

  /**
   * Constructor
   *
   * Initialization code of corresponding package extension 
   * will be executed when an object of this class is created.
   *
   */
  SBMLExtensionRegister() { SBMLExtensionType::init(); };

};

LIBSBML_CPP_NAMESPACE_END

#endif  /* __cplusplus */
#endif  /* SBMLExtensionRegister_h */
