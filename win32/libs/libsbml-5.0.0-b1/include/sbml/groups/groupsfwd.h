/**
 * @file    groupsfwd.h
 * @brief   SBML Groups C structure declarations
 * @author  Akiya Jouraku
 *
 * $Id: groupsfwd.h 10667 2010-01-16 10:20:44Z ajouraku $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/branches/libsbml-5/src/packages/groups/common/groupsfwd.h $
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


#ifndef groupsfwd_h__
#define groupsfwd_h__


/**
 * Forward declaration of all opaque C types.
 *
 * Declaring all types up-front avoids "redefinition of type 'Foo'" compile
 * errors and allows our combined C/C++ headers to depend minimally upon
 * each other.  Put another way, the type definitions below serve the same
 * purpose as "class Foo;" forward declarations in C++ code.
 */

#ifdef __cplusplus
#  define CLASS_OR_STRUCT class
#else
#  define CLASS_OR_STRUCT struct
#endif  /* __cplusplus */

LIBSBML_CPP_NAMESPACE_BEGIN

typedef CLASS_OR_STRUCT Group   Group_t;
typedef CLASS_OR_STRUCT Member  Member_t;

LIBSBML_CPP_NAMESPACE_END

#undef CLASS_OR_STRUCT


#endif  /* groupsfwd_h__ */
