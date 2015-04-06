/**
 * @file    GroupsExtension.h
 * @brief   Definition of GroupsExtension, the core module of groups package. 
 * @author  Akiya Jouraku
 *
 * $Id: GroupsExtension.h 12233 2010-11-26 22:44:21Z fbergmann $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/branches/libsbml-5/src/packages/groups/extension/GroupsExtension.h $
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

#ifndef GroupsExtension_h
#define GroupsExtension_h

#include <sbml/common/extern.h>
#include <sbml/SBMLTypeCodes.h>

#ifdef __cplusplus

#include <sbml/extension/SBMLExtension.h>
#include <sbml/extension/SBMLExtensionNamespaces.h>
#include <sbml/extension/SBMLExtensionRegister.h>

#include <vector>

LIBSBML_CPP_NAMESPACE_BEGIN

class LIBSBML_EXTERN GroupsExtension : public SBMLExtension
{
public:

  //---------------------------------------------------------------
  //
  // Required class methods
  //
  //---------------------------------------------------------------

  /**
   * Returns the package name of this extension.
   */
  static const std::string& getPackageName ();

  /**
   * Returns the default SBML Level this extension.
   */
  static unsigned int getDefaultLevel();

  /**
   * Returns the default SBML Version this extension.
   */
  static unsigned int getDefaultVersion();

  /**
   * Returns the default SBML version this extension.
   */
  static unsigned int getDefaultPackageVersion();

  /**
   * Returns URI of supported versions of this package.
   */
  static const std::string&  getXmlnsL3V1V1();

  //
  // Other URI needed in this package (if any)
  //

  //---------------------------------------------------------------


  /**
   * Constructor
   */
  GroupsExtension ();


  /**
   * Copy constructor.
   */
  GroupsExtension(const GroupsExtension&);


  /**
   * Destroy this object.
   */
  virtual ~GroupsExtension ();


  /**
   * Assignment operator for GroupsExtension.
   */
  GroupsExtension& operator=(const GroupsExtension&);


  /**
   * Creates and returns a deep copy of this GroupsExtension object.
   * 
   * @return a (deep) copy of this SBase object
   */
  virtual GroupsExtension* clone () const;


  /**
   * Returns the name of this package ("groups")
   *
   * @pram the name of this package ("groups")
   */
  virtual const std::string& getName() const;


  /**
   * Returns the URI (namespace) of the package corresponding to the combination of 
   * the given sbml level, sbml version, and package version.
   * Empty string will be returned if no corresponding URI exists.
   *
   * @param sbmlLevel the level of SBML
   * @param sbmlVersion the version of SBML
   * @param pkgVersion the version of package
   *
   * @return a string of the package URI
   */
  virtual const std::string& getURI(unsigned int sbmlLevel, unsigned int sbmlVersion, 
                                    unsigned int pkgVersion) const;


  /**
   * Returns the SBML level with the given URI of this package.
   *
   * @param uri the string of URI that represents one of versions of layout package
   *
   * @return the SBML level with the given URI of this package. 0 will be returned
   * if the given URI is invalid.
   *
   */
  virtual unsigned int getLevel(const std::string &uri) const;


  /**
   * Returns the SBML version with the given URI of this package.
   *
   * @param uri the string of URI that represents one of versions of layout package
   *
   * @return the SBML version with the given URI of this package. 0 will be returned
   * if the given URI is invalid.
   */
  virtual unsigned int getVersion(const std::string &uri) const;


  /**
   * Returns the package version with the given URI of this package.
   *
   * @param uri the string of URI that represents one of versions of layout package
   *
   * @return the package version with the given URI of this package. 0 will be returned
   * if the given URI is invalid.
   */
  virtual unsigned int getPackageVersion(const std::string &uri) const;


  /**
   * Returns an SBMLExtensionNamespaces<GroupsExtension> object whose alias type is 
   * LayoutPkgNamespace.
   * Null will be returned if the given uri is not defined in the layout package.
   *
   * @param uri the string of URI that represents one of versions of layout package
   *
   * @return an LayoutPkgNamespace object corresponding to the given uri. NULL will
   * be returned if the given URI is not defined in layout package.
   */
  virtual SBMLNamespaces* getSBMLExtensionNamespaces(const std::string &uri) const;


  /**
   * This method takes a type code of groups package and returns a string representing 
   * the code.
   */
  virtual const char* getStringFromTypeCode(int typeCode) const;


  /** @cond doxygen-libsbml-internal */
  /**
   * Initializes layout extension by creating an object of this class with 
   * required SBasePlugin derived objects and registering the object 
   * to the SBMLExtensionRegistry class.
   *
   * (NOTE) This function is automatically invoked when creating the following
   *        global object in GroupsExtension.cpp
   *
   *        static SBMLExtensionRegister<GroupsExtension> groupsExtensionRegistry;
   *
   */

  static void init();

  /** @endcond doxygen-libsbml-internal */

};


// --------------------------------------------------------------------
//
// Required typedef definitions 
//
// LayoutPkgNamespaces is derived from the SBMLNamespaces class and
// used when creating an object of SBase derived classes defined in
// layout package.
//
// --------------------------------------------------------------------

//
// (NOTE) 
//
// SBMLExtensionNamespaces<GroupsExtension> must be instantiated
// in GroupsExtension.cpp for DLL.
//
typedef SBMLExtensionNamespaces<GroupsExtension> GroupsPkgNamespaces; 

typedef enum
{
   SBML_GROUPS_GROUP  = 200
 , SBML_GROUPS_MEMBER = 201
} SBMLGroupsTypeCode_t;


LIBSBML_CPP_NAMESPACE_END

#endif  /* __cplusplus */
#endif  /* GroupsExtension_h */
