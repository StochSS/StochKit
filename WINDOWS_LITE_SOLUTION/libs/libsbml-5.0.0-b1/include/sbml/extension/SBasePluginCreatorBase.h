/**
 * @file    SBasePluginCreatorBase.h
 * @brief   Definition of SBasePluginCreatorBase, the base class of 
 *          SBasePlugin creator classes.
 * @author  Akiya Jouraku
 *
 * $Id: SBasePluginCreatorBase.h 10667 2010-01-16 10:20:44Z ajouraku $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/branches/libsbml-5/src/extension/SBasePluginCreatorBase.h $
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

#ifndef SBasePluginCreatorBase_h
#define SBasePluginCreatorBase_h


#include <sbml/SBMLDocument.h>
#include <sbml/SBMLNamespaces.h>
#include <sbml/extension/SBaseExtensionPoint.h>

#ifdef __cplusplus

LIBSBML_CPP_NAMESPACE_BEGIN

class SBasePlugin;

class LIBSBML_EXTERN SBasePluginCreatorBase
{
public:

  typedef std::vector<std::string>           SupportedPackageURIList;
  typedef std::vector<std::string>::iterator SupportedPackageURIListIter;

  /**
   * Destructor
   */
  virtual ~SBasePluginCreatorBase ();


  /**
   * Creates an SBasePlugin with the given uri and the prefix
   * of the target pacakge extension.
   */
  virtual SBasePlugin* createPlugin(const std::string& uri, 
                                    const std::string& prefix,
                                    const XMLNamespaces *xmlns) const = 0;


  /**
   * clone
   */
  virtual SBasePluginCreatorBase* clone() const = 0;


  /**
   * Returns the number of supported packages by this creator object.
   */
  unsigned int getNumOfSupportedPackageURI() const;


  /**
   * Returns the supported package with the given index.
   */
  std::string getSupportedPackageURI(unsigned int) const;


  /**
   * Get an SBMLTypeCode tied with this creator object.
   */
  int getTargetSBMLTypeCode() const;


  /**
   * Get an SBMLTypeCode tied with this creator object.
   */
  const std::string& getTargetPackageName() const;


  /**
   * Get an SBaseExtensionPoint tied with this creator object.
   */
  const SBaseExtensionPoint& getTargetExtensionPoint() const;


  /**
   * Returns true if a package with the given namespace is supported.
   */
  bool isSupported(const std::string& uri) const;

protected:

  /**
   * Constructor
   */
  SBasePluginCreatorBase (const SBaseExtensionPoint& extPoint,
                          const std::vector<std::string>&);


  /**
   * Copy Constructor
   */
  SBasePluginCreatorBase (const SBasePluginCreatorBase&);

  /** @cond doxygen-libsbml-internal */

  SupportedPackageURIList  mSupportedPackageURI;
  SBaseExtensionPoint       mTargetExtensionPoint;

  /** @endcond doxygen-libsbml-internal */


private:
  /** @cond doxygen-libsbml-internal */
  
  SBasePluginCreatorBase& operator=(const SBasePluginCreatorBase&);

  /** @endcond doxygen-libsbml-internal */
};

LIBSBML_CPP_NAMESPACE_END

#endif  /* __cplusplus */
#endif  /* SBasePluginCreatorBase_h */
