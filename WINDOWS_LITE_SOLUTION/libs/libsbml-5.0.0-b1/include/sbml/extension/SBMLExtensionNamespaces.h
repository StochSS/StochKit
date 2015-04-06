/**
 * @file    SBMLExtensionNamespaces.h
 * @brief   SBMLExtensionNamespaces class to store level/version and namespace of
 *          SBML extension package
 * @author  Akiya Jouraku
 *
 * $Id: SBMLExtensionNamespaces.h 10156 2009-09-01 12:09:35Z sarahkeating $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/sbml/SBMLExtensionNamespaces.h $
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
 * @class SBMLExtensionNamespaces
 * @brief Class to store level, version and namespace information of SBML extension
 *        package.
 *
 */

#ifndef SBMLExtensionNamespaces_h
#define SBMLExtensionNamespaces_h

#include <sbml/SBMLNamespaces.h>
#include <sbml/common/common.h>
#include <sbml/extension/SBMLExtensionRegistry.h>
#include <sbml/extension/SBMLExtensionException.h>

#ifdef __cplusplus

#include <string>
#include <stdexcept>

LIBSBML_CPP_NAMESPACE_BEGIN

template<class SBMLExtensionType>
class LIBSBML_EXTERN SBMLExtensionNamespaces : public SBMLNamespaces
{
public:

  /**
   * Creates a new SBMLExtensionNamespaces object corresponding to the given SBML
   * @p level, @p version and @p package version.
   *
   * @note SBMLExtensionException will be thrown if the extension module
   *       that supports the combination of the given sbml level, sbml version, 
   *       package name, and package version has not been registered.
   * 
   * @param level   the SBML level
   * @param version the SBML version
   * @param pkgVersion the package version
   * @param prefix  the prefix of the package namespace (e.g. "layout", "multi") 
   *        to be added. The package's name will be used if the given string is empty 
   *        (default).
   */
  SBMLExtensionNamespaces(unsigned int level        = SBMLExtensionType::getDefaultLevel(), 
                          unsigned int version      = SBMLExtensionType::getDefaultVersion(), 
                          unsigned int pkgVersion   = SBMLExtensionType::getDefaultPackageVersion(), 
                          const std::string& prefix = SBMLExtensionType::getPackageName()) 
#ifndef SWIG
    : SBMLNamespaces(level, version, SBMLExtensionType::getPackageName(), pkgVersion, prefix)
     ,mPackageVersion(pkgVersion), mPackageName(prefix)
   {}
#else
   ;
#endif //SWIG


  /**
   * Destroys this SBMLExtensionNamespaces object.
   */
  virtual ~SBMLExtensionNamespaces() 
#ifndef SWIG
  {}
#else
  ;
#endif //SWIG

  
  /**
   * Copy constructor; creates a copy of a SBMLExtensionNamespaces.
   * 
   * @param orig the SBMLExtensionNamespaces instance to copy.
   */
  SBMLExtensionNamespaces(const SBMLExtensionNamespaces& orig)
#ifndef SWIG
   : SBMLNamespaces(orig)
    ,mPackageVersion(orig.mPackageVersion), mPackageName(orig.mPackageName)
  {}
#else
  ;
#endif //SWIG


  /**
   * Assignment operator for SBMLExtensionNamespaces.
   */
  SBMLExtensionNamespaces& operator=(const SBMLExtensionNamespaces& orig)
#ifndef SWIG
  {
    SBMLNamespaces::operator=(orig);
    mPackageVersion = orig.mPackageVersion;

    return *this;
  }
#else
  ;
#endif //SWIG


  /**
   * Creates and returns a deep copy of this SBMLExtensionNamespaces.
   * 
   * @return a (deep) copy of this SBMLExtensionNamespaces.
   */
  virtual SBMLExtensionNamespaces* clone () const
#ifndef SWIG
  {
    return new SBMLExtensionNamespaces(*this);
  }
#else
  ;
#endif //SWIG


  /**
   * Returns a string representing the Package XML namespace of this
   * object.
   *
   * @return a string representing the SBML namespace that reflects the
   * SBML Level and Version of this object.
   */
  virtual std::string getURI() const
#ifndef SWIG
  {
    const SBMLExtension *sbext = SBMLExtensionRegistry::getInstance().getExtensionInternal(SBMLExtensionType::getPackageName());
    return sbext->getURI(mLevel,mVersion,mPackageVersion);
  }
#else
  ;
#endif //SWIG
  

  /**
   * Get the SBML Package Version of this SBMLExtensionNamespaces object.
   *
   * @return the SBML Package Version of this SBMLExtensionNamespaces object.
   */
  unsigned int getPackageVersion() const
#ifndef SWIG
  {
    return mPackageVersion;
  }
#else
  ;
#endif //SWIG

	/**
	 * Returns the name of the main package for this namespace.
	 *
	 * @return the name of the main package for this namespace.
	 * 
	 */
	virtual const std::string& getPackageName() const
#ifndef SWIG
	{
		return mPackageName;
	}
#else
	;
#endif //SWIG	
	
#ifndef SWIG
  /** @cond doxygen-libsbml-internal */

  void setPackageVersion(unsigned int pkgVersion)
  {
    mPackageVersion = pkgVersion;
  }

  /** @endcond doxygen-libsbml-internal */
#endif //SWIG


protected:  
  /** @cond doxygen-libsbml-internal */

  unsigned int mPackageVersion;
  const std::string& mPackageName;

  /** @endcond doxygen-libsbml-internal */
};

LIBSBML_CPP_NAMESPACE_END

#endif  /* __cplusplus */

#endif  /* SBMLExtensionNamespaces_h */
