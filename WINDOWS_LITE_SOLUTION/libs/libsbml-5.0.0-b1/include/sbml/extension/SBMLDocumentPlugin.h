/**
 * @file    SBMLDocumentPlugin.h
 * @brief   Definition of SBMLDocumentPlugin, the derived class of
 *          SBasePlugin.
 * @author  Akiya Jouraku
 *
 * $Id: SBMLDocumentPlugin.h 10667 2010-01-16 10:20:44Z ajouraku $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/branches/libsbml-5/src/extension/SBMLDocumentPlugin.h $
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

#ifndef SBMLDocumentPlugin_h
#define SBMLDocumentPlugin_h

#include <sbml/common/sbmlfwd.h>
#include <sbml/SBMLTypeCodes.h>
#include <sbml/SBMLErrorLog.h>
#include <sbml/SBMLDocument.h>
#include <sbml/xml/XMLInputStream.h>
#include <sbml/xml/XMLOutputStream.h>
#include <sbml/extension/SBasePlugin.h>

#ifdef __cplusplus

LIBSBML_CPP_NAMESPACE_BEGIN

//
// (NOTE) Plugin objects for the SBMLDocument element must be this class or 
//        a derived class of this class.
//        Package developers should use this class as-is if only "required" 
//        attribute is added in the SBMLDocument element by their packages, 
//        otherwise developers must implement a derived class of this class 
//        and use the class as the plugin object for the SBMLDocument element. 
//

class LIBSBML_EXTERN SBMLDocumentPlugin : public SBasePlugin
{
public:

  /**
   *  Constructor
   *
   * @param uri the URI of package 
   * @param prefix the prefix for the given package
   */
  SBMLDocumentPlugin (const std::string &uri, const std::string &prefix,
                      SBMLNamespaces *sbmlns);


  /**
   * Copy constructor. Creates a copy of this object.
   */
  SBMLDocumentPlugin(const SBMLDocumentPlugin& orig);


  /**
   * Destroy this object.
   */
  virtual ~SBMLDocumentPlugin ();

  /**
   * Assignment operator for SBMLDocumentPlugin.
   */
  SBMLDocumentPlugin& operator=(const SBMLDocumentPlugin& orig);


  /**
   * Creates and returns a deep copy of this SBMLDocumentPlugin object.
   * 
   * @return a (deep) copy of this object
   */
  virtual SBMLDocumentPlugin* clone () const;


  // ----------------------------------------------------------
  //
  // overridden virtual functions for reading/writing/checking
  // attributes
  //
  // ----------------------------------------------------------

#ifndef SWIG

  /** @cond doxygen-libsbml-internal */

  /**
   * Subclasses should override this method to get the list of
   * expected attributes.
   * This function is invoked from corresponding readAttributes()
   * function.
   */
  virtual void addExpectedAttributes(ExpectedAttributes& attributes);


  /**
   * Reads the attributes of corresponding package in SBMLDocument element.
   */
  virtual void readAttributes (const XMLAttributes& attributes,
                               const ExpectedAttributes& expectedAttributes);


  /**
   * Writes the attributes of corresponding package in SBMLDocument element.
   */
  virtual void writeAttributes (XMLOutputStream& stream) const;

  /** @endcond doxygen-libsbml-internal */

#endif //SWIG

  // -----------------------------------------------------------
  //
  // Additional public functions for manipulating attributes of 
  // corresponding package in SBMLDocument element.
  //
  // -----------------------------------------------------------


  /**
   *
   * Returns the bool value of "required" attribute of corresponding 
   * package in SBMLDocument element.
   *
   * @return the bool value of "required" attribute of corresponding
   * package in SBMLDocument element.
   */
  virtual bool getRequired() const;
  

  /**
   *
   * Sets the bool value of "required" attribute of corresponding package
   * in SBMLDocument element.
   *
   * @param value the bool value of "required" attribute of corresponding 
   * package in SBMLDocument element.
   *
   * @return 
   */
  virtual int setRequired(bool value);

protected:
  /** @cond doxygen-libsbml-internal */

  /*-- data members --*/

  bool mRequired;

  /** @endcond doxygen-libsbml-internal */

};

LIBSBML_CPP_NAMESPACE_END

#endif  /* __cplusplus */
#endif  /* SBMLDocumentPlugin_h */
