/**
 * @file    GroupsModelPlugin.h
 * @brief   Definition of GroupsModelPlugin, the plugin class of
 *          groups package for the Model element.
 * @author  Akiya Jouraku
 *
 * $Id: GroupsModelPlugin.h 10673 2010-01-17 07:18:20Z ajouraku $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/branches/libsbml-5/src/packages/groups/extension/GroupsModelPlugin.h $
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

#ifndef GroupsModelPlugin_h
#define GroupsModelPlugin_h


#include <sbml/common/extern.h>
#include <sbml/common/sbmlfwd.h>
#include <sbml/groups/groupsfwd.h>
#include <sbml/SBMLTypeCodes.h>

#ifdef __cplusplus

#include <sbml/SBMLErrorLog.h>
#include <sbml/Model.h>
#include <sbml/xml/XMLInputStream.h>
#include <sbml/xml/XMLOutputStream.h>
#include <sbml/extension/SBasePlugin.h>
#include <sbml/groups/Group.h>

LIBSBML_CPP_NAMESPACE_BEGIN

class LIBSBML_EXTERN GroupsModelPlugin : public SBasePlugin
{
public:

  /**
   * Constructor
   */
  GroupsModelPlugin (const std::string &uri, const std::string &prefix,
                    GroupsPkgNamespaces *groupsns);


  /**
   * Copy constructor. Creates a copy of this SBase object.
   */
  GroupsModelPlugin(const GroupsModelPlugin& orig);


  /**
   * Destroy this object.
   */
  virtual ~GroupsModelPlugin ();


  /**
   * Assignment operator for GroupsModelPlugin.
   */
  GroupsModelPlugin& operator=(const GroupsModelPlugin& orig);


  /**
   * Creates and returns a deep copy of this GroupsModelPlugin object.
   * 
   * @return a (deep) copy of this SBase object
   */
  virtual GroupsModelPlugin* clone () const;


  // --------------------------------------------------------
  //
  // overridden virtual functions for reading/writing/checking 
  // elements
  //
  // --------------------------------------------------------

  /** @cond doxygen-libsbml-internal */

  /**
   * Subclasses must override this method to create, store, and then
   * return an SBML object corresponding to the next XMLToken in the
   * XMLInputStream if they have their specific elements.
   *
   * @return the SBML object corresponding to next XMLToken in the
   * XMLInputStream or NULL if the token was not recognized.
   */
  virtual SBase* createObject (XMLInputStream& stream);


  /**
   * Subclasses must override this method to write out their contained
   * SBML objects as XML elements if they have their specific elements.
   */
  virtual void writeElements (XMLOutputStream& stream) const;


  /**
   * Checks if this plugin object has all the required elements.
   *
   * Subclasses should override this function if they have their specific
   * elements.
   *
   * @return true if this pugin object has all the required elements,
   * otherwise false will be returned.
   */
  virtual bool hasRequiredElements() const ;


  /** ------------------------------------------------------------------
   *
   *  Additional public functions
   *
   * ------------------------------------------------------------------
   */
  
  /**
   * Returns the ListOfGroups in this plugin object.
   *
   * @return ListOfGroups object in this plugin object.
   */
  const ListOfGroups* getListOfGroups () const;


  /**
   * Returns the ListOfGroups in this plugin object.
   *
   * @return ListOfGroups object in this plugin object.
   */
  ListOfGroups* getListOfGroups ();


  /**
   * Returns the Group object that belongs to the given index. If the
   * index is invalid, NULL is returned.
   *
   * @param n the index number of the Group to get.
   *
   * @return the nth Group in the ListOfGroups.
   */
  const Group* getGroup (unsigned int n) const;


  /**
   * Returns the Group object that belongs to the given index. If the
   * index is invalid, NULL is returned.
   *
   * @param n the index number of the Group to get.
   *
   * @return the nth Group in the ListOfGroups.
   */
  Group* getGroup (unsigned int n);


  /**
   * Returns the group object based on its identifier.
   *
   * @param sid a string representing the identifier 
   * of the Group to get.
   * 
   * @return Group in the ListOfGroups with the given id
   * or NULL if no such Group exists.
   *
   * @see get(unsigned int n)
   * @see size()
   */
  Group* getGroup (const std::string& sid);


  /**
   * Returns the group object based on its identifier.
   *
   * @param sid a string representing the identifier 
   * of the Group to get.
   * 
   * @return Group in the ListOfGroups with the given id 
   * or NULL if no such Group exists.
   *
   * @see get(unsigned int n)
   * @see size()
   */
  const Group* getGroup (const std::string& sid) const;

  /**
   * Adds a copy of the given Group object to the list of groups.
   *
   * @param group the Group object to be added to the list of groups.
   *
   * @return integer value indicating success/failure of the
   * function.  @if clike The value is drawn from the
   * enumeration #OperationReturnValues_t. @endif The possible values
   * returned by this function are:
   * @li LIBSBML_OPERATION_SUCCESS
   */ 
  int addGroup (const Group* group);


  /**
   * Creates a new groups object and adds it to the list of groups objects
   * and returns it.
   *
   * @return a newly created Group object
   */
  Group* createGroup();


  /**
   * Removes the nth Group object from this plugin object and
   * returns a pointer to it.
   *
   * The caller owns the returned object and is responsible for
   *  deleting it.
   *
   * @param n the index of the Group object to remove
   *
   * @return the Group object removed.  As mentioned above, the 
   * caller owns the returned object. NULL is returned if the 
   * given index is out of range.
   */
  Group* removeGroup (unsigned int n);


  /**
   * Removes the Group object with the given id attribute from 
   * this plugin object and returns a pointer to it.
   *
   * The caller owns the returned object and is responsible for
   * deleting it.
   *
   * @param sid the id attribute of the Group object to remove
   *
   * @return the Group object removed.  As mentioned above, the 
   * caller owns the returned object. NULL is returned if the 
   * given index is out of range.
   */
  Group* removeGroup (const std::string& sid);


  /**
   * Returns the number of Group object in this plugin object.
   *
   * @return the number of Group object in this plugin object.
   */
  int getNumGroups() const;

  // ---------------------------------------------------------
  //
  // virtual functions (internal implementation) which should
  // be overridden by subclasses.
  //
  // ---------------------------------------------------------

  /** @cond doxygen-libsbml-internal */

  /**
   * Sets the parent SBMLDocument of this plugin object.
   *
   * Subclasses which contain one or more SBase derived elements must
   * override this function.
   *
   * @param d the SBMLDocument object to use
   *
   * @see connectToParent
   * @see enablePackageInternal
   */
  virtual void setSBMLDocument (SBMLDocument* d);


  /**
   * Sets the parent SBML object of this plugin object to
   * this object and child elements (if any).
   * (Creates a child-parent relationship by this plugin object)
   *
   * This function is called when this object is created by
   * the parent element.
   * Subclasses must override this this function if they have one
   * or more child elements.Also, SBasePlugin::connectToParent()
   * must be called in the overridden function.
   *
   * @param sbase the SBase object to use
   *
   * @see setSBMLDocument
   * @see enablePackageInternal
   */
  virtual void connectToParent (SBase *sbase);


  /**
   * Enables/Disables the given package with child elements in this plugin
   * object (if any).
   * (This is an internal implementation invoked from
   *  SBase::enablePakcageInternal() function)
   *
   * @note Subclasses in which one or more SBase derived elements are
   * defined must override this function.
   *
   * @see setSBMLDocument
   * @see connectToParent
   */
  virtual void enablePackageInternal(const std::string& pkgURI,
                                     const std::string& pkgPrefix, bool flag);
  /** @endcond doxygen-libsbml-internal */

protected:
  /** @cond doxygen-libsbml-internal */

  /*-- data members --*/

  ListOfGroups mGroups;

  /** @endcond doxygen-libsbml-internal */
};

LIBSBML_CPP_NAMESPACE_END

#endif  /* __cplusplus */
#endif  /* GroupsModelPlugin_h */
