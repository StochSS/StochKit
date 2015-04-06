/**
 * @file    Group.h
 * @brief   Definition of Group, the SBase derived class of groups package.
 * @author  Akiya Jouraku
 *
 * $Id: Group.h 12233 2010-11-26 22:44:21Z fbergmann $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/branches/libsbml-5/src/packages/groups/sbml/Group.h $
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


#ifndef Group_H__
#define Group_H__


#include <sbml/common/extern.h>
#include <sbml/common/sbmlfwd.h>
#include <sbml/groups/groupsfwd.h>

#ifdef __cplusplus

#include <string>

#include <sbml/SBase.h>
#include <sbml/ListOf.h>
#include <sbml/groups/GroupsExtension.h>
#include <sbml/groups/Member.h>

LIBSBML_CPP_NAMESPACE_BEGIN


class LIBSBML_EXTERN Group : public SBase
{
protected:

  std::string   mId;
  std::string   mName;
  ListOfMembers mMembers;

public:

  /**
   * Creates a new Group with the given level, version, and package version.
   */
   Group(unsigned int level      = GroupsExtension::getDefaultLevel(),
         unsigned int version    = GroupsExtension::getDefaultVersion(),
         unsigned int pkgVersion = GroupsExtension::getDefaultPackageVersion());


  /**
   * Creates a new Group with the given GroupsPkgNamespaces object.
   */
   Group(GroupsPkgNamespaces* groupsns);


  /**
   * Copy constructor.
   */
   Group(const Group& source);

  /**
   * Assignment operator.
   */
   Group& operator=(const Group& source);


  /**
   * Destructor.
   */ 
  virtual ~Group ();


  /**
   * Returns the value of the "id" attribute of this Group.
   *
   * @return the value of the "id" attribute of this Group.
   */
  virtual const std::string& getId () const;


  /**
   * Predicate returning @c true or @c false depending on whether this
   * Group's "id" attribute has been set.
   *
   * @return @c true if this Group's "id" attribute has been set, 
   * otherwise @c false is returned.
   */
  virtual bool isSetId () const;

  
  /**
   * Sets the value of the "id" attribute of this Group.
   *
   * @return integer value indicating success/failure of the
   * function.  @if clike The value is drawn from the
   * enumeration #OperationReturnValues_t. @endif The possible values
   * returned by this function are:
   * @li LIBSBML_OPERATION_SUCCESS
   * @li LIBSBML_INVALID_ATTRIBUTE_VALUE
   */
  virtual int setId (const std::string& id);


  /**
   * Unsets the value of the "id" attribute of this Group.
   *
   * @return integer value indicating success/failure of the
   * function.  @if clike The value is drawn from the
   * enumeration #OperationReturnValues_t. @endif The possible values
   * returned by this function are:
   * @li LIBSBML_OPERATION_SUCCESS
   * @li LIBSBML_OPERATION_FAILED
   */
  virtual int unsetId ();


  /**
   * Sets the value of the "name" attribute of this Group.
   *
   * The string in @p name is copied.
   *
   * @htmlinclude libsbml-comment-set-methods.html
   *
   * @param name the new name for the Group
   *
   * @return integer value indicating success/failure of the
   * function.  @if clike The value is drawn from the
   * enumeration #OperationReturnValues_t. @endif The possible values
   * returned by this function are:
   * @li LIBSBML_OPERATION_SUCCESS
   * @li LIBSBML_INVALID_ATTRIBUTE_VALUE
   */
  virtual int setName (const std::string& name);


  /**
   * Returns the value of the "name" attribute of this Group.
   * 
   * @return the name of this Group.
   */
  virtual const std::string& getName () const;


  /**
   * Predicate returning @c true or @c false depending on whether this
   * Member's "name" attribute has been set.
   *
   * @htmlinclude libsbml-comment-set-methods.html
   * 
   * @return @c true if the "name" attribute of this Member has been
   * set, @c false otherwise.
   */
  virtual bool isSetName () const;


  /**
   * Unsets the value of the "name" attribute of this Group.
   *
   * @htmlinclude libsbml-comment-set-methods.html
   *
   * @return integer value indicating success/failure of the
   * function.  @if clike The value is drawn from the
   * enumeration #OperationReturnValues_t. @endif The possible values
   * returned by this function are:
   * @li LIBSBML_OPERATION_SUCCESS
   * @li LIBSBML_OPERATION_FAILED
   */
  virtual int unsetName ();


  /**
   * Returns the ListOf object that holds all members.
   *
   * @return the ListOf object that holds all members.
   */ 
  const ListOfMembers* getListOfMembers () const;

  /**
   * Returns the member with the given index.
   * If the index is invalid, NULL is returned.
   *
   * @param n the index number of the Member to get.
   *
   * @return the nth Member in the ListOfMembers.
   */ 
  Member* getMember (unsigned int n);

  /**
   * Returns the member with the given index.
   * If the index is invalid, NULL is returned.
   *
   * @param n the index number of the Member to get.
   *
   * @return the nth Member in the ListOfMembers.
   */ 
  const Member* getMember (unsigned int n) const;

  /**
   * Returns the member with the given symbol.
   * If the index is invalid, NULL is returned.
   *
   * @param symbol a string representing the symbol attribute
   * of the Member to get.
   * 
   * @return Member in the ListOfMembers with the given symbol
   * or NULL if no such Member exists.
   */
  Member* getMember (const std::string& symbol);


  /**
   * Returns the member with the given symbol.
   * If the index is invalid, NULL is returned.
   *
   * @param symbol a string representing the symbol attribute
   * of the Member to get.
   * 
   * @return Member in the ListOfMembers with the given symbol
   * or NULL if no such Member exists.
   */
  const Member* getMember (const std::string& symbol) const;


  /**
   * Adds a copy of the given Member objcect to the list of members.
   *
   * @param member the Member object to be added to the list of 
   * members.
   *
   * @return integer value indicating success/failure of the
   * function.  @if clike The value is drawn from the
   * enumeration #OperationReturnValues_t. @endif The possible values
   * returned by this function are:
   * @li LIBSBML_OPERATION_SUCCESS
   */
  int addMember (const Member* member);


  /**
   * Returns the number of members for this group.
   *
   * @return the number of members for this group.
   */
  unsigned int getNumMembers () const;


  /**
   * Creates a Member object, adds it to the end of the
   * member objects list and returns a pointer to the newly
   * created object.
   *
   * @return a newly created Member object
   */
  Member* createMember ();


  /**
   * Removes the member with the given index from the group.
   * A pointer to the member that was removed is returned.
   * If no member has been removed, NULL is returned.
   *
   * @param n the index of the Member object to remove
   *
   * @return the Member object removed.  As mentioned above, 
   * the caller owns the returned object. NULL is returned if 
   * the given index is out of range.
   */
  Member* removeMember(unsigned int index);


  /**
   * Removes the member with the given symbol from the group.
   * A pointer to the member that was removed is returned.
   * If no member has been removed, NULL is returned.
   *
   * @param symbol the symbol attribute of the Member object to remove
   *
   * @return the Member object removed.  As mentioned above, 
   * the caller owns the returned object. NULL is returned if 
   * the given index is out of range.
   */
  Member* removeMember(const std::string& symbol);


  /**
   * Subclasses should override this method to return XML element name of
   * this SBML object.
   *
   * @return the string of the name of this element.
   */
  virtual const std::string& getElementName () const ;


  /**
   * @return a (deep) copy of this Model.
   */
  virtual SBase* clone () const;


  /**
   * @return the typecode (int) of this SBML object or SBML_UNKNOWN
   * (default).
   *
   * @see getElementName()
   */
  int getTypeCode () const;


  /** @cond doxygen-libsbml-internal */
  /**
   * Subclasses should override this method to write out their contained
   * SBML objects as XML elements.  Be sure to call your parents
   * implementation of this method as well.  For example:
   *
   *   SBase::writeElements(stream);
   *   mReactans.write(stream);
   *   mProducts.write(stream);
   *   ...
   */
  virtual void writeElements (XMLOutputStream& stream) const;


  /**
   * Accepts the given SBMLVisitor.
   *
   * @return the result of calling <code>v.visit()</code>, which indicates
   * whether or not the Visitor would like to visit the SBML object's next
   * sibling object (if available).
   */
  virtual bool accept (SBMLVisitor& v) const;
  /** @endcond doxygen-libsbml-internal */


  /** @cond doxygen-libsbml-internal */
  /**
   * Sets the parent SBMLDocument of this SBML object.
   *
   * @param d the SBMLDocument object to use
   */
  virtual void setSBMLDocument (SBMLDocument* d);


  /**
   * Sets this SBML object to child SBML objects (if any).
   * (Creates a child-parent relationship by the parent)
   *
   * Subclasses must override this function if they define
   * one ore more child elements.
   * Basically, this function needs to be called in
   * constructor, copy constructor, assignment operator.
   *
   * @see setSBMLDocument
   * @see enablePackageInternal
   */
  virtual void connectToChild ();


  /**
   * Enables/Disables the given package with this element and child
   * elements (if any).
   * (This is an internal implementation for enablePakcage function)
   *
   * @note Subclasses in which one or more child elements are defined
   * must override this function.
   */
  virtual void enablePackageInternal(const std::string& pkgURI,
                                     const std::string& pkgPrefix, bool flag);
  /** @endcond doxygen-libsbml-internal */


  /** @cond doxygen-libsbml-internal */
  /* function returns true if component has all the required
   * elements
   * needs to be overloaded for each component
   */
  virtual bool hasRequiredElements() const ;
  /** @endcond doxygen-libsbml-internal */

    
protected:
  /**
   * @return the SBML object corresponding to next XMLToken in the
   * XMLInputStream or NULL if the token was not recognized.
   */
  virtual SBase*
  createObject (XMLInputStream& stream);

  /**
   * Subclasses should override this method to get the list of
   * expected attributes.
   * This function is invoked from corresponding readAttributes()
   * function.
   */
  virtual void addExpectedAttributes(ExpectedAttributes& attributes);


  /**
   * Subclasses should override this method to read values from the given
   * XMLAttributes set into their specific fields.  Be sure to call your
   * parents implementation of this method as well.
   */
  virtual void readAttributes (const XMLAttributes& attributes, 
                               const ExpectedAttributes& expectedAttributes);

  /**
   * Subclasses should override this method to write their XML attributes
   * to the XMLOutputStream.  Be sure to call your parents implementation
   * of this method as well.  For example:
   *
   *   SBase::writeAttributes(stream);
   *   stream.writeAttribute( "id"  , mId   );
   *   stream.writeAttribute( "name", mName );
   *   ...
   */
  virtual void writeAttributes (XMLOutputStream& stream) const;

};

class LIBSBML_EXTERN ListOfGroups : public ListOf
{
public:

  /**
   * @return a (deep) copy of this ListOfUnitDefinitions.
   */
  virtual SBase* clone () const;


  /**
   * Creates a new ListOfGroups with the given level, version, and package version.
   */
   ListOfGroups(unsigned int level      = GroupsExtension::getDefaultLevel(), 
                unsigned int version    = GroupsExtension::getDefaultVersion(), 
                unsigned int pkgVersion = GroupsExtension::getDefaultPackageVersion());


  /**
   * Creates a new ListOfGroups with the given GroupsPkgNamespaces object.
   */
   ListOfGroups(GroupsPkgNamespaces* groupsns);


  /**
   * Get a Group from the ListOfGroups.
   *
   * @param n the index number of the Group to get.
   * 
   * @return the nth Group in this ListOfGroups.
   *
   * @see size()
   */
  virtual Group* get(unsigned int n); 


  /**
   * Get a Group from the ListOfGroups.
   *
   * @param n the index number of the Group to get.
   * 
   * @return the nth Group in this ListOfGroups.
   *
   * @see size()
   */
  virtual const Group * get(unsigned int n) const; 


  /**
   * Get a Group from the ListOfGroups
   * based on its identifier.
   *
   * @param sid a string representing the identifier 
   * of the Group to get.
   * 
   * @return Group in this ListOfGroups
   * with the given id or NULL if no such
   * Group exists.
   *
   * @see get(unsigned int n)
   * @see size()
   */
  virtual Group* get (const std::string& sid);


  /**
   * Get a Group from the ListOfGroups
   * based on its identifier.
   *
   * @param sid a string representing the identifier 
   * of the Group to get.
   * 
   * @return Group in this ListOfGroups
   * with the given id or NULL if no such
   * Group exists.
   *
   * @see get(unsigned int n)
   * @see size()
   */
  virtual const Group* get (const std::string& sid) const;


  /**
   * Removes the nth item from this ListOfGroups items and returns a pointer to
   * it.
   *
   * The caller owns the returned item and is responsible for deleting it.
   *
   * @param n the index of the item to remove
   *
   * @see size()
   */
  virtual Group* remove (unsigned int n);


  /**
   * Removes item in this ListOfGroups items with the given identifier.
   *
   * The caller owns the returned item and is responsible for deleting it.
   * If none of the items in this list have the identifier @p sid, then @c
   * NULL is returned.
   *
   * @param sid the identifier of the item to remove
   *
   * @return the item removed.  As mentioned above, the caller owns the
   * returned item.
   */
  virtual Group* remove (const std::string& sid);


  /**
   * @return the typecode (int) of SBML objects contained in this ListOf or
   * SBML_UNKNOWN (default).
   */
  virtual int getItemTypeCode () const;

  /**
   * Subclasses should override this method to return XML element name of
   * this SBML object.
   *
   * @return the string of the name of this element.
   */
  virtual const std::string& getElementName () const;


protected:

  /**
   * @return the SBML object corresponding to next XMLToken in the
   * XMLInputStream or NULL if the token was not recognized.
   */
  virtual SBase* createObject (XMLInputStream& stream);

  virtual void writeXMLNS (XMLOutputStream& stream) const;
};


LIBSBML_CPP_NAMESPACE_END

#endif /* __cplusplus */


#ifndef SWIG

LIBSBML_CPP_NAMESPACE_BEGIN
BEGIN_C_DECLS

//
// C API will be added here.
//

END_C_DECLS
LIBSBML_CPP_NAMESPACE_END


#endif  /* !SWIG */
#endif  /* Group_H__ */
