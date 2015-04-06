/*
 * @file    Member.h
 * @brief   Definition of Member, the SBase derived class of groups package.
 * @author  Akiya Jouraku
 *
 * $Id: Member.h 12233 2010-11-26 22:44:21Z fbergmann $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/branches/libsbml-5/src/packages/groups/sbml/Member.h $
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


#ifndef Member_H__
#define Member_H__

#include <sbml/common/extern.h>
#include <sbml/common/sbmlfwd.h>
#include <sbml/groups/groupsfwd.h>

#ifdef __cplusplus

#include <string>

#include <sbml/SBase.h>
#include <sbml/ListOf.h>
#include <sbml/groups/GroupsExtension.h>

LIBSBML_CPP_NAMESPACE_BEGIN


class LIBSBML_EXTERN Member : public SBase
{
protected:

  std::string mSymbol;

public:

  /**
   * Creates a new Member with the given level, version, and package version.
   */
   Member(unsigned int level      = GroupsExtension::getDefaultLevel(),
          unsigned int version    = GroupsExtension::getDefaultVersion(),
          unsigned int pkgVersion = GroupsExtension::getDefaultPackageVersion());


  /**
   * Creates a new Member with the given GroupsPkgNamespaces object.
   */
   Member(GroupsPkgNamespaces* groupsns);


  /**
   * Copy constructor.
   */
   Member(const Member& source);


  /**
   * Assignment operator.
   */
   Member& operator=(const Member& source);


  /**
   * Destructor.
   */ 
  virtual ~Member ();


  /**
   * Returns the string of the "symbol" attribute of this Member.
   *
   * @return the string of the "symbol" attribute of this Member.
   */
  virtual const std::string& getSymbol () const;


  /**
   * Predicate returning @c true or @c false depending on whether this
   * Member's "symbol" attribute has been set.
   *
   * @return @c true if this Member's "symbol" attribute has been set, 
   * otherwise @c false is returned.
   */
  virtual bool isSetSymbol () const;

  
  /**
   * Sets the SIdRef string of the "symbol" attribute of this Member.
   *
   * @param symbol a SIdRef string to be set.
   *
   * @return integer value indicating success/failure of the
   * function.  @if clike The value is drawn from the
   * enumeration #OperationReturnValues_t. @endif The possible values
   * returned by this function are:
   * @li LIBSBML_OPERATION_SUCCESS
   * @li LIBSBML_INVALID_ATTRIBUTE_VALUE
   */
  virtual int setSymbol (const std::string& symbol);


  /**
   * Unsets the value of the "id" attribute of this Member.
   *
   * @return integer value indicating success/failure of the
   * function.  @if clike The value is drawn from the
   * enumeration #OperationReturnValues_t. @endif The possible values
   * returned by this function are:
   * @li LIBSBML_OPERATION_SUCCESS
   * @li LIBSBML_OPERATION_FAILED
   */
  virtual int unsetSymbol ();


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

class LIBSBML_EXTERN ListOfMembers : public ListOf
{
public:

  /**
   * @return a (deep) copy of this ListOfUnitDefinitions.
   */
  virtual SBase* clone () const;


  /**
   * Creates a new ListOfMembers with the given level, version, and package version.
   */
   ListOfMembers(unsigned int level      = GroupsExtension::getDefaultLevel(), 
                 unsigned int version    = GroupsExtension::getDefaultVersion(), 
                 unsigned int pkgVersion = GroupsExtension::getDefaultPackageVersion());


  /**
   * Creates a new ListOfMembers with the given GroupsPkgNamespaces object.
   */
   ListOfMembers(GroupsPkgNamespaces* groupssns);


  /**
   * Get a Member from the ListOfMembers.
   *
   * @param n the index number of the Member to get.
   * 
   * @return the nth Member in this ListOfMembers.
   *
   * @see size()
   */
  virtual Member * get(unsigned int n); 


  /**
   * Get a Member from the ListOfMembers.
   *
   * @param n the index number of the Member to get.
   * 
   * @return the nth Member in this ListOfMembers.
   *
   * @see size()
   */
  virtual const Member * get(unsigned int n) const; 

  /**
   * Get a Member from the ListOfMembers
   * based on its identifier.
   *
   * @param sid a string representing the identifier 
   * of the Member to get.
   * 
   * @return Member in this ListOfMembers
   * with the given id or NULL if no such
   * Member exists.
   *
   * @see get(unsigned int n)
   * @see size()
   */
  virtual Member* get (const std::string& sid);


  /**
   * Get a Member from the ListOfMembers
   * based on its identifier.
   *
   * @param sid a string representing the identifier 
   * of the Member to get.
   * 
   * @return Member in this ListOfMembers
   * with the given id or NULL if no such
   * Member exists.
   *
   * @see get(unsigned int n)
   * @see size()
   */
  virtual const Member* get (const std::string& sid) const;


  /**
   * Removes the nth item from this ListOfMembers items and returns a pointer to
   * it.
   *
   * The caller owns the returned item and is responsible for deleting it.
   *
   * @param n the index of the item to remove
   * @return the item removed.  As mentioned above, the caller owns the
   * returned item.
   *
   * @see size()
   */
  virtual Member* remove (unsigned int n);


  /**
   * Removes item in this ListOfMembers items with the given identifier.
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
  virtual Member* remove (const std::string& sid);


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
};

/** @cond doxygen-libsbml-internal */
/**
 * Used by ListOfMembers::get() to lookup an SBase based by its 
 * symbol
 */
#ifndef SWIG
template<>
struct IdEq<Member> : public std::unary_function<SBase*, bool>
{
  const std::string& id;

  IdEq (const std::string& id) : id(id) { }
  bool operator() (SBase* sb) 
       { return static_cast <Member*> (sb)->getSymbol() == id; }
};
#endif
/** @endcond doxygen-libsbml-internal */

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
#endif  /* Member_H__ */
